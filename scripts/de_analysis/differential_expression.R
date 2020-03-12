# A pipeline for the differential expression analysis of a dual RNA-seq dataset of chickens and Eimeria tenella parsites

library(edgeR)
library(EnhancedVolcano)
library(org.Gg.eg.db)
library(ggbiplot)
library(KEGGREST)

# Specify the path to the read files and metadata and read them into R
path_to_chicken_reads <- "results/htseq/reverse/in_vitro_fusion/processed_reads/chicken_counts.csv"
path_to_eimeria_reads <- "results/htseq/reverse/in_vitro_fusion/processed_reads/eimeria_counts.csv"
path_to_annotation_file <- "results/htseq/reverse/in_vitro_fusion/processed_reads/metadata_table.csv"

chicken_reads <- read.delim(path_to_chicken_reads, sep = ",", check.names=FALSE, stringsAsFactors=FALSE)
eimeria_reads <- read.delim(path_to_eimeria_reads, sep = ",", check.names=FALSE, stringsAsFactors=FALSE)
annotation <- read.delim(path_to_annotation_file, sep = ",", check.names=FALSE, stringsAsFactors=FALSE)

# Remove uninfected samples from the Eimeria data and the sample with only Eimeria data from the chicken data
infected_samples <- rep("Infected", 39) == annotation$Infection_status
eimeria_annotation <- annotation[infected_samples,]
eimeria_reads <- eimeria_reads[,c("gene_id", eimeria_annotation$File_name)]
rownames(eimeria_annotation) <- eimeria_annotation$File_name

eimeria_sample <- rep("33_S29", 39) == annotation$File_name
chicken_annotation <- annotation[!eimeria_sample,]
chicken_reads$`33_S29` <- NULL
rownames(chicken_annotation) <- chicken_annotation$File_name

# Define the grouping for the chicken and Eimeria samples.  For the chicken, each timepoint and infection status 
# is its own group while only time is a factor for the Eimeria samples
group_chicken <- factor(paste(substr(chicken_annotation[colnames(chicken_reads)[-1],"Infection_status"],1,1),
                              chicken_annotation[colnames(chicken_reads)[-1],"Timepoint_hours"],
                              sep="."))
group_eimeria <- factor(eimeria_annotation[colnames(eimeria_reads)[-1],"Timepoint_hours"])

# Run a PCA on the data to identify and remove outliers and other potential issues
pca_chicken <- prcomp(chicken_reads[,-1])
pca_eimeria <- prcomp(eimeria_reads[,-1])

ggbiplot(pca_chicken, labels = chicken_reads[,1]) + coord_cartesian(xlim = c(0, 160))
ggbiplot(pca_eimeria, labels = eimeria_reads[,1]) + coord_cartesian(xlim = c(0, 60))

chicken_remove = which(chicken_reads[,1] == "LOC112533599") # Talk about why it's removed in the report, find other rRNA genes and remove them

chicken_reads = chicken_reads[-chicken_remove,]

pca_chicken <- prcomp(chicken_reads[,-1])

ggbiplot(pca_chicken, labels = chicken_reads[,1])# + coord_cartesian(xlim = c(0, 160))

# Create the DGElist object containing the count data in a format that edgeR can work with
chicken_dgelist <- DGEList(counts = chicken_reads[,2:length(chicken_reads)], genes = chicken_reads[,1], group = group_chicken)
eimeria_dgelist <- DGEList(counts = eimeria_reads[,2:length(eimeria_reads)], genes = eimeria_reads[,1], group = group_eimeria)

rownames(chicken_dgelist$counts) <- rownames(chicken_dgelist$genes) <- chicken_dgelist$genes[,1]
rownames(eimeria_dgelist$counts) <- rownames(eimeria_dgelist$genes) <- eimeria_dgelist$genes[,1]

# Filter out the genes that lacked annotation
gene_list <- chicken_dgelist$genes[,1]
idfound <- gene_list %in% mappedRkeys(org.Gg.egSYMBOL)
chicken_dgelist <- chicken_dgelist[idfound,]

# Filter out genes that are lowly expressed in all samples
keep_chicken <- filterByExpr(chicken_dgelist)
keep_eimeria <- filterByExpr(eimeria_dgelist)
chicken_dgelist_filt <- chicken_dgelist[keep_chicken, ,] 
eimeria_dgelist_filt <- eimeria_dgelist[keep_eimeria, ,]

# Setting the library size to the right value after filtering
chicken_dgelist_filt$samples$lib.size <- colSums(chicken_dgelist_filt$counts)
eimeria_dgelist_filt$samples$lib.size <- colSums(eimeria_dgelist_filt$counts)

# Normalize gene counts to account for a few genes dominating the counts of the samples
chicken_dgelist_filt_norm <- calcNormFactors(chicken_dgelist_filt)
eimeria_dgelist_filt_norm <- calcNormFactors(eimeria_dgelist_filt)

# Check MDS plots for outliers or other potential issues
plotMDS(chicken_dgelist_filt_norm)
plotMDS(eimeria_dgelist_filt_norm)

# The chicken analysis design
chicken_dgelist_filt_norm$samples$group <- relevel(chicken_dgelist_filt_norm$samples$group, ref="U.0")
design_chicken <- model.matrix(~group, data = chicken_dgelist_filt_norm$samples)
#design_chicken <- model.matrix(~infection_status+timepoint_chicken)
rownames(design_chicken) <- colnames(chicken_dgelist_filt_norm)

# Define the design of the eimeria analysis
design_eimeria <- model.matrix(~group, data = eimeria_dgelist_filt_norm$samples)
rownames(design_eimeria) <- colnames(eimeria_dgelist_filt_norm)

# Estimate the dispersion for each dataset and plot for the biological coefficient of variation (BCV)
chicken_disp <- estimateDisp(chicken_dgelist_filt_norm, design_chicken, robust = TRUE)
eimeria_disp <- estimateDisp(eimeria_dgelist_filt_norm, design_eimeria, robust = TRUE)

plotBCV(chicken_disp)
plotBCV(eimeria_disp)

# Test for differential expression in chicken
fit_chicken <- glmQLFit(chicken_disp, design_chicken)
contrasts_chicken <- makeContrasts(UI.2=groupI.2-groupU.2, UI.4=groupI.4-groupU.4, UI.12=groupI.12-groupU.12, 
                                   UI.24=groupI.24-groupU.24, UI.48=groupI.48-groupU.48, UI.72=groupI.72-groupU.72,
                                   levels = design_chicken)

qlf_chicken_I_vs_U.0 <- glmQLFTest(fit_chicken, coef = 2:7)
qlf_chicken_U_vs_U.0 <- glmQLFTest(fit_chicken, coef = 8:13)
qlf_chicken_UI.2 <- glmQLFTest(fit_chicken, contrast = contrasts_chicken[,"UI.2"])
qlf_chicken_UI.4 <- glmQLFTest(fit_chicken, contrast = contrasts_chicken[,"UI.4"])
qlf_chicken_UI.12 <- glmQLFTest(fit_chicken, contrast = contrasts_chicken[,"UI.12"])
qlf_chicken_UI.24 <- glmQLFTest(fit_chicken, contrast = contrasts_chicken[,"UI.24"])
qlf_chicken_UI.48 <- glmQLFTest(fit_chicken, contrast = contrasts_chicken[,"UI.48"])
qlf_chicken_UI.72 <- glmQLFTest(fit_chicken, contrast = contrasts_chicken[,"UI.72"])

chicken_qlf_list <- list(qlf_chicken_UI.2,qlf_chicken_UI.4,qlf_chicken_UI.12,qlf_chicken_UI.24,
                      qlf_chicken_UI.48,qlf_chicken_UI.72)
topgenes_chicken_list <- lapply(chicken_qlf_list,
                                function(x) topTags(x, n = 12638))

# Test for differential expression in Eimeria
fit_eimeria <- glmQLFit(eimeria_disp, design_eimeria)

qlf_eimeria_2 <- glmQLFTest(fit_eimeria, coef = 3)
qlf_eimeria_4 <- glmQLFTest(fit_eimeria, coef = 5)
qlf_eimeria_12 <- glmQLFTest(fit_eimeria, coef = 2)
qlf_eimeria_24 <- glmQLFTest(fit_eimeria, coef = 4)
qlf_eimeria_48 <- glmQLFTest(fit_eimeria, coef = 6)
qlf_eimeria_72 <- glmQLFTest(fit_eimeria, coef = 7)
qlf_eimeria_all <- glmQLFTest(fit_eimeria, coef = 2:7)

eimeria_qlf_list <- list(qlf_eimeria_2,qlf_eimeria_4,qlf_eimeria_12,qlf_eimeria_24,qlf_eimeria_48,
                      qlf_eimeria_72)
topgenes_eimeria_list <- lapply(eimeria_qlf_list,
                                function(x) topTags(x, n = 1000))

# Volcano plots for all comparisons
timepoints <- c("2h","4h","12h","24h","48h","72h")
i = 1
while (i <= length(chicken_qlf_list)) {
  p <- EnhancedVolcano(chicken_qlf_list[[i]]$table,
                  lab = rownames(chicken_qlf_list[[i]]$table),
                  title = paste('Chicken in vitro', timepoints[i]),
                  x = 'logFC',
                  y = 'PValue')
  print(p)
  i = i + 1
}

i = 1
while (i <= length(eimeria_qlf_list)) {
  p <- EnhancedVolcano(eimeria_qlf_list[[i]]$table,
                  lab = rownames(eimeria_qlf_list[[i]]$table),
                  title = paste('Eimeria in vitro', timepoints[i]),
                  x = 'logFC',
                  y = 'PValue')
  print(p)
  i = i + 1
}

# Get the Entrez Gene identifiers for the chicken genes for GO and KEGG analysis
egSYMBOL <- toTable(org.Gg.egSYMBOL)

gene_list <- chicken_dgelist_filt_norm$genes[,1]

m <- match(gene_list, egSYMBOL$symbol)

entrez_geneIDs <- egSYMBOL$gene_id[m]

# GO category and KEGG pathway analysis of the data.  For each list, the data is in order of increasing time
go_chicken_list <- lapply(chicken_qlf_list, 
                          function(x) goana(x, geneid = entrez_geneIDs, FDR = 0.05, species = "Gg"))
kegg_chicken_list <- lapply(chicken_qlf_list, 
                            function(x) kegga(x, geneid = entrez_geneIDs, FDR = 0.05, species.KEGG = "gga"))
topgo_chicken_up_list <- lapply(go_chicken_list, 
                                function(x) topGO(x, ont="BP", sort="Up", n=30))
topgo_chicken_down_list <- lapply(go_chicken_list, 
                                  function(x) topGO(x, ont="BP", sort="Down", n=30))
topkegg_chicken_up_list <- lapply(kegg_chicken_list,
                                  function(x) topKEGG(x, sort="Up", n=30))
topkegg_chicken_down_list <- lapply(kegg_chicken_list,
                                    function(x) topKEGG(x, sort="Down", n=30))

# Get a lsit of genes that belong to each GO category of interest
get_de_go_cat_genes <- function(de_gene_list, go_id, database_go, database_symbols) {
  cat_genes <- unique(database_go$gene_id[database_go$go_id == go_id])
  cat_gene_ids <- database_symbols[match(cat_genes, database_symbols$gene_id),]$symbol
  cat_de_genes <- lapply(de_gene_list, function(x) x$table[cat_gene_ids,])
  cat_de_genes <- lapply(cat_de_genes, function(x) x[complete.cases(x),])
  return(cat_de_genes)
}

egGO <- toTable(org.Gg.egGO2ALLEGS)
positive_regulation_nitrogen_metabolism_gene_list <- get_de_go_cat_genes(topgenes_chicken_list, "GO:0051173", egGO, egSYMBOL)

# Get a list of genes that belong to each KEGG category of interest
get_de_kegg_cat_genes <- function(de_gene_list, organism_kegg_id, pathway_id, database_symbols) {
  # For a given organism, list of differentially expressed genes in a set of experiments and 
  # KEGG pathway of interest, this function returns a list containing information on the DE of
  # each gene in the pathway for this organism.
  pathway_genes <- keggLink("pathway", organism_kegg_id) # Get all genes connected to KEGG pathways for this organism
  cat_genes <- names(pathway_genes[pathway_genes == pathway_id]) # Get the names of genes connected to the pathway of interest
  cat_genes <- sapply(cat_genes, function(x) substr(x, 5, nchar(x)))
  cat_gene_ids <- database_symbols[match(cat_genes, database_symbols$gene_id),]$symbol
  cat_de_genes <- lapply(de_gene_list, function(x) x$table[cat_gene_ids,])
  cat_de_genes <- lapply(cat_de_genes, function(x) x[complete.cases(x),])
  return(cat_de_genes)
}

toll_like_gene_list <- get_de_kegg_cat_genes(topgenes_chicken_list,"gga", "path:gga04620", egSYMBOL)
























