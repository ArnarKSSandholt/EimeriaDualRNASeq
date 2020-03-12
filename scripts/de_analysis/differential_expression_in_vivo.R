# A pipeline for the differential expression analysis of a dual RNA-seq dataset of chickens and 
# Eimeria tenella parsites from in vivo samples

library(edgeR)
library(EnhancedVolcano)
library(org.Gg.eg.db)
library(ggplot2)
library(ggbiplot)
library(KEGGREST)

# Specify the path to the read files and metadata and read them into R
path_to_chicken_reads <- "results/htseq/reverse/SH-2259/processed_reads/chicken_counts.csv"
path_to_eimeria_reads <- "results/htseq/reverse/SH-2259/processed_reads/eimeria_counts.csv"
path_to_annotation_file <- "results/htseq/reverse/SH-2259/processed_reads/metadata_table.csv"

chicken_reads <- read.delim(path_to_chicken_reads, sep = ",", check.names=FALSE, stringsAsFactors=FALSE)
eimeria_reads <- read.delim(path_to_eimeria_reads, sep = ",", check.names=FALSE, stringsAsFactors=FALSE)
annotation <- read.delim(path_to_annotation_file, sep = ",", check.names=FALSE, stringsAsFactors=FALSE)

# Remove uninfected sample and the samples from the first two days of infection from the Eimeria data.
# The first two days have too few reads to be significant and would likely cause issues with the analysis
low_eimeria_samples <- annotation$File_name[match(c("0","1","2"),annotation$Timepoint_days)]
m = match(low_eimeria_samples,colnames(eimeria_reads))
eimeria_reads <- eimeria_reads[,-m]
eimeria_annotation <- annotation[-match(c("0","1","2"),annotation$Timepoint_days),]
rownames(eimeria_annotation) <- eimeria_annotation$File_name

chicken_annotation <- annotation
rownames(chicken_annotation) <- chicken_annotation$File_name

# Define the grouping for the chicken and Eimeria samples.  While there is an uninfected sample in the
# chicken, that group is defined as day 0 and so that won't be included in the groups.  Only timepoints
# will be used.
group_chicken <- factor(chicken_annotation[colnames(chicken_reads)[-1],"Timepoint_days"])
group_eimeria <- factor(eimeria_annotation[colnames(eimeria_reads)[-1],"Timepoint_days"])

# Run a PCA on the data to identify and remove outliers and other potential issues
pca_chicken <- prcomp(chicken_reads[,-1])
pca_eimeria <- prcomp(eimeria_reads[,-1])

ggbiplot(pca_chicken, labels = chicken_reads[,1]) + coord_cartesian(xlim = c(0, 150))
ggbiplot(pca_eimeria, labels = eimeria_reads[,1]) + coord_cartesian(xlim = c(0, 80))

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
design_chicken <- model.matrix(~group, data = chicken_dgelist_filt_norm$samples)
rownames(design_chicken) <- colnames(chicken_dgelist_filt_norm)

# Define the design of the eimeria analysis
design_eimeria <- model.matrix(~group, data = eimeria_dgelist_filt_norm$samples)
rownames(design_eimeria) <- colnames(eimeria_dgelist_filt_norm)

# Estimate the dispersion for each dataset and plot for the biological coefficient of variation (BCV)
#chicken_disp <- estimateDisp(chicken_dgelist_filt_norm, design_chicken, robust = TRUE)
#eimeria_disp <- estimateDisp(eimeria_dgelist_filt_norm, design_eimeria, robust = TRUE)

#plotBCV(chicken_disp)
#plotBCV(eimeria_disp)

# As there is no biological replication in this dataset, the BCV can't be estimated from it.  Therefore,
# it will be set at 0.2, since the manual indicates that this might be a reasonable value for this type
# of data
bcv <- 0.2

# Test for differential expression in chicken
chicken_1vs0 <- exactTest(chicken_dgelist_filt_norm, c(1,2), dispersion = bcv)
chicken_2vs0 <- exactTest(chicken_dgelist_filt_norm, c(1,3), dispersion = bcv)
chicken_3vs0 <- exactTest(chicken_dgelist_filt_norm, c(1,4), dispersion = bcv)
chicken_4vs0 <- exactTest(chicken_dgelist_filt_norm, c(1,5), dispersion = bcv)
chicken_10vs0 <- exactTest(chicken_dgelist_filt_norm, c(1,6), dispersion = bcv)

chicken_test_list <- list(chicken_1vs0,chicken_2vs0,chicken_3vs0,chicken_4vs0,chicken_10vs0)
topgenes_chicken_list <- lapply(chicken_test_list,
                                function(x) topTags(x, n = 18670))

# Test for differential expression in Eimeria
eimeria_4vs3 <- exactTest(eimeria_dgelist_filt_norm, c(1,2), dispersion = bcv)
eimeria_10vs3 <- exactTest(eimeria_dgelist_filt_norm, c(1,3), dispersion = bcv)

eimeria_test_list <- list(eimeria_4vs3,eimeria_10vs3)
topgenes_eimeria_list <- lapply(eimeria_test_list,
                                function(x) topTags(x, n = 1000))

# Volcano plots for all comparisons
timepoints <- c("1","2","3","4","10")
i = 1
while (i <= length(chicken_test_list)) {
  p <- EnhancedVolcano(chicken_test_list[[i]]$table,
                       lab = rownames(chicken_test_list[[i]]$table),
                       title = paste('Chicken in vivo, day', timepoints[i], 'vs day 0'),
                       x = 'logFC',
                       y = 'PValue')
  print(p)
  i = i + 1
}

i = 1
while (i <= length(eimeria_test_list)) {
  p <- EnhancedVolcano(eimeria_test_list[[i]]$table,
                       lab = rownames(eimeria_test_list[[i]]$table),
                       title = paste('Eimeria in vivo, days', timepoints[i+3], 'vs 3'),
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
go_chicken_list <- lapply(chicken_test_list, 
                          function(x) goana(x, geneid = entrez_geneIDs, FDR = 0.05, species = "Gg"))
kegg_chicken_list <- lapply(chicken_test_list, 
                            function(x) kegga(x, geneid = entrez_geneIDs, FDR = 0.05, species.KEGG = "gga"))
topgo_chicken_up_list <- lapply(go_chicken_list, 
                                function(x) topGO(x, ont="BP", sort="Up", n=30))
topgo_chicken_down_list <- lapply(go_chicken_list, 
                                  function(x) topGO(x, ont="BP", sort="Down", n=30))
topkegg_chicken_up_list <- lapply(kegg_chicken_list,
                                  function(x) topKEGG(x, sort="Up", n=30))
topkegg_chicken_down_list <- lapply(kegg_chicken_list,
                                    function(x) topKEGG(x, sort="Down", n=30))






# Get the CPM values for all genes in all samples
chicken_cpm <- cpm.DGEList(chicken_dgelist_filt_norm)
eimeria_cpm <- cpm.DGEList(eimeria_dgelist_filt_norm)

df1 <- data.frame(genes=colnames(eimeria_cpm),cpm_counts=log(eimeria_cpm["ETH_00034060",], base = 10))
df2 <- data.frame(genes=colnames(eimeria_cpm),cpm_counts=log(eimeria_cpm["ETH_00032220",], base = 10))
df3 <- data.frame(genes=colnames(eimeria_cpm),cpm_counts=log(eimeria_cpm["ETH_00014680",], base = 10))
df4 <- data.frame(genes=colnames(eimeria_cpm),cpm_counts=log(eimeria_cpm["ETH_00007730",], base = 10))

ggplot(data=df1, aes(x=genes, y=cpm_counts)) +
  geom_bar(stat="identity") +
  coord_cartesian(ylim = c(-2, 4)) +
  scale_x_discrete(limits=c("SH-2259-32_S4", "SH-2259-11_S5", "SH-2259-34_S6"))
ggplot(data=df2, aes(x=genes, y=cpm_counts)) +
  geom_bar(stat="identity") +
  coord_cartesian(ylim = c(-2, 4)) +
  scale_x_discrete(limits=c("SH-2259-32_S4", "SH-2259-11_S5", "SH-2259-34_S6"))
ggplot(data=df3, aes(x=genes, y=cpm_counts)) +
  geom_bar(stat="identity") +
  coord_cartesian(ylim = c(-2, 4)) +
  scale_x_discrete(limits=c("SH-2259-32_S4", "SH-2259-11_S5", "SH-2259-34_S6"))
ggplot(data=df4, aes(x=genes, y=cpm_counts)) +
  geom_bar(stat="identity") +
  coord_cartesian(ylim = c(-2, 4)) +
  scale_x_discrete(limits=c("SH-2259-32_S4", "SH-2259-11_S5", "SH-2259-34_S6"))











