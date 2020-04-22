# A pipeline for the differential expression analysis of a dual RNA-seq dataset of chickens and Eimeria tenella parsites

library(edgeR)
library(EnhancedVolcano)
library(org.Gg.eg.db)
library(GO.db)
library(ggplot2)
library(ggbiplot)
library(KEGGREST)
library(metacoder)
library(WGCNA)
library(gplots)
library(pals)
library(ClassDiscovery)
library(RColorBrewer)

# Specify the path to the read files and metadata and read them into R
path_to_chicken_reads <- "results/htseq/reverse/in_vitro_fusion/processed_reads/chicken_counts.csv"
path_to_eimeria_reads <- "results/htseq/reverse/in_vitro_fusion/processed_reads/eimeria_counts.csv"
path_to_annotation_file <- "results/htseq/reverse/in_vitro_fusion/processed_reads/metadata_table.csv"

chicken_reads <- read.delim(path_to_chicken_reads, sep = ",", check.names=FALSE, stringsAsFactors=FALSE)
eimeria_reads <- read.delim(path_to_eimeria_reads, sep = ",", check.names=FALSE, stringsAsFactors=FALSE)
annotation <- read.delim(path_to_annotation_file, sep = ",", check.names=FALSE, stringsAsFactors=FALSE)

# Remove uninfected samples from the Eimeria data and the sample with only Eimeria data from the chicken data
eimeria_annotation <- annotation[annotation$Infection_status == "Infected",]
eimeria_reads <- eimeria_reads[,c("entrez_gene_id", "locus_tag", eimeria_annotation$File_name)]
rownames(eimeria_annotation) <- eimeria_annotation$File_name

eimeria_sample <- annotation$File_name == "33_S29"
chicken_annotation <- annotation[!eimeria_sample,]
chicken_reads$`33_S29` <- NULL
rownames(chicken_annotation) <- chicken_annotation$File_name

# Define the grouping for the chicken and Eimeria samples.  For the chicken, each timepoint and infection status 
# is its own group while only time is a factor for the Eimeria samples
group_chicken <- factor(paste(substr(chicken_annotation[colnames(chicken_reads)[-c(1,2)],"Infection_status"],1,1),
                              chicken_annotation[colnames(chicken_reads)[-c(1,2)],"Timepoint_hours"],
                              sep="."))
group_eimeria <- factor(eimeria_annotation[colnames(eimeria_reads)[-c(1,2)],"Timepoint_hours"])

# Run a PCA on the data to identify, examine and deali with outliers and other potential issues
chicken_cpm_raw <- cpm(chicken_reads[,-c(1,2)])
eimeria_cpm_raw <- cpm(eimeria_reads[,-c(1,2)])

rownames(chicken_cpm_raw) <- chicken_reads[,2]
rownames(eimeria_cpm_raw) <- eimeria_reads[,1]

pca_chicken <- prcomp(chicken_cpm_raw)
pca_eimeria <- prcomp(eimeria_cpm_raw)

ggbiplot(pca_chicken, labels = chicken_reads[,2]) + coord_cartesian(xlim = c(0, 160)) + 
  labs(title = "PCA of normalized chicken read counts per gene") +
  theme(plot.title = element_text(hjust = 0.5))
ggbiplot(pca_eimeria, labels = eimeria_reads[,1]) + coord_cartesian(xlim = c(0, 60)) + 
  labs(title = "PCA of normalized E. tenella read counts per gene") +
  theme(plot.title = element_text(hjust = 0.5))

# Investigate eimeria outliers
df1 <- data.frame(groups=group_eimeria,log10_cpm_counts=log(eimeria_cpm_raw["25249494",], base = 10))
df1_mean <- aggregate(x = df1$log10_cpm_counts, by = list(df1$groups), mean)
colnames(df1_mean) <- colnames(df1)
df2 <- data.frame(groups=group_eimeria,log10_cpm_counts=log(eimeria_cpm_raw["25252472",], base = 10))
df2_mean <- aggregate(x = df2$log10_cpm_counts, by = list(df2$groups), mean)
colnames(df2_mean) <- colnames(df2)
df3 <- data.frame(groups=group_eimeria,log10_cpm_counts=log(eimeria_cpm_raw["25250679",], base = 10))
df3_mean <- aggregate(x = df3$log10_cpm_counts, by = list(df3$groups), mean)
colnames(df3_mean) <- colnames(df3)
df4 <- data.frame(groups=group_eimeria,log10_cpm_counts=log(eimeria_cpm_raw["25253446",], base = 10))
df4_mean <- aggregate(x = df4$log10_cpm_counts, by = list(df4$groups), mean)
colnames(df4_mean) <- colnames(df4)
df5 <- data.frame(groups=group_chicken,log10_cpm_counts=log(chicken_cpm_raw["112533599",], base = 10))
df5_mean <- aggregate(x = df5$log10_cpm_counts, by = list(df5$groups), mean)
colnames(df5_mean) <- colnames(df5)
df6 <- data.frame(groups=group_chicken,log10_cpm_counts=log(chicken_cpm_raw["770705",], base = 10))
df6_mean <- aggregate(x = df6$log10_cpm_counts, by = list(df6$groups), mean)
colnames(df6_mean) <- colnames(df6)

ggplot(data=df1_mean, aes(x=groups, y=log10_cpm_counts)) +
  geom_bar(stat="identity") +
  coord_cartesian(ylim = c(0, 6))
ggplot(data=df2_mean, aes(x=groups, y=log10_cpm_counts)) +
  geom_bar(stat="identity") +
  coord_cartesian(ylim = c(0, 6))
ggplot(data=df3_mean, aes(x=groups, y=log10_cpm_counts)) +
  geom_bar(stat="identity") +
  coord_cartesian(ylim = c(0, 6))
ggplot(data=df4_mean, aes(x=groups, y=log10_cpm_counts)) +
  geom_bar(stat="identity") +
  coord_cartesian(ylim = c(0, 6))
ggplot(data=df5_mean, aes(x=groups, y=log10_cpm_counts)) +
  geom_bar(stat="identity") +
  coord_cartesian(ylim = c(0, 6)) +
  scale_x_discrete(limits=c("U.0","U.2","I.2","U.4","I.4","U.12","I.12","U.24","I.24","U.48","I.48","U.72","I.72"))
ggplot(data=df6_mean, aes(x=groups, y=log10_cpm_counts)) +
  geom_bar(stat="identity") +
  coord_cartesian(ylim = c(0, 6)) +
  scale_x_discrete(limits=c("U.0","U.2","I.2","U.4","I.4","U.12","I.12","U.24","I.24","U.48","I.48","U.72","I.72"))

# Create the DGElist object containing the count data in a format that edgeR can work with
chicken_dgelist <- DGEList(counts = chicken_reads[,3:length(chicken_reads)], genes = chicken_reads[,c(1,2)], group = group_chicken)
eimeria_dgelist <- DGEList(counts = eimeria_reads[,3:length(eimeria_reads)], genes = eimeria_reads[,c(1,2)], group = group_eimeria)

rownames(chicken_dgelist$counts) <- rownames(chicken_dgelist$genes) <- chicken_dgelist$genes[,1]
rownames(eimeria_dgelist$counts) <- rownames(eimeria_dgelist$genes) <- eimeria_dgelist$genes[,2]

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

# Check MDS plots of logFC values for outliers or other potential issues
plotMDS(chicken_dgelist_filt_norm)
plotMDS(eimeria_dgelist_filt_norm)

# Check PCA plots of CPM values for each sample for outliers and other potential issues
chicken_cpm <- cpm.DGEList(chicken_dgelist_filt_norm)
eimeria_cpm <- cpm.DGEList(eimeria_dgelist_filt_norm)

pca_cpm_chicken <- prcomp(t(chicken_cpm))
pca_cpm_eimeria <- prcomp(t(eimeria_cpm[,-11]))

ggbiplot(pca_cpm_chicken, var.axes = FALSE, choices = 1:2, alpha = 1) +
  geom_text(aes(label=chicken_dgelist_filt_norm$samples$group),hjust=0.5, vjust=-0.5) +
  ggtitle("PCA of normalized chicken read counts per sample") +
  theme(plot.title = element_text(hjust = 0.5))
ggbiplot(pca_cpm_eimeria, var.axes = FALSE, choices = 1:2, alpha = 1) +
  geom_text(aes(label=eimeria_dgelist_filt_norm$samples$group[-11]),hjust=0.5, vjust=-0.5) +
  ggtitle("PCA of normalized E. tenella read counts per sample") +
  theme(plot.title = element_text(hjust = 0.5))

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
fit_chicken <- glmQLFit(chicken_disp, design_chicken, prior.count = 0.125)
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
                                function(x) topTags(x, n = dim(chicken_qlf_list[[1]]$table)[1]))

# Test for differential expression in Eimeria
fit_eimeria <- glmQLFit(eimeria_disp, design_eimeria)

qlf_eimeria_2 <- glmQLFTest(fit_eimeria, coef = 2)
qlf_eimeria_4 <- glmQLFTest(fit_eimeria, coef = 3)
qlf_eimeria_12 <- glmQLFTest(fit_eimeria, coef = 4)
qlf_eimeria_24 <- glmQLFTest(fit_eimeria, coef = 5)
qlf_eimeria_48 <- glmQLFTest(fit_eimeria, coef = 6)
qlf_eimeria_72 <- glmQLFTest(fit_eimeria, coef = 7)
qlf_eimeria_all <- glmQLFTest(fit_eimeria, coef = 2:7)

eimeria_qlf_list <- list(qlf_eimeria_2,qlf_eimeria_4,qlf_eimeria_12,qlf_eimeria_24,qlf_eimeria_48,
                      qlf_eimeria_72)
topgenes_eimeria_list <- lapply(eimeria_qlf_list,
                                function(x) topTags(x, n = dim(eimeria_qlf_list[[1]]$table)[1]))

# Define FDR threshold for significance
fdr_threshold <- 0.05
logfc_threshold <- 1

# Volcano plots for all comparisons
timepoints <- c("2h","4h","12h","24h","48h","72h")
i = 1
while (i <= length(topgenes_chicken_list)) {
  p <- EnhancedVolcano(topgenes_chicken_list[[i]]$table,
                  lab = rownames(topgenes_chicken_list[[i]]$table),
                  title = paste('Chicken in vitro', timepoints[i]),
                  x = 'logFC',
                  y = 'FDR',
                  xlim = c(-10,10),
                  ylim = c(0,10),
                  pCutoff = fdr_threshold,
                  FCcutoff = logfc_threshold)
  print(p)
  i = i + 1
}

i = 1
while (i <= length(topgenes_eimeria_list)) {
  p <- EnhancedVolcano(topgenes_eimeria_list[[i]]$table,
                       lab = rownames(topgenes_eimeria_list[[i]]$table),
                       title = paste('Eimeria in vitro', timepoints[i]),
                       x = 'logFC',
                       y = 'FDR',
                       xlim = c(-13,13),
                       ylim = c(0,40),
                       pCutoff = fdr_threshold,
                       FCcutoff = logfc_threshold)
  print(p)
  i = i + 1
}


# Compute the number of unique genes differentially expressed at each time point for both datasets, which genes they are and
# which genes are differentially expressed at all time points
get_unique_and_shared_de_genes <- function(de_gene_list, id_col_num = 1) {
  unique_de_gene_list <- de_gene_list
  shared_genes <- c()
  num_shared_genes <- 0
  
  i <- 0
  while (i < length(de_gene_list)) {
    i <- i + 1
    k <- 0
    while (k < dim(unique_de_gene_list[[i]])[1]) {
      j <- 0
      num_samples <- 0
      while (j < length(de_gene_list)) {
        j <- j + 1
        if (j == i) {
          next()
        } else {
          name_in_k <- unique_de_gene_list[[i]][k+1,id_col_num] %in% de_gene_list[[j]][,id_col_num]
        }
        if (name_in_k) {
          num_samples <- num_samples + 1
        }
      }
      if (num_samples == 5 && !(unique_de_gene_list[[i]][k+1,id_col_num] %in% shared_genes)) {
        num_shared_genes <- num_shared_genes + 1
        shared_genes[num_shared_genes] <- unique_de_gene_list[[i]][k+1,id_col_num]
      }
      if (num_samples > 0) {
        unique_de_gene_list[[i]] <- unique_de_gene_list[[i]][-(k+1),]
      } else {
        k <- k + 1
      }
    }
  }
  return(list(unique_de_gene_list, shared_genes))
}


de_genes_chicken <- lapply(topgenes_chicken_list, function(x) x$table[x$table$FDR < fdr_threshold,])
de_genes_chicken_pos <- lapply(de_genes_chicken, function(x) x[x$logFC > logfc_threshold,])
de_genes_chicken_neg <- lapply(de_genes_chicken, function(x) x[x$logFC < -logfc_threshold,])
de_genes_chicken_all <- lapply(de_genes_chicken, function(x) x[abs(x$logFC) > logfc_threshold,])

unique_chicken_de_genes_pos <- get_unique_and_shared_de_genes(de_genes_chicken_pos, 1)
unique_chicken_de_genes_neg <- get_unique_and_shared_de_genes(de_genes_chicken_neg, 1)
unique_chicken_de_genes_all <- get_unique_and_shared_de_genes(de_genes_chicken_all, 1)

de_genes_eimeria <- lapply(topgenes_eimeria_list, function(x) x$table[x$table$FDR < fdr_threshold,])
de_genes_eimeria_pos <- lapply(de_genes_eimeria, function(x) x[x$logFC > logfc_threshold,])
de_genes_eimeria_neg <- lapply(de_genes_eimeria, function(x) x[x$logFC < -logfc_threshold,])
de_genes_eimeria_all <- lapply(de_genes_eimeria, function(x) x[abs(x$logFC) > logfc_threshold,])

unique_eimeria_de_genes_pos <- get_unique_and_shared_de_genes(de_genes_eimeria_pos, 2)
unique_eimeria_de_genes_neg <- get_unique_and_shared_de_genes(de_genes_eimeria_neg, 2)
unique_eimeria_de_genes_all <- get_unique_and_shared_de_genes(de_genes_eimeria_all, 2)

# Get a list of all genes that are significantly DE at any timepoint and a list of those that are DE and have a logFC > 1
egGENENAME <- toTable(org.Gg.egGENENAME)
i <- 0
while (i < length(de_genes_chicken)) {
  i <- i + 1
  df_temp <- de_genes_chicken[[i]][,c("gene_name", "entrez_gene_id", "FDR")]
  if (i == 1) {
    all_de_genes_chicken <- df_temp
  } else {
    all_de_genes_chicken <- rbind(all_de_genes_chicken, df_temp)
  }
}
all_de_genes_chicken <- all_de_genes_chicken[order(all_de_genes_chicken$FDR),]
all_de_genes_chicken <- all_de_genes_chicken[!duplicated(all_de_genes_chicken[,c("gene_name","entrez_gene_id")]),]
all_de_genes_chicken$description <- egGENENAME[match(all_de_genes_chicken$entrez_gene_id, egGENENAME$gene_id),]$gene_name

write.table(all_de_genes_chicken, "results/de_analysis/all_de_genes_chicken.csv", sep = ",", row.names = FALSE)

i <- 0
while (i < length(de_genes_eimeria)) {
  i <- i + 1
  df_temp <- de_genes_eimeria[[i]][,c("locus_tag", "entrez_gene_id", "FDR")]
  if (i == 1) {
    all_de_genes_eimeria <- df_temp
  } else {
    all_de_genes_eimeria <- rbind(all_de_genes_eimeria, df_temp)
  }
}
all_de_genes_eimeria <- all_de_genes_eimeria[order(all_de_genes_eimeria$FDR),]
all_de_genes_eimeria <- all_de_genes_eimeria[!duplicated(all_de_genes_eimeria[,c("locus_tag","entrez_gene_id")]),]

write.table(all_de_genes_eimeria, "results/de_analysis/all_de_genes_eimeria.csv", sep = ",", row.names = FALSE)

i <- 0
while (i < length(de_genes_chicken_all)) {
  i <- i + 1
  df_temp <- de_genes_chicken_all[[i]][,c("gene_name", "entrez_gene_id", "FDR")]
  if (i == 1) {
    logfc_filt_de_genes_chicken <- df_temp
  } else {
    logfc_filt_de_genes_chicken <- rbind(logfc_filt_de_genes_chicken, df_temp)
  }
}
logfc_filt_de_genes_chicken <- logfc_filt_de_genes_chicken[order(logfc_filt_de_genes$FDR),]
logfc_filt_de_genes_chicken <- logfc_filt_de_genes_chicken[!duplicated(logfc_filt_de_genes_chicken[,c("gene_name","entrez_gene_id")]),]
logfc_filt_de_genes_chicken$description <- egGENENAME[match(logfc_filt_de_genes_chicken$entrez_gene_id, egGENENAME$gene_id),]$gene_name

write.table(logfc_filt_de_genes_chicken, "results/de_analysis/logfc_filt_de_genes_chicken.csv", sep = ",", row.names = FALSE)


i <- 0
while (i < length(de_genes_eimeria_all)) {
  i <- i + 1
  df_temp <- de_genes_eimeria_all[[i]][,c("locus_tag", "entrez_gene_id", "FDR")]
  if (i == 1) {
    logfc_filt_de_genes_eimeria <- df_temp
  } else {
    logfc_filt_de_genes_eimeria <- rbind(logfc_filt_de_genes_eimeria, df_temp)
  }
}
logfc_filt_de_genes_eimeria <- logfc_filt_de_genes_eimeria[order(logfc_filt_de_genes$FDR),]
logfc_filt_de_genes_eimeria <- logfc_filt_de_genes_eimeria[!duplicated(logfc_filt_de_genes_eimeria[,c("locus_tag","entrez_gene_id")]),]

write.table(logfc_filt_de_genes_eimeria, "results/de_analysis/logfc_filt_de_genes_eimeria.csv", sep = ",", row.names = FALSE)


# Get the top 30 most significantly DE genes from each timepoint and export to a table
top_de_genes_path <- "results/de_analysis/top_de_gene_tables"
top_de_gene_file_name <- "top_de_genes"
logfc_filtered_topgenes_chicken_list <- lapply(topgenes_chicken_list, function(x) x[abs(x$table$logFC) > logfc_threshold, ])
logfc_filtered_topgenes_eimeria_list <- lapply(topgenes_eimeria_list, function(x) x[abs(x$table$logFC) > logfc_threshold, ])

i <- 0
while (i < length(timepoints)) {
  i <- i + 1
  top_de_genes_file_path_chicken <- paste(top_de_genes_path, "/", top_de_gene_file_name, "_chicken_", timepoints[i], ".csv", sep = "")
  top_de_genes_file_path_eimeria <- paste(top_de_genes_path, "/", top_de_gene_file_name, "_eimeria_", timepoints[i], ".csv", sep = "")
  write.table(logfc_filtered_topgenes_chicken_list[[i]][1:30,], top_de_genes_file_path_chicken, sep = ",", row.names = FALSE)
  write.table(logfc_filtered_topgenes_eimeria_list[[i]][1:30,], top_de_genes_file_path_eimeria, sep = ",", row.names = FALSE)
}

# Produce global heatmaps of gene expression for both chicken and E. tenella using both p-values and log fold change
create_gene_heatmap_df <- function(de_gene_list, heatmap_value, gene_names, timepoints) {
  # de_gene_list is a list of DEGlist objects after performing a DE analysis
  # heatmap_value is a string with a value of either "FDR" or "logFC" depending on what the heatmap should show
  # gene_names is a vector of gene names in the analysis
  # timepoints is a vector of time points at which the samples were taken
  mapped_val_df <- as.data.frame(matrix(0,
                                        ncol = length(de_gene_list),
                                        nrow = dim(de_gene_list[[1]])[1]))
  
  colnames(mapped_val_df) <- timepoints
  
  rownames(mapped_val_df) <- gene_names
  if (heatmap_value == "FDR") {
    i <- 0
    while (i < length(de_gene_list)) {
      i <- i + 1
      m <- match(gene_names, rownames(de_gene_list[[i]]$table))
      mapped_val_df[,i] <- -log10(de_gene_list[[i]]$table[m,]$FDR)
    }
  } else if(heatmap_value == "logFC") {
    i <- 0
    while (i < length(de_gene_list)) {
      i <- i + 1
      m <- match(gene_names, rownames(de_gene_list[[i]]$table))
      mapped_val_df[,i] <- de_gene_list[[i]]$table[m,]$logFC
    }
  }
  return(mapped_val_df)
}

filter_gene_heatmap_df <- function(gene_heatmap_df, de_gene_list) {
  m <- match(de_gene_list, rownames(gene_heatmap_df))
  m <- m[complete.cases(m)]
  return(gene_heatmap_df[m,])
}

chicken_fdr_df <- create_gene_heatmap_df(topgenes_chicken_list, "FDR", rownames(topgenes_chicken_list[[1]]$table), timepoints)
chicken_logfc_df <- create_gene_heatmap_df(topgenes_chicken_list, "logFC", rownames(topgenes_chicken_list[[1]]$table), timepoints)
eimeria_fdr_df <- create_gene_heatmap_df(topgenes_eimeria_list, "FDR", rownames(topgenes_eimeria_list[[1]]$table), timepoints)
eimeria_logfc_df <- create_gene_heatmap_df(topgenes_eimeria_list, "logFC", rownames(topgenes_eimeria_list[[1]]$table), timepoints)

chicken_fdr_df_filt <- filter_gene_heatmap_df(chicken_fdr_df, rownames(logfc_filt_de_genes_chicken))
chicken_logfc_df_filt <- filter_gene_heatmap_df(chicken_logfc_df, rownames(logfc_filt_de_genes_chicken))
eimeria_fdr_df_filt <- filter_gene_heatmap_df(eimeria_fdr_df, rownames(logfc_filt_de_genes_eimeria))
eimeria_logfc_df_filt <- filter_gene_heatmap_df(eimeria_logfc_df, rownames(logfc_filt_de_genes_eimeria))

heatmap_path <- "figures/heatmaps/chicken_genes_pval.pdf"
pdf(heatmap_path, width = 30, height = 240, pointsize = 60)
aspectHeatmap(as.matrix(chicken_fdr_df_filt), Colv = NA,
              xlab = "Timepoints", main = "Chicken gene p-values", 
              col = colorRampPalette(brewer.pal(9, "YlOrRd"))(500),
              hExp = 20, wExp = 1)
dev.off()

heatmap_path <- "figures/heatmaps/chicken_genes_logfc.pdf"
pdf(heatmap_path, width = 140, height = 180, pointsize = 100)
aspectHeatmap(as.matrix(chicken_logfc_df_filt), Colv = NA,
              xlab = "Timepoints", main = "Chicken gene log fold change", 
              col = colorRampPalette(brewer.pal(9, "RdYlGn"))(500),
              hExp = 5, wExp = 0.8)
dev.off()

heatmap_path <- "figures/heatmaps/eimeria_genes_pval.pdf"
pdf(heatmap_path, width = 130, height = 120, pointsize = 120)
aspectHeatmap(as.matrix(eimeria_fdr_df_filt), Colv = NA,
              xlab = "Timepoints", main = "E. tenella gene p-values", 
              col = colorRampPalette(brewer.pal(9, "YlOrRd"))(500),
              hExp = 5, wExp = 0.8)
dev.off()

heatmap_path <- "figures/heatmaps/eimeria_genes_logfc.pdf"
pdf(heatmap_path, width = 100, height = 100, pointsize = 100)
aspectHeatmap(as.matrix(eimeria_logfc_df_filt), Colv = NA,
              xlab = "Timepoints", main = "E. tenella gene log fold change", 
              col = colorRampPalette(brewer.pal(9, "RdYlGn"))(500),
              hExp = 5, wExp = 0.8)
dev.off()



# GO category and KEGG pathway analysis of the data.  For each list, the data is in order of increasing time
go_chicken_list <- lapply(chicken_qlf_list, 
                          function(x) goana(x, geneid = chicken_dgelist_filt_norm$genes$entrez_gene_id, FDR = fdr_threshold, species = "Gg"))
kegg_chicken_list <- lapply(chicken_qlf_list, 
                            function(x) kegga(x, geneid = chicken_dgelist_filt_norm$genes$entrez_gene_id, FDR = fdr_threshold, species.KEGG = "gga"))
topgo_chicken_up_list <- lapply(go_chicken_list, 
                                function(x) topGO(x, ont="BP", sort="Up", n=250))
topgo_chicken_down_list <- lapply(go_chicken_list, 
                                  function(x) topGO(x, ont="BP", sort="Down", n=250))
topkegg_chicken_up_list <- lapply(kegg_chicken_list,
                                  function(x) topKEGG(x, sort="Up", n=250))
topkegg_chicken_down_list <- lapply(kegg_chicken_list,
                                    function(x) topKEGG(x, sort="Down", n=250))


# GO and KEGG category heatmaps
create_cat_pval_df <- function(cat_list, cat_name, timepoints) {
  # cat_name_col is a vector containing names of categories
  cat_pval_up_df <- as.data.frame(matrix(0,
                                         ncol = length(cat_list),
                                         nrow = dim(cat_list[[1]])[1]))
  cat_pval_down_df <- as.data.frame(matrix(0,
                                           ncol = length(cat_list),
                                           nrow = dim(cat_list[[1]])[1]))
  cat_pval_mix_df <- as.data.frame(matrix(0,
                                          ncol = length(cat_list),
                                          nrow = dim(cat_list[[1]])[1]))
  colnames(cat_pval_up_df) <- timepoints
  colnames(cat_pval_down_df) <- timepoints
  colnames(cat_pval_mix_df) <- timepoints
  
  rownames(cat_pval_up_df) <- cat_name
  rownames(cat_pval_down_df) <- cat_name
  rownames(cat_pval_mix_df) <- cat_name
  
  i <- 0
  while (i < length(cat_list)) {
    i <- i + 1
    cat_pval_up_df[,i] <- -log10(cat_list[[i]]$P.Up)
    cat_pval_down_df[,i] <- -log10(cat_list[[i]]$P.Down)
  }
  
  i <- 0
  while (i < dim(cat_pval_up_df)[1]) {
    i <- i + 1
    j <- 0
    while (j < dim(cat_pval_up_df)[2]) {
      j <- j + 1
      if (cat_pval_up_df[i,j] >= cat_pval_down_df[i,j]) {
        cat_pval_mix_df[i,j] <- cat_pval_up_df[i,j]
      } else {
        cat_pval_mix_df[i,j] <- -cat_pval_down_df[i,j]
      }
    }
  }
  return(list(cat_pval_up_df, cat_pval_down_df, cat_pval_mix_df))
}

remove_high_pval_rows <- function(pval_df, pval_thresh) {
  pval_thresh <- -log10(pval_thresh)
  i <- 1
  while (i <= dim(pval_df)[1]) {
    num_high_pval <- 0
    j <- 0
    while (j < dim(pval_df)[2]) {
      j <- j + 1
      if (pval_df[i,j] < pval_thresh) {
        num_high_pval <- num_high_pval + 1
      }
    }
    if (num_high_pval == j) {
      pval_df <- pval_df[-i,]
    } else {
      i <- i + 1
    }
  }
  return(pval_df)
}

go_chicken_pval_list <- create_cat_pval_df(go_chicken_list, go_chicken_list[[1]]$Term, timepoints)
kegg_chicken_pval_list <- create_cat_pval_df(kegg_chicken_list, kegg_chicken_list[[1]]$Pathway, timepoints)

go_chicken_pval_list_filtered <- lapply(go_chicken_pval_list, function(x) remove_high_pval_rows(x, fdr_threshold/10))
kegg_chicken_pval_list_filtered <- lapply(kegg_chicken_pval_list, function(x) remove_high_pval_rows(x, fdr_threshold))

heatmap_path <- "figures/heatmaps/GO_cat_upregulated.pdf"
pdf(heatmap_path, width = 60, height = 180, pointsize = 60)
aspectHeatmap(as.matrix(go_chicken_pval_list_filtered[[1]]), Colv = NA,
              xlab = "Timepoints", main = "Upregulated GO categories", 
              hExp = 5, wExp = 0.8, col = colorRampPalette(brewer.pal(9, "YlOrRd"))(50))
dev.off()
heatmap_path <- "figures/heatmaps/GO_cat_downregulated.pdf"
pdf(heatmap_path, width = 140, height = 180, pointsize = 100)
aspectHeatmap(as.matrix(go_chicken_pval_list_filtered[[2]]), Colv = NA,
        xlab = "Timepoints", main = "Downregulated GO categories", 
        hExp = 2, wExp = 0.7, col = colorRampPalette(brewer.pal(9, "YlOrRd"))(50))
dev.off()

col_breaks <- c(seq(-10,-3,length = 15), seq(-2.95, 2.95, length = 119), seq(3,10,length = 15))
col_palette <- colorRampPalette(c("red", "white", "blue"))(n = 148)

heatmap_path <- "figures/heatmaps/GO_cat_mixed.pdf"
pdf(heatmap_path, width = 110, height = 70, pointsize = 60)
heatmap.2(as.matrix(go_chicken_pval_list_filtered[[3]]), Colv = FALSE, dendrogram = "row", col = col_palette, breaks = col_breaks,
        xlab = "Timepoints", main = "GO categories",cexRow=0.7,cexCol=1,margins=c(5,12),trace="none",srtCol=45,key = FALSE)
dev.off()

heatmap_path <- "figures/heatmaps/KEGG_cat_upregulated.pdf"
pdf(heatmap_path, width = 130, height = 120, pointsize = 120)
aspectHeatmap(as.matrix(kegg_chicken_pval_list_filtered[[1]]), Colv = NA,
              xlab = "Timepoints", main = "Upregulated KEGG categories", 
              wExp = 0.7, col = colorRampPalette(brewer.pal(9, "YlOrRd"))(50))
dev.off()
heatmap_path <- "figures/heatmaps/KEGG_cat_downregulated.pdf"
pdf(heatmap_path, width = 100, height = 100, pointsize = 100)
aspectHeatmap(as.matrix(kegg_chicken_pval_list_filtered[[2]]), Colv = NA,
              xlab = "Timepoints", main = "Downregulated KEGG categories", 
              hExp = 2, col = colorRampPalette(brewer.pal(9, "YlOrRd"))(50))
dev.off()

col_breaks <- c(seq(-10,-5,length = 11), seq(-4.95, 4.95, length = 199), seq(5,10,length = 11))
col_palette <- colorRampPalette(c("red", "white", "blue"))(n = 220)

heatmap_path <- "figures/heatmaps/KEGG_cat_mixed.pdf"
pdf(heatmap_path, width = 110, height = 70, pointsize = 100)
heatmap.2(as.matrix(kegg_chicken_pval_list_filtered[[3]]), Colv = FALSE, dendrogram = "row", col = col_palette, breaks = col_breaks,
          xlab = "Timepoints", main = "KEGG categories",cexRow=0.7,cexCol=1,margins=c(5,12),trace="none",srtCol=45)
dev.off()

# Get a lsit of genes that belong to each GO category of interest
get_de_go_cat_genes <- function(de_gene_list, go_id, database_go, database_symbols) {
  cat_genes <- unique(database_go$gene_id[database_go$go_id == go_id])
  cat_gene_ids <- database_symbols[match(cat_genes, database_symbols$gene_id),]$symbol
  cat_de_genes <- lapply(de_gene_list, function(x) x$table[cat_gene_ids,])
  cat_de_genes <- lapply(cat_de_genes, function(x) x[complete.cases(x),])
  return(cat_de_genes)
}

egGO <- toTable(org.Gg.egGO2ALLEGS)
egSYMBOL <- toTable(org.Gg.egSYMBOL)
cell_cell_signaling_by_wnt_gene_list <- get_de_go_cat_genes(topgenes_chicken_list, "GO:0198738", egGO, egSYMBOL)

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


plot_de_cat_genes <- function(de_cat_genes, experimental_conditions, plot_title, fdr_thresh = 0.05, 
                              logfc_thresh = 1, num_samples_fdr_thresh = 1, num_samples_logfc_thresh = 1) {
  # Produces a line plot of log2 fold change across a set of samples for genes in a specific
  # GO or KEGG category
  i <- 1
  j <- 1
  de_cat_gene_df <- as.data.frame(matrix(0, 
                                            ncol = dim(de_cat_genes[[1]])[2], 
                                            nrow = dim(de_cat_genes[[1]])[1]*length(de_cat_genes)))
  colnames(de_cat_gene_df) <- colnames(de_cat_genes[[1]])
  num_conditions <- length(de_cat_genes)
  num_genes <- dim(de_cat_genes[[1]])[1]
  while (i <= num_genes) {
    k <- 1
    while (k <= num_conditions) {
      de_cat_gene_df[j,] <- de_cat_genes[[k]][i,]
      j <- j + 1
      k <- k + 1
    }
    i <- i + 1
  }
  
  below_fdr_threshold <- de_cat_gene_df$FDR < fdr_thresh
  de_cat_gene_df$below_fdr_threshold <- below_fdr_threshold
  curr_num_genes <- dim(de_cat_genes[[1]])[1]
  neg_logfc_thresh <- -logfc_thresh
  
  i <- 0
  while (i < curr_num_genes) {
    gene_start <- num_conditions*i + 1
    gene_end <- num_conditions*(i + 1)
    gene_below_thresh <- FALSE
    j <- num_conditions*i + 1
    num_samples_below_fdr_thresh = 0
    num_samples_above_logfc_thresh = 0
    while (j <= num_conditions*(i+1)) {
      if (de_cat_gene_df$below_fdr_threshold[j]) {
        num_samples_below_fdr_thresh <- num_samples_below_fdr_thresh + 1
      }
      if (de_cat_gene_df$logFC[j] >= logfc_thresh || de_cat_gene_df$logFC[j] <= neg_logfc_thresh) {
        num_samples_above_logfc_thresh <- num_samples_above_logfc_thresh + 1
      }
      j <- j + 1
    }
    if (num_samples_below_fdr_thresh >= num_samples_fdr_thresh && num_samples_above_logfc_thresh >= num_samples_logfc_thresh) {
      i <- i + 1
    } 
    else {
      gene_start <- num_conditions*i + 1
      gene_end <- num_conditions*(i + 1)
      de_cat_gene_df <- de_cat_gene_df[-c(gene_start:gene_end),]
    }
    curr_num_genes <- dim(de_cat_gene_df)[1]/num_conditions
  }
  
  if (curr_num_genes > 30) {
    # To increase readability, the plot is split into two if there are too many genes being plotted
    split_vec_1 <- rep(1, curr_num_genes/2)
    split_vec_2 <- rep(2, curr_num_genes/2)
    split_vec <- c(split_vec_1, split_vec_1, split_vec_1, split_vec_1, split_vec_1, split_vec_1, 
                   split_vec_2, split_vec_2, split_vec_2, split_vec_2, split_vec_2, split_vec_2)
    if (curr_num_genes*6 > length(split_vec)) {
      split_vec <- c(split_vec, rep(2, num_conditions))
    }
    de_cat_gene_df_list <- split(de_cat_gene_df, split_vec)
    de_cat_gene_df <- de_cat_gene_df_list[[1]]
    curr_num_genes <- dim(de_cat_gene_df)[1]/num_conditions
    
    de_cat_gene_df$gene_timepoints <- rep(experimental_conditions, curr_num_genes)
    p <- ggplot(de_cat_gene_df, aes(x=gene_timepoints, y=logFC, group=gene_name)) +
      geom_line(aes(color = gene_name)) +
      scale_color_manual(values = as.vector(polychrome(curr_num_genes))) +
      geom_point(aes(shape=below_fdr_threshold, color = gene_name), size = 3) +
      scale_shape_manual(values = c(1,17)) +
      scale_x_discrete(limits=experimental_conditions) +
      theme_bw() +
      ggtitle(plot_title) +
      theme(plot.title = element_text(hjust = 0.5, size = 18)) +
      coord_cartesian(ylim = c(-5, 7.5), xlim = c(0, 80)) +
      labs(shape = "Below FDR threshold", color = "Gene symbol") +
      xlab("Sample time points") + ylab("log2(Fold change)") +
      geom_hline(yintercept = 0)
    print(p)
    
    de_cat_gene_df <- de_cat_gene_df_list[[2]]
    curr_num_genes <- dim(de_cat_gene_df)[1]/num_conditions
  }
    
  de_cat_gene_df$gene_timepoints <- rep(experimental_conditions, curr_num_genes)
  p <- ggplot(de_cat_gene_df, aes(x=gene_timepoints, y=logFC, group=gene_name)) +
    geom_line(aes(color = gene_name)) +
    scale_color_manual(values = as.vector(polychrome(curr_num_genes))) +
    geom_point(aes(shape=below_fdr_threshold, color = gene_name), size = 3) +
    scale_shape_manual(values = c(1,17)) +
    scale_x_discrete(limits=experimental_conditions) +
    theme_bw() +
    ggtitle(plot_title) +
    theme(plot.title = element_text(hjust = 0.5, size = 18)) +
    coord_cartesian(ylim = c(-5, 7.5), xlim = c(0, 80)) +
    labs(shape = "Below FDR threshold", color = "Gene symbol") +
    xlab("Sample time points") + ylab("log2(Fold change)") +
    geom_hline(yintercept = 0)
  print(p)
}

# Plot a lineplot showing log fold change in each sample for the genes in the toll like KEGG pathway
plot_de_cat_genes(toll_like_gene_list, c(2,4,12,24,48,72), "Toll-like receptor signaling pathway gene expression", fdr_thresh = fdr_threshold, 
                  logfc_thresh = logfc_threshold, num_samples_fdr_thresh = 1, num_samples_logfc_thresh = 1)
tlr15_list <- lapply(topgenes_chicken_list, function(x) x$table["TLR15",])
tlr21_list <- lapply(topgenes_chicken_list, function(x) x$table["TLR21",])
plot_de_cat_genes(tlr15_list, c(2,4,12,24,48,72), "TLR15 expression", fdr_thresh = fdr_threshold, 
                  logfc_thresh = logfc_threshold, num_samples_fdr_thresh = 0, num_samples_logfc_thresh = 0)
plot_de_cat_genes(tlr21_list, c(2,4,12,24,48,72), "TLR21 expression", fdr_thresh = fdr_threshold, 
                  logfc_thresh = logfc_threshold, num_samples_fdr_thresh = 0, num_samples_logfc_thresh = 0)


path_to_immuno_gene_list <- "immune_classified_genes.csv"
immune_classified_genes <- read.delim(path_to_immuno_gene_list, check.names=FALSE, stringsAsFactors=FALSE)
immune_cats <- unique(immune_classified_genes$Category)

i <- 0
while (i < length(immune_cats)) {
  i <- i + 1
  cat_gene_ids <- immune_classified_genes[immune_cats[i] == immune_classified_genes$Category,]$entrez_gene_id
  cat_gene_list <- lapply(topgenes_chicken_list, function(x) x$table[match(cat_gene_ids, x$table$entrez_gene_id),])
  cat_title <- paste("Immune category", immune_cats[i], "gene expression")
  cat_plot_path <- paste("figures/expression_line_plots/immune_cat_", gsub(" ", "_", gsub("/", "_", immune_cats[i])), ".pdf", sep = "")
  pdf(cat_plot_path, width = 19.2, height = 10.8)
  plot_de_cat_genes(cat_gene_list, c(2,4,12,24,48,72), cat_title, fdr_thresh = fdr_threshold, 
                    logfc_thresh = logfc_threshold, num_samples_fdr_thresh = 1, num_samples_logfc_thresh = 1)
  dev.off()
}


# Eimeria GO and KEGG analysis

# Import the file containing GO annotations downloaded from ToxoDB and KEGG annotations found through the KEGG KAAS service
path_to_eimeria_go_file <- "data/eimeria_go_kegg_annotation/GenesByGoTerm_Summary.csv"
path_to_eimeria_kegg_file <- "data/eimeria_go_kegg_annotation/eimeria_kegg_annotation.tsv"
path_to_transcript_to_gene_file <- "data/reference_genomes/eimeria_transcript_to_gene.tsv"
eimeria_go_annotation <- read.delim(path_to_eimeria_go_file, sep = ",", check.names=FALSE, stringsAsFactors=FALSE)
eimeria_kegg_best_cat_list <- read.delim(path_to_eimeria_kegg_file, check.names=FALSE, stringsAsFactors=FALSE, header = FALSE)
eimeria_transcript_to_gene <- read.delim(path_to_transcript_to_gene_file, check.names=FALSE, stringsAsFactors=FALSE, header = FALSE)

# Create a new dataframe with seperate entries for each gene for each GO category it belongs to
eimeria_full_go <- as.data.frame(matrix(0, 
                                        ncol = ncol(eimeria_go_annotation), 
                                        nrow = nrow(eimeria_go_annotation)))
colnames(eimeria_full_go) <- colnames(eimeria_go_annotation)

i <- 0
j <- 1
len_eimeria_go <- dim(eimeria_go_annotation)[1]
while (i < len_eimeria_go) {
  go_cats <- strsplit(eimeria_go_annotation$`Computed GO Process IDs`[i+1], ";")[[1]]
  go_cats_description <- strsplit(eimeria_go_annotation$`Computed GO Processes`[i+1], ";")[[1]]
  k <- 1
  len_go_cats <- length(go_cats)
  while (k <= len_go_cats) {
    eimeria_full_go[j,] <- eimeria_go_annotation[i+1,]
    eimeria_full_go[j,"Computed GO Process IDs"] <- go_cats[k]
    eimeria_full_go[j,"Computed GO Processes"] <- go_cats_description[k]
    k <- k + 1
    j <- j + 1
  }
  i <- i + 1
}

eimeria_go_analysis <- function(topgenes_eimeria, eimeria_full_go, eimeria_go_cats, fdr_thresh, logfc_thresh, ont) {
  # Function for computing gene membership in GO categories
  FDR_significant_genes <- topgenes_eimeria$table[topgenes_eimeria$table$FDR < fdr_thresh,]
  num_significant_genes <- c(dim(FDR_significant_genes[FDR_significant_genes$logFC > logfc_thresh,])[1],
                             dim(FDR_significant_genes[FDR_significant_genes$logFC < -logfc_thresh,])[1])
  eimeria_go_results <- as.data.frame(matrix(0, ncol = 8, nrow = length(eimeria_go_cats)))
  colnames(eimeria_go_results) <- c("GO_ID", "Term", "Ont", "N", "Up", "Down", "P.Up", "P.Down")
  i <- 1
  neg_logfc_thresh <- -logfc_thresh
  for (eimeria_go_cat in eimeria_go_cats) {
    genes_in_go <- eimeria_full_go[eimeria_full_go$`Computed GO Process IDs` == eimeria_go_cat,]
    eimeria_go_cat_term <- genes_in_go[1,"Computed GO Processes"]
    N <- dim(genes_in_go)[1]
    m <- match(genes_in_go$`Gene ID`, topgenes_eimeria$table$locus_tag)
    if (!is.na(m[1])) {
      m <- m[complete.cases(m)]
      topgenes_eimeria_go_cat <- topgenes_eimeria$table[m,]
      topgenes_eimeria_go_cat <- topgenes_eimeria_go_cat[topgenes_eimeria_go_cat$FDR < fdr_thresh,]
      if (dim(topgenes_eimeria_go_cat)[1] != 0) {
        num_up_genes <- nrow(topgenes_eimeria_go_cat[topgenes_eimeria_go_cat$logFC >= logfc_thresh,])
        num_down_genes <- nrow(topgenes_eimeria_go_cat[topgenes_eimeria_go_cat$logFC <= neg_logfc_thresh,])
      } else {
        num_up_genes <- 0
        num_down_genes <- 0
      }
    } else {
      num_up_genes <- 0
      num_down_genes <- 0
    }
    x <- c(num_up_genes, num_down_genes)
    n <- dim(topgenes_eimeria)[1] - N
    k <- num_significant_genes
    fisher_up <- fisher.test(matrix(c(x[1],k[1]-x[1],N-x[1],n-(k[1]-x[1])),nrow=2,ncol=2),alternative="greater")
    fisher_down <- fisher.test(matrix(c(x[2],k[2]-x[2],N-x[2],n-(k[2]-x[2])),nrow=2,ncol=2),alternative="greater")
    go_cat_results <- data.frame(eimeria_go_cat, eimeria_go_cat_term, ont, N, num_up_genes, num_down_genes, 
                                 fisher_up$p.value, fisher_down$p.value, stringsAsFactors = FALSE)
    eimeria_go_results[i,] <- go_cat_results[1,]
    i <- i + 1
  }
  rownames(eimeria_go_results) <- eimeria_go_results$GO_ID
  eimeria_go_results <- eimeria_go_results[,-1]
  return(eimeria_go_results)
}

# Get the number of genes that are up or down regulated in each GO category
eimeria_go_cats <- unique(eimeria_full_go$`Computed GO Process IDs`)
allgo_eimeria_list <- lapply(topgenes_eimeria_list, 
                             function(x) eimeria_go_analysis(x, eimeria_full_go, eimeria_go_cats, fdr_threshold, fdr_threshold, "BP"))




# Preprocessing of KEGG data for the KEGG enrichment analaysis
genes_in_analysis <- topgenes_eimeria_list[[1]]$table$entrez_gene_id
m <- match(genes_in_analysis, eimeria_transcript_to_gene[,1])
transcripts_in_analysis <- eimeria_transcript_to_gene[m,]
m <- match(transcripts_in_analysis[,2], eimeria_kegg_best_cat_list[,1])
kegg_ids_in_analysis <- eimeria_kegg_best_cat_list[m,2]
eimeria_kegg_genes <- data.frame(transcripts_in_analysis, kegg_ids_in_analysis)
eimeria_kegg_genes <- eimeria_kegg_genes[eimeria_kegg_genes[,3] != "",]
eimeria_kegg_id_universe <- eimeria_kegg_best_cat_list[eimeria_kegg_best_cat_list[,2] != "",2]

gene_to_pathway <- keggLink("pathway", "ko")
pathway_to_term <- keggList("pathway")

m <- match(eimeria_kegg_genes[,3], substr(names(gene_to_pathway), 4, nchar(gene_to_pathway[1])))
eimeria_kegg_genes <- eimeria_kegg_genes[!is.na(m),]
gene_pathways <- gene_to_pathway[m[complete.cases(m)]]
m <- match(gene_pathways, names(pathway_to_term))
pathway_terms <- pathway_to_term[m]
eimeria_kegg_annotation <- data.frame(genes = eimeria_kegg_genes[,1], kegg_pathway = gene_pathways, 
                                      kegg_term = pathway_terms, stringsAsFactors = FALSE)

m <- match(eimeria_kegg_id_universe, substr(names(gene_to_pathway), 4, nchar(gene_to_pathway[1])))
gene_universe_pathways <- gene_to_pathway[m[complete.cases(m)]]
eimeria_kegg_cats <- unique(gene_universe_pathways)

eimeria_kegg_analysis <- function(topgenes_eimeria, eimeria_kegg_annotation, eimeria_kegg_cats, 
                                  eimeria_universe, fdr_thresh, logfc_thresh) {
  # Function for computing gene membership in KEGG categories
  # Uses KEGGREST to access data from the KEGG database
  FDR_significant_genes <- topgenes_eimeria$table[topgenes_eimeria$table$FDR < fdr_thresh,]
  num_significant_genes <- c(dim(FDR_significant_genes[FDR_significant_genes$logFC > logfc_thresh,])[1],
                             dim(FDR_significant_genes[FDR_significant_genes$logFC < -logfc_thresh,])[1])
  eimeria_kegg_results <- as.data.frame(matrix(0, ncol = 7, nrow = length(eimeria_kegg_cats)))
  colnames(eimeria_kegg_results) <- c("KEGG_ID", "Term", "N", "Up", "Down", "P.Up", "P.Down")
  i <- 1
  neg_logfc_thresh <- -logfc_thresh
  for (eimeria_kegg_cat in eimeria_kegg_cats) {
    genes_in_kegg <- eimeria_kegg_annotation[eimeria_kegg_annotation$kegg_pathway == eimeria_kegg_cat,]
    eimeria_kegg_cat_term <- genes_in_kegg[1,"kegg_term"]
    N <- length(eimeria_universe[eimeria_universe == eimeria_kegg_cat])
    m <- match(genes_in_kegg$genes, topgenes_eimeria$table$entrez_gene_id)
    if (!is.na(m[1])) {
      m <- m[complete.cases(m)]
      topgenes_eimeria_kegg_cat <- topgenes_eimeria$table[m,]
      topgenes_eimeria_kegg_cat <- topgenes_eimeria_kegg_cat[topgenes_eimeria_kegg_cat$FDR < fdr_thresh,]
      if (dim(topgenes_eimeria_kegg_cat)[1] != 0) {
        num_up_genes <- nrow(topgenes_eimeria_kegg_cat[topgenes_eimeria_kegg_cat$logFC >= logfc_thresh,])
        num_down_genes <- nrow(topgenes_eimeria_kegg_cat[topgenes_eimeria_kegg_cat$logFC <= neg_logfc_thresh,])
      } else {
        num_up_genes <- 0
        num_down_genes <- 0
      }
    } else {
      num_up_genes <- 0
      num_down_genes <- 0
    }
    x <- c(num_up_genes, num_down_genes)
    n <- dim(topgenes_eimeria)[1] - N
    k <- num_significant_genes
    fisher_up <- fisher.test(matrix(c(x[1],k[1]-x[1],N-x[1],n-(k[1]-x[1])),nrow=2,ncol=2),alternative="greater")
    fisher_down <- fisher.test(matrix(c(x[2],k[2]-x[2],N-x[2],n-(k[2]-x[2])),nrow=2,ncol=2),alternative="greater")
    kegg_cat_results <- data.frame(eimeria_kegg_cat, eimeria_kegg_cat_term, N, num_up_genes, num_down_genes,
                                   fisher_up$p.value, fisher_down$p.value, stringsAsFactors = FALSE)
    eimeria_kegg_results[i,] <- kegg_cat_results[1,]
    i <- i + 1
  }
  rownames(eimeria_kegg_results) <- eimeria_kegg_results$KEGG_ID
  eimeria_kegg_results <- eimeria_kegg_results[,-1]
  return(eimeria_kegg_results)
}


allkegg_eimeria_list <- lapply(topgenes_eimeria_list, 
                             function(x) eimeria_kegg_analysis(x, eimeria_kegg_annotation, eimeria_kegg_cats,
                                                               gene_universe_pathways, fdr_threshold, 1))






# WGCNA analysis of the results
# Analysis adapted from the WGCNA tutorial

options(stringsAsFactors = FALSE)

chicken_expr <- t(cpm(chicken_dgelist_filt_norm))

gsg = goodSamplesGenes(chicken_expr, verbose = 3)
gsg$allOK

if (!gsg$allOK) {
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = hclust(dist(chicken_expr), method = "average")
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

# Plot a line to show the cut
abline(h = 50000, col = "red")
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 50000, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
chicken_expr_cut = chicken_expr[keepSamples, ]
nGenes = ncol(chicken_expr_cut)
nSamples = nrow(chicken_expr_cut)

chicken_expr_trait <- chicken_annotation[-38,2:3]
i <- 0
while (i < length(chicken_expr_trait$Infection_status)) {
  i <- i + 1
  if (chicken_expr_trait$Infection_status[i] == "Infected") {
    chicken_expr_trait$Infection_status[i] = 1
  } else {
    chicken_expr_trait$Infection_status[i] = 0
  }
}

# Choose a set of soft-thresholding powers
powers = c(c(1:9), seq(from = 10, to=100, by=10))
# Call the network topology analysis function
sft = pickSoftThreshold(chicken_expr_cut, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.70,col="red")
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


#net6 = blockwiseModules(chicken_expr_cut, power = 6,
#                       TOMType = "unsigned", minModuleSize = 30,
#                       reassignThreshold = 0, mergeCutHeight = 0.25,
#                       numericLabels = TRUE, pamRespectsDendro = FALSE,
#                       saveTOMs = TRUE,
#                       saveTOMFileBase = "chickenTOM6", 
#                       verbose = 3,
#                       maxBlockSize = 14000)

#net30 = blockwiseModules(chicken_expr_cut, power = 30,
#                        TOMType = "unsigned", minModuleSize = 30,
#                        reassignThreshold = 0, mergeCutHeight = 0.25,
#                        numericLabels = TRUE, pamRespectsDendro = FALSE,
#                        saveTOMs = TRUE,
#                        saveTOMFileBase = "chickenTOM30", 
#                        verbose = 3,
#                        maxBlockSize = 14000)

table(net6$colors)
table(net30$colors)

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net6$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net6$dendrograms[[1]], mergedColors[net6$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net6$colors
moduleColors = labels2colors(net6$colors)
MEs = net6$MEs
geneTree = net6$dendrograms[[1]]
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "chickenR6-networkConstruction-auto.RData")

# Define numbers of genes and samples
nGenes = ncol(chicken_expr_cut)
nSamples = nrow(chicken_expr_cut)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(chicken_expr_cut, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, chicken_expr_trait, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(chicken_expr_trait),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


# Define variable infection_status containing the Infection_status column of chicken_expr_trait
infection_status = as.data.frame(chicken_expr_trait$Infection_status);
names(infection_status) = "infection_status"

# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(chicken_expr_cut, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(chicken_expr_cut, infection_status, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(infection_status), sep="");
names(GSPvalue) = paste("p.GS.", names(infection_status), sep="");

module = "midnightblue"
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for infection status",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


colnames(chicken_expr_cut)[moduleColors=="brown"]

symbol2entrez = match(chicken_dgelist_filt_norm$genes$gene_name, colnames(chicken_expr_cut))

# Create the starting data frame
geneInfo0 <- data.frame(geneSymbol = colnames(chicken_expr_cut),
                       EntrezID = chicken_dgelist_filt_norm$genes$entrez_gene_id[],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for infection_status
modOrder <- order(-abs(cor(MEs, infection_status, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.infection_status));
geneInfo = geneInfo0[geneOrder, ]

write.csv(geneInfo, file = "geneInfo.csv")



# Get the corresponding Locuis Link IDs
allLLIDs <- geneInfo$EntrezID

GOenr = GOenrichmentAnalysis(moduleColors, allLLIDs, organism = "chicken", nBestP = 10)

tab = GOenr$bestPTerms[[4]]$enrichment


  
  
  
  
# Analysis of a fusion of Eimeria and chicken data

options(stringsAsFactors = FALSE)

# Load the chicken and E. tenella reads, concatenate them and remove uninfected samples
chicken_reads <- read.delim(path_to_chicken_reads, sep = ",", check.names=FALSE, stringsAsFactors=FALSE)
eimeria_reads <- read.delim(path_to_eimeria_reads, sep = ",", check.names=FALSE, stringsAsFactors=FALSE)
annotation <- read.delim(path_to_annotation_file, sep = ",", check.names=FALSE, stringsAsFactors=FALSE)
colnames(eimeria_reads) <- c("gene_name", colnames(eimeria_reads)[-1])

fusion_reads <- rbind(chicken_reads, eimeria_reads)

fusion_annotation <- annotation[annotation$Infection_status == "Infected",]
fusion_reads <- fusion_reads[,c("gene_name", "entrez_gene_id", fusion_annotation$File_name)]
rownames(fusion_annotation) <- fusion_annotation$File_name

group_fusion <- factor(fusion_annotation[colnames(fusion_reads)[-c(1,2)],"Timepoint_hours"])

# Create the DGElist object containing the count data in a format that edgeR can work with
fusion_dgelist <- DGEList(counts = fusion_reads[,3:length(fusion_reads)], genes = fusion_reads[,c(1,2)], group = group_fusion)

rownames(fusion_dgelist$counts) <- rownames(fusion_dgelist$genes) <- fusion_dgelist$genes[,1]

# Filter out genes that are lowly expressed in all samples
keep_fusion <- filterByExpr(fusion_dgelist)
fusion_dgelist_filt <- fusion_dgelist[keep_fusion, ,] 

# Setting the library size to the right value after filtering
fusion_dgelist_filt$samples$lib.size <- colSums(fusion_dgelist_filt$counts)

# Normalize gene counts to account for a few genes dominating the counts of the samples
fusion_dgelist_filt_norm <- calcNormFactors(fusion_dgelist_filt)

fusion_expr <- t(cpm(fusion_dgelist_filt_norm))

gsg = goodSamplesGenes(fusion_expr, verbose = 3)
gsg$allOK

if (!gsg$allOK) {
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = hclust(dist(fusion_expr), method = "average")
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

# Plot a line to show the cut
abline(h = 50000, col = "red")
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 50000, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
fusion_expr_cut = fusion_expr[keepSamples, ]
nGenes = ncol(fusion_expr_cut)
nSamples = nrow(fusion_expr_cut)

fusion_expr_trait <- fusion_annotation[-11,3]

# Choose a set of soft-thresholding powers
powers = c(c(1:9), seq(from = 10, to=100, by=5))
# Call the network topology analysis function
sft = pickSoftThreshold(fusion_expr_cut, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.70,col="red")
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


#net9 = blockwiseModules(fusion_expr_cut, power = 9,
#                        TOMType = "unsigned", minModuleSize = 30,
#                        reassignThreshold = 0, mergeCutHeight = 0.25,
#                        numericLabels = TRUE, pamRespectsDendro = FALSE,
#                        saveTOMs = TRUE,
#                        saveTOMFileBase = "fusionTOM9", 
#                        verbose = 3,
#                        maxBlockSize = 21000)

net18 = blockwiseModules(fusion_expr_cut, power = 18,
                        TOMType = "unsigned", minModuleSize = 30,
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs = TRUE,
                        saveTOMFileBase = "fusionTOM18", 
                        verbose = 3,
                        maxBlockSize = 21000)

table(net9$colors)
table(net18$colors)

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net9$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net9$dendrograms[[1]], mergedColors[net9$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net9$colors
moduleColors = labels2colors(net9$colors)
MEs = net9$MEs
geneTree = net9$dendrograms[[1]]
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "fusionR9-networkConstruction-auto.RData")

# Define numbers of genes and samples
nGenes = ncol(fusion_expr_cut)
nSamples = nrow(fusion_expr_cut)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(fusion_expr_cut, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, fusion_expr_trait, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = "Time point",
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


# Define variable infection_status containing the Infection_status column of fusion_expr_trait
time_point = as.data.frame(fusion_expr_trait);
names(time_point) = "Time_points"

# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(fusion_expr_cut, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(fusion_expr_cut, time_point, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(time_point), sep="");
names(GSPvalue) = paste("p.GS.", names(time_point), sep="");

module = "black"
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for infection status",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


colnames(fusion_expr_cut)[moduleColors=="brown"]

symbol2entrez = match(fusion_dgelist_filt_norm$genes$gene_name, colnames(fusion_expr_cut))

# Create the starting data frame
geneInfo0 <- data.frame(geneSymbol = colnames(fusion_expr_cut),
                        EntrezID = fusion_dgelist_filt_norm$genes$entrez_gene_id[],
                        moduleColor = moduleColors,
                        geneTraitSignificance,
                        GSPvalue)
# Order modules by their significance for time_point
modOrder <- order(-abs(cor(MEs, time_point, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership)) {
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Time_points));
geneInfo = geneInfo0[geneOrder, ]

write.csv(geneInfo, file = "geneInfo_fusion.csv", row.names = FALSE)

# GO and KEGG analyses of the modules

# Chicken gene to GO
# GO annotation file from the Gene Ontology Consortium
path_to_chicken_GO_annotation_file <- "data/chicken_go_annotation/goa_chicken.gaf"
chicken_go_annotation <- read.delim(path_to_chicken_GO_annotation_file, stringsAsFactors=FALSE, skip = 31, header = FALSE,
                                    col.names = c("DB", "DB_ID", "Symbol", "Qualifier", "GO_ID", "DB:Reference",
                                                  "Evidence_code", "With_or_from", "Aspect", "DB_Object_Name",
                                                  "DB_Object_Synonyms", "DB_Object_type", "Taxon", "Date", "Assigned_by",
                                                  "Annotation_extension", "Gene_product_form_ID"))
chicken_go <- chicken_go_annotation[,c(3,5)]
go_terms <- toTable(GOTERM)
m <- match(chicken_go$GO_ID, go_terms$go_id)
chicken_go$Term <- go_terms[m,"Term"]

# Eimeria gene to GO from ToxoDB
eimeria_go <- eimeria_full_go[,c(1,4,5)]
colnames(eimeria_go) <- colnames(chicken_go)
go_annotation <- rbind(chicken_go, eimeria_go)

module_cat_enrichment <- function(geneInfo, annotation, cat_id, num_genes) {
  # Function for finding the GO enrichment of a module from a WGCNA analysis of chicken and E. tenella data
  # geneInfo should contain gene symbols and module membership
  # annotation should contain gene symbols and  either GO category or KEGG pathway membershp for both chicken and E. tenella genes
  # cat_id is a string that tells what type of categories are being analysed
  # num_genes is the total number of genes in the analysis
  categories <- unique(annotation[,2])
  cat_analysis_results <- as.data.frame(matrix(0, ncol = 5, nrow = length(categories)))
  colnames(cat_analysis_results) <- c(paste(cat_id, "_", geneInfo[1,2], sep = ""), "Term", "N", "num_in_cat", "P_val")
  i <- 0
  while (i < length(categories)) {
    i <- i + 1
    genes_in_cat <- annotation[annotation[,2] == categories[i],]
    N <- dim(genes_in_cat)[1]
    m <- match(genes_in_cat[,1], geneInfo[,1])
    x <- length(m[complete.cases(m)])
    n <- num_genes - N
    k <- dim(geneInfo)[1]
    fisher_result <- fisher.test(matrix(c(x, k-x, N-x, n-(k-x)),nrow=2,ncol=2),alternative="greater")
    cat_results <- data.frame(categories[i], genes_in_cat[1,3], N, x, fisher_result$p.value, stringsAsFactors = FALSE)
    cat_analysis_results[i,] <- cat_results[1,]
  }
  print("Hello, world")
  return(cat_analysis_results)
}

geneInfo_mod <- list()
i <- 0
while (i < length(modNames)) {
  i <- i + 1
  geneInfo_mod[[i]] <- geneInfo[geneInfo$moduleColor == modNames[i],c(1,3)]
}

module_go_cats <- lapply(geneInfo_mod, function(x) module_cat_enrichment(x, go_annotation, "GO_ID", dim(geneInfo)[1]))

i <- 0
while (i < length(module_go_cats)) {
  i <- i + 1
  go_file_name <- paste("results/de_analysis/wgcna_module_go_kegg/", modNames[i], "_module_go_cats.csv", sep = "")
  write.csv(module_go_cats[[i]], file = go_file_name, row.names = FALSE)
}

# Chicken gene to KEGG
chicken_pathway_genes <- keggLink("pathway", "gga")
pathways <- keggList("pathway")

m <- match(substr(names(chicken_pathway_genes), 5, nchar(names(chicken_pathway_genes))), chicken_dgelist_filt_norm$genes$entrez_gene_id)
chicken_pathway_genes <- gsub("gga", "map", chicken_pathway_genes)

chicken_kegg_annotation <- data.frame(gene_symbol = chicken_dgelist$genes$gene_name[m],
                                      kegg_pathway = chicken_pathway_genes,
                                      kegg_term = pathways[match(substr(chicken_pathway_genes, 9, nchar(chicken_pathway_genes)), substr(names(pathways), 9, nchar(names(pathways))))])

chicken_kegg_annotation <- chicken_kegg_annotation[complete.cases(chicken_kegg_annotation),]

# Eimeria gene to KEGG from KAAS
eimeria_kegg_annotation

m <- match(eimeria_kegg_annotation$genes, eimeria_dgelist_filt_norm$genes$entrez_gene_id)

eimeria_kegg_annotation_symbol <- data.frame(gene_symbol = eimeria_dgelist$genes$locus_tag[m],
                                             kegg_pathway = eimeria_kegg_annotation$kegg_pathway,
                                             kegg_term = eimeria_kegg_annotation$kegg_term)

eimeria_kegg_annotation_symbol <- eimeria_kegg_annotation_symbol[complete.cases(eimeria_kegg_annotation_symbol),]

kegg_annotation <- rbind(chicken_kegg_annotation, eimeria_kegg_annotation_symbol)

module_kegg_cats <- lapply(geneInfo_mod, function(x) module_cat_enrichment(x, kegg_annotation, "KEGG_ID", dim(geneInfo)[1]))

i <- 0
while (i < length(module_go_cats)) {
  i <- i + 1
  go_file_name <- paste("results/de_analysis/wgcna_module_go_kegg/", modNames[i], "_module_kegg_cats.csv", sep = "")
  write.csv(module_go_cats[[i]], file = go_file_name, row.names = FALSE)
}






























