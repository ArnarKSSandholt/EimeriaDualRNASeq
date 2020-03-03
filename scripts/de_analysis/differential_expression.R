# A pipeline for the differential expression analysis of a dual RNA-seq dataset of chickens and Eimeria tenella parsites

library(edgeR)
library(EnhancedVolcano)
library(org.Gg.eg.db)

# Specify the files with the count data, the path to them and the labels
chicken_files <- c("1_S17_chicken_table.csv","2_S18_chicken_table.csv","3_S19_chicken_table.csv","4_S20_chicken_table.csv",
                   "5_S21_chicken_table.csv","6_S22_chicken_table.csv","7_S23_chicken_table.csv","8_S24_chicken_table.csv",
                   "9_S25_chicken_table.csv","10_S26_chicken_table.csv","1_S1_chicken_table.csv","2_S2_chicken_table.csv",
                   "3_S3_chicken_table.csv","4_S4_chicken_table.csv","5_S5_chicken_table.csv","6_S6_chicken_table.csv",
                   "7_S7_chicken_table.csv","8_S8_chicken_table.csv","9_S9_chicken_table.csv","11_S10_chicken_table.csv",
                   "12_S11_chicken_table.csv","13_S12_chicken_table.csv","14_S13_chicken_table.csv","15_S14_chicken_table.csv",
                   "18_S15_chicken_table.csv","19_S16_chicken_table.csv","20_S17_chicken_table.csv","21_S18_chicken_table.csv",
                   "23_S19_chicken_table.csv","24_S20_chicken_table.csv","25_S21_chicken_table.csv","26_S22_chicken_table.csv",
                   "27_S23_chicken_table.csv","28_S24_chicken_table.csv","29_S25_chicken_table.csv","30_S26_chicken_table.csv",
                   "31_S27_chicken_table.csv","32_S28_chicken_table.csv")
eimeria_files <- c("33_S29_eimeria_table.csv","2_S18_eimeria_table.csv","4_S20_eimeria_table.csv","6_S22_eimeria_table.csv",
                   "8_S24_eimeria_table.csv","10_S26_eimeria_table.csv","2_S2_eimeria_table.csv","4_S4_eimeria_table.csv",
                   "6_S6_eimeria_table.csv","8_S8_eimeria_table.csv","12_S11_eimeria_table.csv","14_S13_eimeria_table.csv",
                   "18_S15_eimeria_table.csv","20_S17_eimeria_table.csv","23_S19_eimeria_table.csv","24_S20_eimeria_table.csv",
                   "27_S23_eimeria_table.csv","28_S24_eimeria_table.csv","29_S25_eimeria_table.csv")

path_to_chicken_files <- "results/htseq/reverse/in_vitro_fusion/processed_reads/chicken"
path_to_eimeria_files <- "results/htseq/reverse/in_vitro_fusion/processed_reads/eimeria"

chicken_labels <- c("1_S17","2_S18","3_S19","4_S20","5_S21","6_S22","7_S23","8_S24","9_S25","10_S26","1_S1","2_S2","3_S3","4_S4",
                   "5_S5","6_S6","7_S7","8_S8","9_S9","11_S10","12_S11","13_S12","14_S13","15_S14","18_S15","19_S16","20_S17","21_S18",
                   "23_S19","24_S20","25_S21","26_S22","27_S23","28_S24","29_S25","30_S26","31_S27","32_S28")
eimeria_labels <- c("33_S29","2_S18","4_S20","6_S22","8_S24","10_S26","2_S2","4_S4","6_S6","8_S8","12_S11","14_S13","18_S15","20_S17",
                    "23_S19","24_S20","27_S23","28_S24","29_S25")

# Define the design parameters, here infection status (I=Infected, U=Uninfected) and time since infection
infection_status <- c("U","I","U","I","U","I","U","I","U","I","U","I","U","I","U","I","U","I","U","U","I","U","I","U",
                      "I","U","I","U","I","I","U","U","I","I","I","U","U","U")
timepoint_chicken <- c("4","4","24","24","2","2","4","4","24","24","0","4","4","24","24","48","48","72","72","0","48","48",
                       "72","72","48","48","72","72","2","2","2","2","12","12","12","12","12","12")
timepoint_eimeria <- c("0","4","24","2","4","24","4","24","48","72","48","72","48","72","2","2","12","12","12")
#batch_chicken <- c("160308","160308","160308","160308","160314","160314","160314","160314","160314","160314",
#                   "170131","170131","170131","170131","170131","170131","170131","170131","170131","170207",
#                   "170207","170207","170207","170207","170221","170221","170221","170221","170324","170324",
#                   "170324","170324","170410","170410","170410","170410","170410","170410")
#batch_eimeria <- c("170427","160308","160308","160314","160314","160314","170131","170131","170131","170131","170207",
#                   "170207","170221","170221","170324","170324","170410","170410","170410")
group_chicken <- factor(paste(infection_status,timepoint_chicken,sep="."))

# Create the DGElist object containing the count data in a format that edgeR can work with
chicken_dgelist <- readDGE(chicken_files, path = path_to_chicken_files, labels = chicken_labels, group = group_chicken, sep = ",")
eimeria_dgelist <- readDGE(eimeria_files, path = path_to_eimeria_files, labels = eimeria_labels, group = timepoint_eimeria, sep = ",")

# Add annotation information for the chicken data
pre_gene_list <- rownames(chicken_dgelist$counts)
gene_list <- sapply(pre_gene_list, function(x) substr(x, 6, nchar(x)))

idfound <- gene_list %in% mappedRkeys(org.Gg.egSYMBOL)
gene_list <- gene_list[idfound]

egSYMBOL <- toTable(org.Gg.egSYMBOL)

m <- match(gene_list, egSYMBOL$symbol)

entrez_geneIDs <- egSYMBOL$gene_id[m]
annotation_df <- data.frame(gene_list, egSYMBOL$gene_id[m])
colnames(annotation_df) <- c("Symbol","EntrezGene")

#Add the annotation information to the chicken_dgelist and filter out the genes that lacked annotation
#chicken_dgelist$genes <- annotation_df
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

pre_gene_list <- rownames(chicken_dgelist_filt_norm$counts)
gene_list <- sapply(pre_gene_list, function(x) substr(x, 6, nchar(x)))

m <- match(gene_list, egSYMBOL$symbol)

entrez_geneIDs <- egSYMBOL$gene_id[m]

# Check MDS plots for outliers or other potential issues
plotMDS(chicken_dgelist_filt_norm_genefilt)
plotMDS(eimeria_dgelist_filt_norm)

# The chicken analysis design
chicken_dgelist_filt_norm$samples$group <- relevel(chicken_dgelist_filt_norm$samples$group, ref="U.0")
design_chicken <- model.matrix(~group, data = chicken_dgelist_filt_norm$samples)
#design_chicken <- model.matrix(~infection_status+timepoint_chicken)
rownames(design_chicken) <- colnames(chicken_dgelist_filt_norm)

# Define the design of the eimeria analysis
design_eimeria <- model.matrix(~timepoint_eimeria)
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

# Test for differential expression in Eimeria
fit_eimeria <- glmQLFit(eimeria_disp, design_eimeria)

qlf_eimeria_12 <- glmQLFTest(fit_eimeria, coef = 2)
qlf_eimeria_2 <- glmQLFTest(fit_eimeria, coef = 3)
qlf_eimeria_24 <- glmQLFTest(fit_eimeria, coef = 4)
qlf_eimeria_4 <- glmQLFTest(fit_eimeria, coef = 5)
qlf_eimeria_48 <- glmQLFTest(fit_eimeria, coef = 6)
qlf_eimeria_72 <- glmQLFTest(fit_eimeria, coef = 7)
qlf_eimeria_all <- glmQLFTest(fit_eimeria, coef = 2:7)
#qlf_eimeria_batch <- glmQLFTest(fit_eimeria, coef = 8:14)

# Volcano plots for all comparisons
EnhancedVolcano(qlf_chicken_UI.2$table,
                lab = rownames(qlf_chicken_UI.2$table),
                title = 'Chicken in vitro 2h',
                x = 'logFC',
                y = 'PValue')
EnhancedVolcano(qlf_chicken_UI.4$table,
                lab = rownames(qlf_chicken_UI.4$table),
                title = 'Chicken in vitro 4h',
                x = 'logFC',
                y = 'PValue')
EnhancedVolcano(qlf_chicken_UI.12$table,
                lab = rownames(qlf_chicken_UI.12$table),
                title = 'Chicken in vitro 12h',
                x = 'logFC',
                y = 'PValue')
EnhancedVolcano(qlf_chicken_UI.24$table,
                lab = rownames(qlf_chicken_UI.24$table),
                title = 'Chicken in vitro 24h',
                x = 'logFC',
                y = 'PValue')
EnhancedVolcano(qlf_chicken_UI.48$table,
                lab = rownames(qlf_chicken_UI.48$table),
                title = 'Chicken in vitro 48h',
                x = 'logFC',
                y = 'PValue')
EnhancedVolcano(qlf_chicken_UI.72$table,
                lab = rownames(qlf_chicken_UI.72$table),
                title = 'Chicken in vitro 72h',
                x = 'logFC',
                y = 'PValue')

EnhancedVolcano(qlf_eimeria_2$table,
                lab = rownames(qlf_eimeria_2$table),
                title = 'Eimeria in vitro 2h',
                x = 'logFC',
                y = 'PValue')
EnhancedVolcano(qlf_eimeria_4$table,
                lab = rownames(qlf_eimeria_4$table),
                title = 'Eimeria in vitro 4h',
                x = 'logFC',
                y = 'PValue')
EnhancedVolcano(qlf_eimeria_12$table,
                lab = rownames(qlf_eimeria_12$table),
                title = 'Eimeria in vitro 12h',
                x = 'logFC',
                y = 'PValue')
EnhancedVolcano(qlf_eimeria_24$table,
                lab = rownames(qlf_eimeria_24$table),
                title = 'Eimeria in vitro 24h',
                x = 'logFC',
                y = 'PValue')
EnhancedVolcano(qlf_eimeria_48$table,
                lab = rownames(qlf_eimeria_48$table),
                title = 'Eimeria in vitro 48h',
                x = 'logFC',
                y = 'PValue')
EnhancedVolcano(qlf_eimeria_72$table,
                lab = rownames(qlf_eimeria_72$table),
                title = 'Eimeria in vitro 72h',
                x = 'logFC',
                y = 'PValue')
  
go_chicken_UI.2 <- goana(qlf_chicken_UI.2, geneid = entrez_geneIDs, FDR = 0.05, species = "Gg", species.KEGG = "gga")
go_chicken_UI.4 <- goana(qlf_chicken_UI.4, geneid = entrez_geneIDs, FDR = 0.05, species = "Gg", species.KEGG = "gga")
go_chicken_UI.12 <- goana(qlf_chicken_UI.12, geneid = entrez_geneIDs, FDR = 0.05, species = "Gg", species.KEGG = "gga")
go_chicken_UI.24 <- goana(qlf_chicken_UI.24, geneid = entrez_geneIDs, FDR = 0.05, species = "Gg", species.KEGG = "gga")
go_chicken_UI.48 <- goana(qlf_chicken_UI.48, geneid = entrez_geneIDs, FDR = 0.05, species = "Gg", species.KEGG = "gga")
go_chicken_UI.72 <- goana(qlf_chicken_UI.72, geneid = entrez_geneIDs, FDR = 0.05, species = "Gg", species.KEGG = "gga")

topgo_chicken_UI.2_up <- topGO(go_chicken_UI.2, ont="BP", sort="Up", n=30, truncate=30)
topgo_chicken_UI.4_up <- topGO(go_chicken_UI.4, ont="BP", sort="Up", n=30, truncate=30)
topgo_chicken_UI.12_up <- topGO(go_chicken_UI.12, ont="BP", sort="Up", n=30, truncate=30)
topgo_chicken_UI.24_up <- topGO(go_chicken_UI.24, ont="BP", sort="Up", n=30, truncate=30)
topgo_chicken_UI.48_up <- topGO(go_chicken_UI.48, ont="BP", sort="Up", n=30, truncate=30)
topgo_chicken_UI.72_up <- topGO(go_chicken_UI.72, ont="BP", sort="Up", n=30, truncate=30)

topgo_chicken_UI.2_down <- topGO(go_chicken_UI.2, ont="BP", sort="Down", n=30, truncate=30)
topgo_chicken_UI.4_down <- topGO(go_chicken_UI.4, ont="BP", sort="Down", n=30, truncate=30)
topgo_chicken_UI.12_down <- topGO(go_chicken_UI.12, ont="BP", sort="Down", n=30, truncate=30)
topgo_chicken_UI.24_down <- topGO(go_chicken_UI.24, ont="BP", sort="Down", n=30, truncate=30)
topgo_chicken_UI.48_down <- topGO(go_chicken_UI.48, ont="BP", sort="Down", n=30, truncate=30)
topgo_chicken_UI.72_down <- topGO(go_chicken_UI.72, ont="BP", sort="Down", n=30, truncate=30)















