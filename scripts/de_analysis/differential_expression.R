# A pipeline for the differential expression analysis of a dual RNA-seq dataset of chickens and Eimeria tenella parsites

library(edgeR)

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
batch_chicken <- c("160308","160308","160308","160308","160314","160314","160314","160314","160314","160314",
                   "170131","170131","170131","170131","170131","170131","170131","170131","170131","170207",
                   "170207","170207","170207","170207","170221","170221","170221","170221","170324","170324",
                   "170324","170324","170410","170410","170410","170410","170410","170410")
batch_eimeria <- c("160308","160308","160314","160314","160314","170131","170131","170131","170131","170207",
                   "170207","170221","170221","170324","170324","170410","170410","170410","170427")
group_chicken <- factor(paste(infection_status,timepoint_chicken,sep="."))

# Create the DGElist object containing the count data in a format that edgeR can work with
chicken_dgelist <- readDGE(chicken_files, path = path_to_chicken_files, labels = chicken_labels, group = group_chicken, sep = ",")
eimeria_dgelist <- readDGE(eimeria_files, path = path_to_eimeria_files, labels = eimeria_labels, group = timepoint_eimeria, sep = ",")

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

# TODO: The chicken analysis design
#data.frame(Sample=colnames(chicken_dgelist),infection_status,timepoint,group_chicken)
#design <- model.matrix(~0+group_chicken)
#rownames(design) <- colnames(chicken_dgelist)

# Define the design of the eimeria analysis
design_eimeria <- model.matrix(~timepoint_eimeria+batch_eimeria)
rownames(design_eimeria) <- colnames(eimeria_dgelist_filt_norm)

# Estimate the dispersion for each dataset and plot for the biological coefficient of variation (BCV)
#chicken_disp <- estimateDisp(chicken_dgelist, design, robust = TRUE)
eimeria_disp <- estimateDisp(eimeria_dgelist_filt_norm, design_eimeria, robust = TRUE)

#plotBCV(chicken_disp)
plotBCV(eimeria_disp)

# Test for differential expression in Eimeria
fit_eimeria <- glmQLFit(eimeria_disp, design_eimeria)

qlf_eimeria_12 <- glmQLFTest(fit_eimeria, coef = 2)
qlf_eimeria_2 <- glmQLFTest(fit_eimeria, coef = 3)
qlf_eimeria_24 <- glmQLFTest(fit_eimeria, coef = 4)
qlf_eimeria_4 <- glmQLFTest(fit_eimeria, coef = 5)
qlf_eimeria_48 <- glmQLFTest(fit_eimeria, coef = 6)
qlf_eimeria_72 <- glmQLFTest(fit_eimeria, coef = 7)
qlf_eimeria_all <- glmQLFTest(fit_eimeria, coef = 2:7)
qlf_eimeria_batch <- glmQLFTest(fit_eimeria, coef = 8:14)



