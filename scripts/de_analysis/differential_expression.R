# A pipeline for the differential expression analysis of a dual RNA-seq dataset of chickens and Eimeria tenella parsites

library(edgeR)

pilot_chicken_files <- c("1_S17_chicken_table.csv","2_S18_chicken_table.csv","3_S19_chicken_table.csv","4_S20_chicken_table.csv",
                         "5_S21_chicken_table.csv","6_S22_chicken_table.csv","7_S23_chicken_table.csv","8_S24_chicken_table.csv",
                         "9_S25_chicken_table.csv","10_S26_chicken_table.csv")
pilot_eimeria_files <- c("1_S17_eimeria_table.csv","2_S18_eimeria_table.csv","3_S19_eimeria_table.csv","4_S20_eimeria_table.csv",
                         "5_S21_eimeria_table.csv","6_S22_eimeria_table.csv","7_S23_eimeria_table.csv","8_S24_eimeria_table.csv",
                         "9_S25_eimeria_table.csv","10_S26_eimeria_table.csv")
path_to_files <- "results/htseq/reverse/in_vitro_pilot/processed_reads"
sample_labels <- c("1_S17","2_S18","3_S19","4_S20","5_S21","6_S22","7_S23","8_S24","9_S25","10_S26")

pilot_chicken_merge <- readDGE(pilot_chicken_files, path = path_to_files, labels = sample_labels, sep = ",")
pilot_eimeria_merge <- readDGE(pilot_eimeria_files, path = path_to_files, labels = sample_labels, sep = ",")
