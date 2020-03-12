# Concatenates the count files for each organism into a single file and fixes the gene names for both
# Usage: python read_count_concat.py /Path/to/chicken/read/files /Path/to/Eimeria/read/files /Path/to/Eimeria/gene/list.tsv /Path/to/output/files

import pandas as pd
from os import listdir, mkdir
import sys
import re

chicken_path = sys.argv[1]
eimeria_path = sys.argv[2]
eimeria_gene_path = sys.argv[3]
output_path = sys.argv[4]

chicken_file_names = listdir(chicken_path)
eimeria_file_names = listdir(eimeria_path)
first_table = True

# Concatenate the chicken files and remove the gene- prefix in the gene names
for file_name in chicken_file_names:
    sample_name = file_name[0:re.search("_chicken", file_name).start()]
    current_table = pd.read_csv(chicken_path+"/"+file_name, header = 0, names = ["gene_name",sample_name])
    if first_table:
        old_table = current_table
        first_table = False
    else:
        old_table = pd.merge(old_table, current_table, on = "gene_name")

i = 0
while i < len(old_table):
    old_table.iloc[i,0] = old_table.iloc[i,0][5:]
    i += 1

old_table.to_csv(output_path+"/chicken_counts.csv", sep = ",", index = False)

# Concatenate the Eimeria files and replace the gene identifiers with proper ones
first_table = True

for file_name in eimeria_file_names:
    sample_name = file_name[0:re.search("_eimeria", file_name).start()]
    current_table = pd.read_csv(eimeria_path+"/"+file_name, header = 0, names = ["gene_name",sample_name])
    if first_table:
        old_table = current_table
        first_table = False
    else:
        old_table = pd.merge(old_table, current_table, on = "gene_name")

gene_table = pd.read_csv(eimeria_gene_path, sep = "\t", header = None, names = ["gene_name","gene_id"])
old_table = pd.merge(old_table, gene_table, on = "gene_name")
old_table = pd.concat([old_table.iloc[:,-1], old_table.iloc[:,1:-1]], axis = 1)
old_table.to_csv(output_path+"/eimeria_counts.csv", sep = ",", index = False)