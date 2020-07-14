# Read counting pipeline

This document contains all code used for the read counting portion of the project. It is presented in the order of execution. Information on software used and versions can be found in the text and in the final section of the document.

## Trimming

The raw reads had been quality checked by the sequencing service, SNP&SEQ Technology Platform. Their reports were used as a basis to determine the settings required for trimming. The software used to trim the raw read data was Trimmomatic (v. 0.36). For both the *in vitro* and *in vivo* datasets, the data was of high quality but the reverse reads had a large amount of adapter content. <br>
As such, only a SLIDINGWINDOW:4:20 was used for quality, i.e. cutting reads anywhere Trimmomatic finds a window of 4 bases with a mean quality lower than 20, and removing reads shorter than 50 bp after trimming. The adapters used by the sequencing lab were standard Illumina adapters and were therefore included in the adapter library provided with Trimmomatic, allowing it to be used to identify and remove adapters. <br>
The script below was used with these settings to trim the raw reads from both the *in vitro* and *in vivo* datasets:

```bash
# Short bash script for using trimmomatic with user input for the location of the read files
# Usage: bash trimmomatic_run.sh /Path/to/folder/containing/read/files1 /Path/to/folder/containing/read/files2 ...

for FILE_PATH in "$@"
do
    GROUP_NAME=$(echo ${FILE_PATH} | cut -d '/' -f 2)
    mkdir -p results/trimmomatic/${GROUP_NAME}
    for f in $(ls ${FILE_PATH} | sed 's/_.._001.fastq.gz//' | sort -u)
    do
        trimmomatic PE -threads 8 \
        ${FILE_PATH}/${f}_R1_001.fastq.gz \
        ${FILE_PATH}/${f}_R2_001.fastq.gz \
        results/trimmomatic/${GROUP_NAME}/${f}_R1_paired.fastq.gz \
        results/trimmomatic/${GROUP_NAME}/${f}_R1_unpaired.fastq.gz \
        results/trimmomatic/${GROUP_NAME}/${f}_R2_paired.fastq.gz \
        results/trimmomatic/${GROUP_NAME}/${f}_R2_unpaired.fastq.gz \
        ILLUMINACLIP:data/adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads \
        SLIDINGWINDOW:4:20 \
        MINLEN:50
    done
done
```

_____________________________

## Quality checking

In order to check that the trimming had been successful, the reads were quality checked using FastQC (v. 0.11.8) and MultiQC (v. 1.8) was used to collate the results using the following script:

```bash
# Script for running FastQC on a set of Trimmomatic trimmed read files and outputting both FastQC and MultiQC reports for them
# Usage: bash fastqc_run.sh /Path/to/folder/containing/trimmed/read/files1 /Path/to/folder/containing/trimmed/read/files2 ...

for FILE_PATH in "$@"
do
    GROUP_NAME=$(echo ${FILE_PATH} | cut -d '/' -f 3)
    mkdir -p results/fastqc/${GROUP_NAME} 
    fastqc -o results/fastqc/${GROUP_NAME} ${FILE_PATH}/*.fastq.gz
    multiqc -o results/fastqc/${GROUP_NAME} results/fastqc/${GROUP_NAME}
done
```

___________________________

## Mapping

The mapping was the first stage of the analysis where the "dual RNA-seq" analysis differed from a regular RNA-seq analysis. The reads needed to be mapped to both genomes. This was done by concatenating both the fna sequence files and gtf annotation files for each organism, effectively treating them as a single genome. This was done using the following script:

```bash
# Script for merging two reference genomes for dual RNA-seq analysis.  Each reference should have a fna and gff file
# Usage: bash reference_merger.sh /Path/to/reference1.fna /Path/to/reference1.gff /Path/to/reference2.fna /Path/to/reference2.gff /Path/to/merged.fna /Path/to/merged.gff

cat ${1} ${3} > ${5}

touch ${6}
head -n -1 ${2} >> ${6}
tail -n +7 ${4} >> ${6}
```

The reads were then mapped to this pseudo-genome using STAR (v. 2.7.2b), a splice-aware mapper. To do so, the pseudo-genome first needed to be indexed using the following script:

```bash
# Script for producing an index for STAR
# Usage:  bash star_index.sh /Path/to/reference/genome.fna /Path/to/reference/genome.gff

mkdir -p results/star/index

STAR --runMode genomeGenerate --genomeDir results/star/index --genomeFastaFiles ${1} --sjdbGTFtagExonParentTranscript Parent \
	--sjdbGTFtagExonParentGene gene --runThreadN 8 --sjdbGTFfile ${2} --sjdbOverhang 150 --genomeChrBinNbits 18
```

Finally, the mapping was run on all trimmed reads using the following script:

```bash
# Script for running the STAR aligner for a set of read files and reference genomes
# Usage: bash star_run.sh /Path/to/read/files/folder

GROUP_NAME=$(echo ${1} | cut -d '/' -f 3)
mkdir -p results/star/mapped_reads/${GROUP_NAME}

for f in $(ls ${1}/*_paired.fastq.gz | sed 's/_R._paired.fastq.gz//' | sort -u)
do
    OUTPUT_NAME=$(echo ${f} | rev | cut -d '/' -f 1 | rev)
    STAR --genomeDir /proj/snic2020-16-20/EimeriaDualRNASeq/results/star/index --readFilesIn ${f}_R1_paired.fastq.gz ${f}_R2_paired.fastq.gz --readFilesCommand gunzip -c \
		--runThreadN 8 --outFileNamePrefix results/star/mapped_reads/${GROUP_NAME}/${OUTPUT_NAME}_ --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate
    STAR --genomeDir /proj/snic2020-16-20/EimeriaDualRNASeq/results/star/index --readFilesIn ${f}_R1_unpaired.fastq.gz --readFilesCommand gunzip -c \
		--runThreadN 8 --outFileNamePrefix results/star/mapped_reads/${GROUP_NAME}/${OUTPUT_NAME}_R1_unpaired_ --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate
    STAR --genomeDir /proj/snic2020-16-20/EimeriaDualRNASeq/results/star/index --readFilesIn ${f}_R2_unpaired.fastq.gz --readFilesCommand gunzip -c \
		--runThreadN 8 --outFileNamePrefix results/star/mapped_reads/${GROUP_NAME}/${OUTPUT_NAME}_R2_unpaired_ --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate
done
```

__________________________________

## Read counting

With the reads mapped to both genomes, the number of reads mapping to each gene needed to be counted.  This was done using HTSeq (v. 0.9.1).  As the read data was reverse stranded, the -s strandedness setting was set on reverse to make use of the extra information available with strandedness.  The following code was used to run HTSeq:

```bash
# Script for running htseq-count on a set of mapped read files in bam format
# Usage: bash htseq_reverse_run.sh /Path/to/mapped/read/file/folder

GROUP_NAME=$(echo ${1} | cut -d '/' -f 4)
mkdir -p results/htseq/reverse/${GROUP_NAME}

for f in $(ls ${1}/*.bam | sed 's/_Aligned.sortedByCoord.out.bam//' | sort -u)
do
    OUTPUT_NAME=$(echo ${f} | rev | cut -d '/' -f 1 | rev)
    htseq-count -f bam -r pos -s reverse ${f}_Aligned.sortedByCoord.out.bam results/star/merged_reference/eimeria_chicken_merge_gene_id.gtf \
        > results/htseq/reverse/${GROUP_NAME}/${OUTPUT_NAME}_reverse.counts
done
```
______________________

## Read count processing

Before the read counts could be analyzed for differential expression, they had to be processed to make them easier to work with.  The first step was to sum up the paired-end read counts with the read counts for reads that had lost their read pair in the trimming step.  The statistics at the bottom of the HTSeq output files were also extracted and formatted in a more readable way.

```python
# Script that processes the read count data from htseq-count and splits the read counts by source as well
# as computing statistics about them, such as how large a fraction of the mapped reads map to the parasite
# Usage: python read_count_process.py /Path/to/metadata/file.tsv /Path/to/read/count/folder /Path/to/output/folder

import pandas as pd
from os import listdir, mkdir
import sys
import re

# Read in the input and output paths from stdin and define folders for the output
metadata_path = sys.argv[1]
read_count_path = sys.argv[2]
output_path = sys.argv[3]
mkdir(output_path+"/chicken")
mkdir(output_path+"/eimeria")
mkdir(output_path+"/stats")

# Read in the data
read_count_filenames = listdir(read_count_path)
read_count_filenames.sort()
metadata_table = pd.read_csv(metadata_path, sep="\t")
metadata_table = metadata_table.sort_values("File_name")
result_list = []
old_name = ""

i = 0
j = 0
while i < len(read_count_filenames):
    # The loop goes through each read count file, separates the chicken and E. tenella reads and sums the 
    # read counts from paired and unpaired reads for a given sample.  It also extracts statistics from the
    # bottom of the HTSeq output file
    if re.search("L00", read_count_filenames[i]):
        # Make sure only read count files are read in
        l_pos = re.search("L00", read_count_filenames[i]).start()
        curr_name = read_count_filenames[i][0:l_pos-1]

        if curr_name != old_name:
            sum_table = pd.read_csv(read_count_path+"/"+read_count_filenames[i], sep="\t", header=None, names=["gene_name","gene_count"], index_col=0)
        else:
            temp_table = pd.read_csv(read_count_path+"/"+read_count_filenames[i], sep="\t", header=None, names=["gene_name","gene_count"], index_col=0)
            sum_table += temp_table

            if i+1 == len(read_count_filenames) or not re.search(curr_name, read_count_filenames[i+1]):
                stat_table = sum_table.iloc[-5:]
                eimeria_table = sum_table.iloc[-8660:-5]
                chicken_table = sum_table.iloc[:-8660]

                total_read_num_sum = sum(sum_table.gene_count.iloc[:-5])
                eimeria_read_num = sum(eimeria_table.gene_count)
                chicken_read_num = sum(chicken_table.gene_count)

                eimeria_read_perc = eimeria_read_num/total_read_num_sum * 100
                no_feature_num = stat_table.gene_count.iloc[0]
                ambiguous_num = stat_table.gene_count.iloc[1]

                result_list.append([curr_name, total_read_num_sum, chicken_read_num, eimeria_read_num, eimeria_read_perc, no_feature_num, ambiguous_num])
                
                stat_table.to_csv(output_path+"/stats/"+curr_name+"_stat_table.csv", sep = ",")
                eimeria_table.to_csv(output_path+"/eimeria/"+curr_name+"_eimeria_table.csv", sep = ",")
                chicken_table.to_csv(output_path+"/chicken/"+curr_name+"_chicken_table.csv", sep = ",")
                j += 1
        old_name = curr_name
    i += 1

header_list = ["File_name","Total_number_of_mapped_reads", "Number_of_reads_mapped_to_chicken", "Number_of_reads_mapped_to_Eimeria", \
    "Percentage_of_Eimeria_reads", "Number_of_read_not_mapped_to_feature", "Number_of_reads_mapped_to_multiple_features"]
data_table = pd.DataFrame(result_list, columns=header_list)
out_table = pd.merge(metadata_table, data_table, on = "File_name")
out_table.to_csv(output_path+"/metadata_table.csv", sep = ",", index = False)
```

The next step was to concatenate the read counts from different samples to produce a single file containing all the read counts for each species and experiment.  The data from the two *in vitro* experiments was combined as the analysis was meant to be conducted on both together.  The gene identifiers also needed to be replaced with Entrez gene identifiers and locus tags to facilitate the use of certain packages in the downstream analysis.

```python
# Concatenates the count files for each organism into a single file and fixes the gene names for both.  Also fixes the
# gene names of an Eimeria annotation file.
# Usage: python read_count_concat.py /Path/to/chicken/read/files /Path/to/Eimeria/read/files /Path/to/chicken/gene/list.tsv \
# /Path/to/Eimeria/gene/list.tsv /Path/to/Eimeria/product/list /Path/to/output/files

import pandas as pd
from os import listdir, mkdir
import sys
import re

chicken_path = sys.argv[1]
eimeria_path = sys.argv[2]
chicken_gene_path = sys.argv[3]
eimeria_gene_path = sys.argv[4]
eimeria_product_path = sys.argv[5]
output_path = sys.argv[6]

chicken_file_names = listdir(chicken_path)
eimeria_file_names = listdir(eimeria_path)
first_table = True

# Concatenate the chicken files, remove the gene- prefix in the gene symbols and replace the gene identifiers with Entrez Gene IDs
for file_name in chicken_file_names:
    sample_name = file_name[0:re.search("_chicken", file_name).start()]
    current_table = pd.read_csv(chicken_path+"/"+file_name, header = 0, names = ["gene_name",sample_name])
    if first_table:
        old_table = current_table
        first_table = False
    else:
        old_table = pd.merge(old_table, current_table, on = "gene_name")

chicken_gene_table = pd.read_csv(chicken_gene_path, sep = "\t", header = None, names = ["gene_name","entrez_gene_id"])
old_table = pd.merge(old_table, chicken_gene_table, on = "gene_name")
temp_table = pd.concat([old_table.iloc[:,-1], old_table.iloc[:,1:-1]], axis = 1)
old_table = pd.concat([old_table.iloc[:,0], temp_table], axis = 1)

i = 0
while i < len(old_table):
    old_table.iloc[i,0] = old_table.iloc[i,0][5:]
    i += 1

old_table.to_csv(output_path+"/chicken_counts.csv", sep = ",", index = False)

# Concatenate the Eimeria files and replace the gene identifiers with Entrez Gene IDs and locus tags
first_table = True

for file_name in eimeria_file_names:
    sample_name = file_name[0:re.search("_eimeria", file_name).start()]
    current_table = pd.read_csv(eimeria_path+"/"+file_name, header = 0, names = ["gene_name",sample_name])
    if first_table:
        old_table = current_table
        first_table = False
    else:
        old_table = pd.merge(old_table, current_table, on = "gene_name")

eimeria_gene_table = pd.read_csv(eimeria_gene_path, sep = "\t", header = None, names = ["gene_name","entrez_gene_id","locus_tag"])
old_table = pd.merge(old_table, eimeria_gene_table, on = "gene_name")
old_table = pd.concat([old_table.iloc[:,-1], old_table.iloc[:,-2], old_table.iloc[:,1:-2]], axis = 1)
old_table.to_csv(output_path+"/eimeria_counts.csv", sep = ",", index = False)

# Replace Eimeria gene identifiers with Entrez Gene IDs and locus tags for an annotation file
eimeria_product_table = pd.read_csv(eimeria_product_path, sep = "\t", header = None, names = ["gene_name","product"])
eimeria_product_table = pd.merge(eimeria_product_table, eimeria_gene_table, on = "gene_name")
eimeria_product_table = pd.concat([eimeria_product_table.iloc[:,-1], eimeria_product_table.iloc[:,-2], eimeria_product_table.iloc[:,1]], axis = 1)
eimeria_product_table.to_csv(output_path+"/eimeria_gene_products.csv", sep = "\t", index = False)
```

Some bash commands were also used to extract mappings of gene identifiers to each other and to annotation information:

```bash
# Commanda used to create the eimeria and chicken gene name lists from the eimeria and chicken gff files from NCBI
cat GCF_000499545.2_ETH001_genomic_eimeria.gff | tail -n +9 | grep gbkey=Gene | cut -f 9 | cut -d ';' -f 1,2,3 | \
sed -e 's/ID=//' -e 's/;Dbxref=GeneID:/\t/' -e 's/;Name=/\t/' > eimeria_genes.tsv

cat GCF_000002315.6_GRCg6a_genomic_chicken.gff | tail -n +9 | grep gbkey=Gene | cut -f 9 | cut -d ';' -f 1,2 | \
sed -r -e 's/ID=//' -e 's/Dbxref=GeneID://' -e 's/Dbxref=CGNC:.+,GeneID://' -e 's/,miRBase.+//' -e 's/;/\t/' > chicken_genes.tsv

# Command used to do the same as above for matching transcript identifiers for a KEGG analysis with gene identifiers
cat GCF_000499545.2_ETH001_genomic_eimeria.gff | tail -n +9 | grep gbkey=CDS | cut -f 9 | cut -d ';' -f 3 | \
sed -e 's/Dbxref=GeneID://' -e 's/,Genbank:/\t/' -e 's/Dbxref=Genbank://' -e 's/,GeneID:/\t/' | sort -u > eimeria_transcript_to_gene.tsv

# Commands used to extract gene identifiers and product names from the gff annotation files
cat GCF_000499545.2_ETH001_genomic_eimeria.gff | tail -n +9 | grep gbkey=mRNA | grep ID=rna | sed 's/exception=low-quality sequence region;//' | \
cut -f 9 | cut -d ";" -f 2,8 | sed -e 's/Parent=//' -e 's/;product=/\t/' > eimeria_gene_products.tsv

cat GCF_000002315.6_GRCg6a_genomic_chicken.gff | tail -n +9 | grep ID=rna | grep mRNA | cut -f 9 | \
sed -e 's/;end_range//' -e 's/partial=true;//' -e 's/;Note=//' -e 's/;exception=//' -e 's/;model_evidence=//' -e 's/;inference=//' | \
cut -d ";" -f 2,7 | sed -e 's/Parent=gene-//' -e 's/;product=/\t/' > chicken_gene_products.tsv
```

With this, the data was ready to be analyzed in R using the edgeR package for differential expression analysis.  This is the subject of the second markdown file.

____________________________________

## Software source and versions

Trimmomatic v. 0.36

FastQC v. 0.11.8

MultiQC v. 1.8

STAR v. 2.7.2b

HTSeq v. 0.9.1
