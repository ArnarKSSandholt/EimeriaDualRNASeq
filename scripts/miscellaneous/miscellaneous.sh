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

# Command used to trim the short reads for Eimeria tenella genome assembly
trimmomatic PE -threads 4 \
data/assembly_reads/short_reads/Eimeria_S2_L001_R1_001.fastq.gz \
data/assembly_reads/short_reads/Eimeria_S2_L001_R2_001.fastq.gz \
results/assembly/trimmomatic/Eimeria_R1_paired.fastq.gz \
results/assembly/trimmomatic/Eimeria_R1_unpaired.fastq.gz \
results/assembly/trimmomatic/Eimeria_R2_paired.fastq.gz \
results/assembly/trimmomatic/Eimeria_R2_unpaired.fastq.gz \
ILLUMINACLIP:data/adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads \
SLIDINGWINDOW:4:20