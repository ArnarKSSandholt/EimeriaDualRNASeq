# Script for running the STAR aligner for a set of read files and reference genomes
# Usage: bash star_run.sh /Path/to/read/files/folder

GROUP_NAME=$(echo ${1} | cut -d '/' -f 3)
mkdir -p results/star/mapped_reads/${GROUP_NAME}

for f in $(ls ${1}/*_paired.fastq.gz | sed 's/${1}//' | sed 's/_R._paired.fastq.gz//' | sort -u)
do
    echo ${f}
    #STAR --genomeDir /proj/snic2020-16-20/EimeriaDualRNASeq/results/star/index --readFilesIn ${1}/${f}_R1_paired.fastq.gz ${1}/${f}_R2_paired.fastq.gz --readFilesCommand gunzip -c \
		#--runThreadN 8 --outFileNamePrefix results/star/mapped_reads/${GROUP_NAME}/${f}_ --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate
    #STAR --genomeDir /proj/snic2020-16-20/EimeriaDualRNASeq/results/star/index --readFilesIn ${1}/${f}_R1_unpaired.fastq.gz --readFilesCommand gunzip -c \
		#--runThreadN 8 --outFileNamePrefix results/star/mapped_reads/${GROUP_NAME}/${f}_R1_unpaired_ --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate
    #STAR --genomeDir /proj/snic2020-16-20/EimeriaDualRNASeq/results/star/index --readFilesIn ${1}/${f}_R2_unpaired.fastq.gz --readFilesCommand gunzip -c \
		#--runThreadN 8 --outFileNamePrefix results/star/mapped_reads/${GROUP_NAME}/${f}_R2_unpaired_ --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate
done


