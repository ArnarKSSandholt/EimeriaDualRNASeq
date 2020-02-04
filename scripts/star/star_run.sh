# Script for running the STAR aligner for a set of read files and reference genomes
# Usage: bash star_run.sh /Path/to/read/files/folder

GROUP_NAME=$(echo ${1} | cut -d '/' -f 3)
mkdir -p results/star/mapped_reads/${GROUP_NAME}

for f in $(ls ${1}/*_paired.fastq.gz | sed 's/_R._paired.fastq.gz//' | sort -u)
do
    STAR --genomeDir results/star/index --readFilesIn ${1}/${f}_R1_paired.fastq.gz ${1}/${f}_R2_paired.fastq.gz --readFilesCommand gunzip -c \
		--runThreadN 16 --outFileNamePrefix results/star/mapped_reads/${GROUP_NAME}/${f}_ --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate
done


