# Script for running the STAR aligner for a set of read files and reference genomes


for f in $(ls ${1} | sed 's/_R._paired.fastq.gz//' | sort -u)
do
    STAR --genomeDir results/star/index --readFilesIn ${1}/${f}_R1_paired.fastq.gz ${1}/${f}_R2_paired.fastq.gz --readFilesCommand gunzip -c \
		--runThreadN 16 --outFileNamePrefix results/star/_ --quantMode GeneCounts
done


