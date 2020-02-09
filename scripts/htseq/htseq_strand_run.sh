# Script for running htseq-count on a set of mapped read files in bam format
# Usage: bash htseq_strand_run.sh /Path/to/mapped/read/file/folder

GROUP_NAME=$(echo ${1} | cut -d '/' -f 4)
mkdir -p results/htseq/${GROUP_NAME}

for f in $(ls ${1}/*.bam | sed 's/_Aligned.sortedByCoord.out.bam//' | sort -u)
do
    OUTPUT_NAME=$(echo ${f} | rev | cut -d '/' -f 1 | rev)
    htseq-count -f bam -r pos -s yes -i Dbxref ${f}_Aligned.sortedByCoord.out.bam results/star/merged_reference/eimeria_chicken_merge.gff \
        > results/htseq/${GROUP_NAME}/${OUTPUT_NAME}_stranded.counts
done