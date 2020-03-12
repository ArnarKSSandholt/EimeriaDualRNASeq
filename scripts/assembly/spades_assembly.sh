#!/bin/bash -l

#SBATCH -A snic2020-15-16
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 48:00:00
#SBATCH -J eimeria_spades_assembly
#SBATCH --mail-type=ALL
#SBATCH --mail-user arnarkari.sigurearsonsandholt.9531@student.uu.se

module load bioinfo-tools
module load spades

# Your commands
spades.py -o results/assembly/spades \
    -1 results/assembly/trimmomatic/Eimeria_R1_paired.fastq.gz \
    -2 results/assembly/trimmomatic/Eimeria_R1_paired.fastq.gz
