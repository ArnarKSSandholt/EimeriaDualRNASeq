#!/bin/bash -l

#SBATCH -A snic2020-15-16
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 10:00:00
#SBATCH -J trimmomatic_run
#SBATCH --mail-type=ALL
#SBATCH --mail-user arnarkari.sigurearsonsandholt.9531@student.uu.se

module load conda
module load bioinfo-tools
module load snakemake

# Your commands
source conda_init.sh

snakemake --cores 8 --use-conda results/trimmomatic/SH-2259/SH-2259-11_S5_L001_R1_paired.fastq.gz
snakemake --cores 8 --use-conda results/trimmomatic/SI-2311/SI-2311-CT11_S1_L001_R1_001.fastq.gz