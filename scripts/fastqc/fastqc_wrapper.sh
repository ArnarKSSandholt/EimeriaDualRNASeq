#!/bin/bash -l

#SBATCH -A snic2020-15-16
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 2:00:00
#SBATCH -J trimmomatic_run
#SBATCH --mail-type=ALL
#SBATCH --mail-user arnarkari.sigurearsonsandholt.9531@student.uu.se

module load conda
module load bioinfo-tools
module load snakemake

# Your commands
source conda_init.sh

snakemake --cores 8 --use-conda results/fastqc/in_vitro_pilot/1_S17_L002_R1.html