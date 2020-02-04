#!/bin/bash -l

#SBATCH -A snic2020-15-16
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 2:00:00
#SBATCH -J trimmomatic_run
#SBATCH --mail-type=ALL
#SBATCH --mail-user arnarkari.sigurearsonsandholt.9531@student.uu.se

module load bioinfo-tools
module load FastQC/0.11.8
module load MultiQC/1.8

# Your commands
bash scripts/fastqc/fastqc_run.sh results/trimmomatic/in_vitro_pilot