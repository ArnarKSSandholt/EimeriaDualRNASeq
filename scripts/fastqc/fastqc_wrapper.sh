#!/bin/bash -l

#SBATCH -A snic2020-15-16
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH -J fastqc_run
#SBATCH --mail-type=ALL
#SBATCH --mail-user arnarkari.sigurearsonsandholt.9531@student.uu.se

module load bioinfo-tools
module load FastQC/0.11.8
module load MultiQC/1.8

# Your commands
bash scripts/fastqc/fastqc_run.sh results/trimmomatic/in_vitro_pilot
bash scripts/fastqc/fastqc_run.sh results/trimmomatic/in_vitro_complementary
bash scripts/fastqc/fastqc_run.sh results/trimmomatic/SH-2259
bash scripts/fastqc/fastqc_run.sh results/trimmomatic/SI-2311