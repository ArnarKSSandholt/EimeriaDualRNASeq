#!/bin/bash -l

#SBATCH -A snic2020-15-16
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 2:00:00
#SBATCH -J trimmomatic_run
#SBATCH --mail-type=ALL
#SBATCH --mail-user arnarkari.sigurearsonsandholt.9531@student.uu.se

module load bioinfo-tools
module load trimmomatic

# Your commands
bash scripts/trimmomatic/trimmomatic_run.sh data/in_vitro_pilot/161012_D00457_0163_AC9TWTANXX/Sample_1