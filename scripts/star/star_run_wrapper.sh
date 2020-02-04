#!/bin/bash -l

#SBATCH -A snic2020-15-16
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 1:00:00
#SBATCH -J trimmomatic_run
#SBATCH --mail-type=ALL
#SBATCH --mail-user arnarkari.sigurearsonsandholt.9531@student.uu.se

module load bioinfo-tools
module load star/2.7.2b

# Your commands
bash scripts/star/star_run.sh results/trimmomatic/in_vitro_pilot