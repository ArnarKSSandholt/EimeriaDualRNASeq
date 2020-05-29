#!/bin/bash -l

#SBATCH -A snic2020-15-16
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 96:00:00
#SBATCH -J star_run_TC-2483
#SBATCH --mail-type=ALL
#SBATCH --mail-user arnarkari.sigurearsonsandholt.9531@student.uu.se

module load bioinfo-tools
module load star/2.7.2b

# Your commands
bash scripts/star/star_run.sh results/trimmomatic/TC-2483