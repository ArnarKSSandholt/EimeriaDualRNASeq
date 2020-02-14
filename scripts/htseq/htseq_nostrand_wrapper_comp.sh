#!/bin/bash -l

#SBATCH -A snic2020-15-16
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 48:00:00
#SBATCH -J htseq_run_nostrand_comp
#SBATCH --mail-type=ALL
#SBATCH --mail-user arnarkari.sigurearsonsandholt.9531@student.uu.se

module load bioinfo-tools
module load htseq/0.9.1

# Your commands
bash scripts/htseq/htseq_nostrand_run.sh results/star/mapped_reads/in_vitro_complementary