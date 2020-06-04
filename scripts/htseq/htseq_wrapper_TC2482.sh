#!/bin/bash -l

#SBATCH -A snic2020-15-16
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 120:00:00
#SBATCH -J htseq_run_rev_TC-2482
#SBATCH --mail-type=ALL
#SBATCH --mail-user arnarkari.sigurearsonsandholt.9531@student.uu.se

module load bioinfo-tools
module load htseq/0.9.1

# Your commands
bash scripts/htseq/htseq_rev_run.sh results/star/mapped_reads/TC-2482