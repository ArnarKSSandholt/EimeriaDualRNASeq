#!/bin/bash -l

#SBATCH -A snic2020-15-16
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 48:00:00
#SBATCH -J trimmomatic_run
#SBATCH --mail-type=ALL
#SBATCH --mail-user arnarkari.sigurearsonsandholt.9531@student.uu.se

module load bioinfo-tools
module load trimmomatic

# Your commands
bash scripts/trimmomatic/trimmomatic_run.sh data/in_vitro_pilot/161012_D00457_0163_AC9TWTANXX/Sample_*
bash scripts/trimmomatic/trimmomatic_run.sh data/in_vitro_complementary/170830_D00457_0216_ACB7DNANXX/Sample_*
bash scripts/trimmomatic/trimmomatic_run.sh data/SH-2259/191210_A00605_0093_BHNFNHDSXX/Sample_SH-2259-*
bash scripts/trimmomatic/trimmomatic_run.sh data/SI-2311/191105_A00605_0082_AHGV3LDRXX/Sample_SI-2311-*