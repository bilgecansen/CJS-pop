#!/bin/bash
#
#SBATCH --job-name=brcr_models
#SBATCH --output=tempout_brcr.txt
#SBATCH --ntasks-per-node=28
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH -p long-28core

module load shared
module load R/3.6.0
module load gcc-stack
module load openblas/dynamic/0.2.18
module load lapack/gcc/64/3.6.0
module load JAGS/4.3.0
export R_LIBS=/gpfs/home/bsen/R_libs

cd /gpfs/home/bsen/CJS-pop
Rscript brcr_models.R