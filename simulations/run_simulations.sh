#!/bin/bash -l

#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nathan.wisnoski@msstate.edu
#SBATCH --job-name=metacom_dispersal_kernel

cd /mnt/scratch/wisnoskilab/GitHub/metacom-dispersal

module load R/4.2.2

R CMD BATCH --no-restore --no-save simulations/metacom_diversity_variability.R
