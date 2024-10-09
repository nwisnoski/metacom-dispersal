#!/bin/bash -l

#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64GB
#SBATCH --job-name=metacom_dispersal_dormancy

cd /mnt/scratch/wisnoskilab/GitHub/mc-seedbank-dispersal

module load R/4.4.0

R CMD BATCH --no-restore --no-save simulations/metacom_dispersal-dormancy.R
