#!/bin/bash -l

#SBATCH --account=coexistence
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nathan.wisnoski@uwyo.edu
#SBATCH --job-name=metacom_dispersal_kernel

cd /gscratch/nwisnosk/GitHub/metacom-dispersal

module load swset/2018.05  gcc/7.3.0 r/4.0.5-py27

R CMD BATCH --no-restore --no-save simulations/metacom_diversity_variability.R
