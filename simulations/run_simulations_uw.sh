#!/bin/bash -l

#SBATCH --account=coexistence
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=64GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nathan.wisnoski@msstate.edu
#SBATCH --job-name=metacom_dispersal_kernel

cd /gscratch/nwisnosk/GitHub/metacom-dispersal

module load gcc/12.2.0 openmpi/4.1.4 r/4.2.2

R CMD BATCH --no-restore --no-save simulations/metacom_diversity_variability.R
