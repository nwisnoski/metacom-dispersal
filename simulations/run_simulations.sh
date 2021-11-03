#!/bin/bash -l

#SBATCH --account=coexistence
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --mem=124GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nathan.wisnoski@uwyo.edu
#SBATCH --job-name=metacom

cd /gscratch/nwisnosk/GitHub/metacom-dispersal

module load swset/2018.05  gcc/7.3.0 r/4.0.5-py27

R CMD BATCH --no-restore --no-save mc_sim.R mc_sim.log
