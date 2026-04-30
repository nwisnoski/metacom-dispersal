# metacom-dispersal: Exploring the metacommunity ecology of dispersal kernels.

## Corresponding paper
Code to reproduce the analyses and figures in the manuscript: 
Wisnoski NI, Szojka MC, Germain RM, Fukami T, Shoemaker LG. "Dispersal kernels regulate environmental, biotic, and stochastic effects on the maintenance of metacommunity diversity"

## Repository organization
The repository contains the following folders:
- `simulations` = This folder contains the scripts necessary to run the simulations and generate simulation data.
- `sim_output` = This folder is the destination where simulation scripts will deposit output files, organized by subdirectories labeled by  date. The files will be timestamped as well to avoid overwriting. 
- `analysis` = This folder contains the scripts necessary for analyses performed on simulated data. 
- `figures` = This folder is where generated figures will be located, labeled by figure number.
- `appendix` = This folder is self-contained and contains the scripts necessary to generate the analyses and figures in the Appendix S1.

## File descriptions
There are two main scripts in this project: 
- `./simulations/metacom-dispersal-kernel.R` is the primary simulation model. Parameters are set in the main script, and the script by default performs a full parameter sweep for a single landscape replicate. Replication can be accomplished by running multiple instances of this script (which will deposit multiple files to the `sim_output` folder). It is also possible to have the single script perform replicates using the replicate flag. All other parameters are described in comments in the script. 
- `./analyses/analyze_diversity.R` is the primary analysis script, which orchestrates the figure generation for the main manuscript. To run this file, you must tell it the location where the focal simulation data is located (i.e., which subdirectory of `sim_output` to read files from). Then, all figures are generated in order. 

Additionally, the following scripts generate figures in the analysis:
- `./appendix/metacom_dispersal-kernel_20sp.R` (and `_75sp.R`) perform the simulations with differing sizes in the regional species pool, which revealed consistent patterns across a gradient in regional diversity. `/appendix/analyze_diversity_specnum-variation.R` is the corresponding analysis of these simulations with different regional species pools. 
- `./appendix/dispersal_distance_analysis.R` is a script that compares the distribution of all pairwise landscape distances with the distribution of dispersal distances under dispersal kernels with different kernel exponents. 
- `./appendix/visualize_connectivity.R` generates network figures showing the connectivity impacts of changes in kernel shape, to further demonstrate that kernels regulate the spatial span/configuration of patches in the metacommunity.

## Generate figures
1. Run `./simulations/metacom_dispersal-kernel.R` locally or deploy `./simulations/run_simulations.sh` to a slurm-based cluster modified for your particular needs, with your preferred level of replication.
2. Output will be placed in `./sim_output` organized by date.
3. Run `./analysis/analyze_diversity.R` to generate figures.
4. Figures will populate in `./figures`.

