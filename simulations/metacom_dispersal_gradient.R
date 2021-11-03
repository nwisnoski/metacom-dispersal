library(tidyverse)
library(data.table)
library(progress)
library(RandomFields)
library(vegan)
library(foreach)
library(doParallel)
library(here)
source(here("simulations/metacom_functions.R"))


# define parameters
nreps <- 1
x_dim <- 100
y_dim <- 100
patches <- 100
species <- 40
#extirp_prob <- 0.000

conditions <- c("equal", "stable", "priority")

timesteps <- 2000
initialization <- 200
burn_in <- 800

# run sim
set.seed(82072)

disp_rates <- 10^seq(-5, 0, length.out = 10)
kernel_vals <- seq(0, 1, length.out = 10)
disturbance_rates <- seq(0, 0.1, length.out = 5)

# remove seed bank
germ <- 1
surv <- 0
# germ_fracs <- seq(.1,1, length.out = 10)
# surv_fracs <- c(.1, .5, 1)

params <- expand.grid(disp_rates, kernel_vals, disturbance_rates)

cl <- parallel::makeCluster(8)
registerDoParallel()
start_sim <- Sys.time()

for(x in conditions){
  
  if(x == "equal"){
    intra = 1 
    min_inter = 1 
    max_inter = 1 
    comp_scaler = 0.05
  }
  
  if(x == "stable"){
    intra = 1 
    min_inter = 0 
    max_inter = 1 
    comp_scaler = 0.05
  }
  
  if(x == "priority"){
    intra = 1
    min_inter = 1
    max_inter = 2
    comp_scaler = 0.05
  }
  
  # for each replicate, rerun parameter sweep
  for(rep in 1:nreps){
    
    # make new landscape, environmental data, and draw new competition coefficients
    landscape <- init_landscape(patches = patches, x_dim = x_dim, y_dim = y_dim)
    env_df <- env_generate(landscape = landscape, env1Scale = 500, 
                           timesteps = timesteps+burn_in, plot = FALSE)
    int_mat <- species_int_mat(species = species, intra = intra,
                               min_inter = min_inter, max_inter = max_inter,
                               comp_scaler = comp_scaler, plot = FALSE)
    
    
    dynamics_list <- foreach(p = 1:nrow(params), .inorder = FALSE,
                             .packages = c("tidyverse", "data.table", "stats")) %dopar% {
                               
                               # extract params
                               disp <- params[p,1]
                               kernel_exp <- params[p,2]
                               extirp_prob <- params[p,3]
                               
                               
                               dynamics_out <- data.table()
                               
                               # init_community
                               species_traits <- init_species(species, 
                                                              dispersal_rate = disp,
                                                              kernel_exp = kernel_exp,
                                                              germ = germ,
                                                              survival = surv,
                                                              env_niche_breadth = 0.5, 
                                                              env_niche_optima = "even")
                               
                               disp_array <- generate_dispersal_matrices(landscape, species, patches, species_traits, torus = FALSE)
                               # int_mat <- species_int_mat(species = species, intra = intra,
                               #                            min_inter = min_inter, max_inter = max_inter,
                               #                            comp_scaler = comp_scaler, plot = TRUE)
                               
                               
                               N <- init_community(initialization = initialization, species = species, patches = patches)
                               N <- N + 1
                               D <- N*0
                               
                               for(i in 1:(initialization + burn_in + timesteps)){
                                 if(i <= initialization){
                                   if(i %in% seq(10, 100, by = 10)){
                                     N <- N + matrix(rpois(n = species*patches, lambda = 0.5), nrow = patches, ncol = species)
                                     D <- D + matrix(rpois(n = species*patches, lambda = 0.5), nrow = patches, ncol = species)
                                   }
                                   env <- env_df$env1[env_df$time == 1]
                                 } else {
                                   env <- env_df$env1[env_df$time == (i - initialization)]
                                 }
                                 
                                 # compute r
                                 r <- compute_r_xt(species_traits, env = env, species = species)
                                 
                                 # compute growth
                                 #N_hat <- N*r / (1 + N%*%int_mat)
                                 
                                 # who germinates? Binomial distributed
                                 N_germ <- germination(N + D, species_traits)
                                 
                                 # of germinating fraction, grow via BH model
                                 N_hat <- growth(N_germ, species_traits, r, int_mat)
                                 
                                 # of those that didn't germinate, compute seed bank survival via binomial draw
                                 D_hat <- survival((N + D - N_germ), species_traits) 
                                 
                                 N_hat[N_hat < 0] <- 0
                                 N_hat <- matrix(rpois(n = species*patches, lambda = N_hat), ncol = species, nrow = patches) # poisson draw on aboveground
                                 
                                 # determine emigrants from aboveground community
                                 E <- matrix(nrow = patches, ncol = species)
                                 disp_rates <- species_traits$dispersal_rate
                                 for(s in 1:species){
                                   E[,s] <- rbinom(n = patches, size = (N_hat[,s]), prob = disp_rates[s])
                                 }
                                 
                                 dispSP <- colSums(E)
                                 
                                 # determine immigrants to each patch
                                 # I_hat_raw <- disp_array[,,1]%*%E
                                 # I_hat <- t(t(I_hat_raw)/colSums(I_hat_raw))
                                 # I_hat[is.nan(I_hat)] <- 1
                                 # I <- sapply(1:species, function(x) {
                                 #   if(dispSP[x]>0){
                                 #     table(factor(sample(x = patches, size = dispSP[x], replace = TRUE, prob = I_hat[,x]), levels = 1:patches))
                                 #   } else {rep(0, patches)}
                                 # })
                                 
                                 I_hat_raw <- matrix(nrow = patches, ncol = species)
                                 
                                 for(s in 1:species){
                                   I_hat_raw[,s] <- disp_array[,,s] %*% E[,s]
                                 }
                                 
                                 # standardize so colsums = 1
                                 I_hat <- t(t(I_hat_raw)/colSums(I_hat_raw))
                                 I_hat[is.nan(I_hat)] <- 1
                                 
                                 I <- sapply(1:species, function(x) {
                                   if(dispSP[x]>0){
                                     table(factor(
                                       sample(x = patches, size = dispSP[x], replace = TRUE, prob = I_hat[,x]),
                                       levels = 1:patches))
                                   } else {rep(0, patches)}
                                 })
                                 
                                 N <- N_hat - E + I
                                 D <- D_hat 
                                 
                                 N[rbinom(n = species * patches, size = 1, prob = extirp_prob) > 0] <- 0
                                 
                                 dynamics_i <- data.table(N = c(N),
                                                          patch = 1:patches,
                                                          species = rep(1:species, each = patches),
                                                          env = env,
                                                          time = i-initialization-burn_in, 
                                                          dispersal = disp,
                                                          kernel_exp = kernel_exp,
                                                          extirp_prob = extirp_prob,
                                                          rep = rep,
                                                          comp = x) %>% 
                                   filter(time %in% seq(2000, timesteps, by = 20))
                                 
                                 dynamics_out <- rbind(dynamics_out, 
                                                       dynamics_i)
                               }
                               
                               # every 20 tsteps?
                               # dynamics_subset <- dynamics_out %>% 
                               #   filter(time %in% seq(0, timesteps, by = 20))
                               
                               #dynamics_total <- rbind(dynamics_total, dynamics_subset)
                               
                               #saveRDS(dynamics_out, file = paste0("sim_output/sim_disp",disp,"_germ_",germ,"_surv_",surv,"_maxinter_",max_inter,"_mininter_",min_inter,".rds"))
                               
                               
                               return(dynamics_out)
                             }
    
    dynamics_total <- rbindlist(dynamics_list)
    end_sims <- Sys.time()
    tstamp <- str_replace_all(end_sims, " ", "_") %>% 
      str_replace_all(":", "")
    
    write_csv(x = dynamics_total, col_names = TRUE, 
              file = here(paste0("sim_output/",x,"_", tstamp ,".csv")))
    rm(dynamics_total)
    gc()
  }
}

