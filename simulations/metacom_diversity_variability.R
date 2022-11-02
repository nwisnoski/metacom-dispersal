library(tidyverse)
library(data.table)
library(progress)
library(fields)
library(multitaper)
library(vegan)
library(foreach)
library(doParallel)
library(here)
source(here("simulations/metacom_functions.R"))
source(here("analysis/metacommunity_variability_partitioning.R"))
source(here("analysis/diversity_partitioning.R"))


# define parameters
nreps <- 10
x_dim <- 100
y_dim <- 100
patches <- 100
species <- 40
full_grid <- FALSE
#extirp_prob <- 0.000

# environmental var params
temp_auto_vec = c(-0.8, 0, 0.8) # temporal autocorrelation, can range from -1 (blue noise) to 0 (white noise), to +1 (red noise)
spat_auto_vec = c(0.0001, 5, 10) # spatial autocorrelation, less clustered toward 0, more clustered toward 1

# traveling sine wave params
A <- 0.5
k <- 1/20 # wave number, wavelengths per unit distance
w <- 1/100 # angular frequency, relates to speed of propagation, v = w/k
phi <- 0 # phase, where in the cycle oscillation is at t=0  


conditions <- c("stable")

timesteps <- 100
initialization <- 200
burn_in <- 800

# run sim
set.seed(82072)

disp_rates <- 10^seq(-5, 0, length.out = 20)
kernel_vals <- seq(0, 1, length.out = 10)
disturbance_rates <- seq(0, 0.5, length.out = 10)

# remove seed bank
germ <- 1
surv <- 0
# germ_fracs <- seq(.1,1, length.out = 10)
# surv_fracs <- c(.1, .5, 1)

params <- expand.grid(disp_rates, kernel_vals, disturbance_rates)

cl <- parallel::makeCluster(16)
registerDoParallel()
start_sim <- Sys.time()

dynamics_total <- data.table()

# for each replicate, rerun parameter sweep
for(rep in 1:nreps){
  
  for(temp_auto in temp_auto_vec){
    
    for(spat_auto in spat_auto_vec){
      
      
      
      # make new landscape, environmental data, and draw new competition coefficients
      landscape <- if(full_grid == TRUE){
        expand.grid(x = 1:x_dim, y = 1:y_dim)
      } else init_landscape(patches = patches, x_dim = x_dim, y_dim = y_dim)
      
      env_df <- env_generate(landscape = landscape, x_dim = x_dim, y_dim = y_dim, 
                             spat_auto = spat_auto, # spatial autocorrelated, less clustered toward 0, more clustered toward 1
                             temp_auto = temp_auto, # temporal autocorrelation, can range from -1 (blue noise) to 0 (white noise), to +1 (red noise)
                             timesteps = timesteps+burn_in, A=A, k=k, w=w, phi=phi)
      
      print(paste("Running rep", rep, "of", nreps))
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
        int_mat <- species_int_mat(species = species, intra = intra,
                                   min_inter = min_inter, max_inter = max_inter,
                                   comp_scaler = comp_scaler, plot = FALSE)
        
        print(paste("Simulating condition:", x))
        # up until this point, parameters are getting set up for this run 
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
                                       env <- env_df$env[env_df$time == 1]
                                     } else {
                                       env <- env_df$env[env_df$time == (i - initialization)]
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
                                     
                                     # At this point, N contains a snapshot of the community at time i
                                     dynamics_i <- data.table(N = c(N),
                                                              patch = 1:patches,
                                                              species = rep(1:species, each = patches),
                                                              env = env,
                                                              time = i-initialization-burn_in, 
                                                              dispersal = disp,
                                                              kernel_exp = kernel_exp,
                                                              extirp_prob = extirp_prob,
                                                              rep = rep,
                                                              comp = x,
                                                              temp_auto = temp_auto,
                                                              spat_auto = spat_auto) 
                                     
                                     # add time step i to the full dynamics "dynamics_out"
                                     dynamics_out <- rbind(dynamics_out, 
                                                           filter(dynamics_i, time > 0))
                                   } # end loop calculating N responses for a timeseries
                                   
                                   # here is where do temporal beta and variability paritioning
                                   
                                   
                                   new <- dynamics_out %>% 
                                     dplyr::select(N, patch, species, time) %>%  
                                     dplyr::filter(time > 0) %>%
                                     pivot_wider(., names_from = "time", values_from = "N") # this form means species are first dim and time is second dim
                                   
                                   metacomm_tsdata <- array(NA, c(species, timesteps, patches)) # array N*T*M where N = abundance of each species, T = timeseries, M= patch
                                   for(m in 1:patches){
                                     temp <- new %>%
                                       dplyr::filter(patch == m) %>%
                                       dplyr::select(-species, -patch)
                                     metacomm_tsdata[,,m] <- as.matrix(temp)
                                   }
                                   
                                   # partition variability
                                   par <- var.partition(metacomm_tsdata)
                                   
                                   
                                   # partition diversity 
                                   div <- diversity.partition(metacomm_tsdata)
                                   
                                   # beta div calculations
                                   beta <- beta.div.calc(metacomm_tsdata)
                                   
                                   # local diversity-stability calculations
                                   local_dsr <- div.stab.comp(metacomm_tsdata)
                                   
                                   # make an output table
                                   output_summary <- data.table(
                                     rep = rep,
                                     condition = x,
                                     disp_rate = params[p,1],
                                     kernel_exp = params[p,2],
                                     disturb_rate = params[p,3],
                                     temp_auto = temp_auto,
                                     spat_auto = spat_auto,
                                     CV_S_L = par[1],
                                     CV_C_L = par[2],
                                     CV_S_R = par[3],
                                     CV_C_R = par[4],
                                     phi_S_L2R = par[5],
                                     phi_C_L2R = par[6],
                                     phi_S2C_L = par[7],
                                     phi_S2C_R = par[8],
                                     alpha_div = div[1],
                                     beta_div = div[2], 
                                     gamma_div = div[3],
                                     beta_spatial = beta$mean_beta_spatial,
                                     beta_temporal = beta$mean_beta_temporal,
                                     local_dsr_richness = local_dsr$local_mean_richness,
                                     local_dsr_cv = local_dsr$local_mean_cv
                                     # add spatial div
                                     # add beta
                                   )
                                   
                                   # output should contain:
                                   # all the parameters and settings and identifiable for this run
                                   # rep specific = rep
                                   # condition specific = x 
                                   # responses that are condition specific and calculated using time_series_i: var metrics, synchrony metrics, spatial diversity
                                   # patch specific things: dispersal_rate = disp, disturbance_rate = disturbance_rates, kernel = kernal_vals (don't want to use these because too specific)
                                   # spatial diversity at final time point (average of last few points?), would use time_series_i
                                   # the four variability metrics: population, metapop, community, metacommunity? 
                                   # Do we need to save the synchrony metrics too? - they are outputs so yes
                                   # temporal beta diversity - BD de caceres and legendre paper
                                   
                                   # check that "output_summary" is 1 row, and a lot of columns
                                   
                                   return(output_summary) # this return means this is final information taken into dynamics_list
                                 }
        
        dynamics_total_i <- rbindlist(dynamics_list)
        dynamics_total <- bind_rows(dynamics_total, dynamics_total_i)
        end_sims <- Sys.time()
        tstamp <- str_replace_all(end_sims, " ", "_") %>% 
          str_replace_all(":", "")
        gc()
      }
    }
  }
}

write_csv(x = dynamics_total, col_names = TRUE, 
          file = here(paste0("sim_output/variability_partitioning_", tstamp ,".csv")))
