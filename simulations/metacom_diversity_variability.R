library(tidyverse)
library(data.table)
library(progress)
library(vegan)
library(foreach)
library(doParallel)
library(here)
source(here("simulations/metacom_functions.R"))
source(here("analysis/metacommunity_variability_partitioning.R"))
source(here("analysis/diversity_partitioning.R"))


# define parameters
nreps <- 1
x_dim <- 100
y_dim <- 100
patches <- 100
species <- 40
full_grid <- FALSE

# environmental var params
temp_noise_color_vec = 0 #c(-0.8, 0, 0.8) # temporal autocorrelation, can range from -1 (blue noise) to 0 (white noise), to +1 (red noise)
temp_noise_sd_vec = 0 #c(0, 0.1, 0.5)

# temporal trends (set cycles very high to have no trend)
cycles <- 100 # how long are interannual cycles in years/timesteps?
w <- (2*pi) / cycles # angular frequency, oscillations per time interval in radians per second
A <- 5 # amplitude, peak deviations
phi <- 0 # phase, where in the cycle oscillation is at t=0  

# spatial heterogeneity, higher numbers, steeper slopes, 
# more spatial relative to temporal. 
# lower numbers have temporal but not spatial var
# higher numbers have spatial but not temporal var
spat_hetero_vec = c(0, .1, 1, 1000) 

conditions <- c("stable", "equal", "priority") # can also be "equal" or "priority"

timesteps <- 100
initialization <- 200
burn_in <- 800

# run sim
disp_rates <- 10^seq(-5, 0, length.out = 20)
kernel_vals <- seq(0, 1, length.out = 10)
disturbance_rates <- seq(0.00, 0.04, length.out = 5)

# remove seed bank
germ_fracs <- 1
surv_fracs <- 0

params <- expand.grid(disp_rates, kernel_vals, disturbance_rates)

cl <- parallel::makeCluster(10)
registerDoParallel()
start_sim <- Sys.time()

dynamics_total <- data.table()
#set.seed(9468751)
# for each replicate, rerun parameter sweep
for(rep in 1:nreps){
  print(paste("Running rep", rep, "of", nreps))
  landscape <- if(full_grid == TRUE){
    expand.grid(x = 1:x_dim, y = 1:y_dim)
  } else init_landscape(patches = patches, x_dim = x_dim, y_dim = y_dim)
  
  for(temp_noise_color in temp_noise_color_vec){
    for(temp_noise_sd in temp_noise_sd_vec){
      for(spat_heterogeneity in spat_hetero_vec){
        
        
        # make new landscape, environmental data, and draw new competition coefficients
        
        
        env_df <- env_generate(landscape = landscape, 
                               spat_heterogeneity = spat_heterogeneity, # smaller this is compared to A, the more overlap there is among patches
                               temp_noise_color = temp_noise_color, # temporal noise autocorrelation, can range from -1 (blue noise) to 0 (white noise), to +1 (red noise)
                               temp_noise_sd = temp_noise_sd, # temp noise standard deviation, more or less variable
                               timesteps = timesteps+burn_in, A=A, w=w, phi=phi)
        env_plot <- ggplot(env_df, aes(x = time, y = env, color = patch)) + 
          geom_line(alpha = 0.5, show.legend = FALSE) + 
          labs(title = paste0("rep = ",rep,", noise color = ", temp_noise_color, ", noise sd = ", temp_noise_sd,", spat het = ", spat_heterogeneity)) +
          theme_minimal()
        #ggsave(plot = env_plot, filename = paste0("figures/env_plots/env_plot_", rep, temp_noise_color, temp_noise_sd, spat_heterogeneity,".pdf"), width = 6, height = 4)
        
        
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
          
          for(germ in germ_fracs){
            for(surv in surv_fracs){
              
              print(paste("germ =", germ, "; surv =", surv, "; spat_het =", spat_heterogeneity))
              
              # up until this point, parameters are getting set up for this run 
              
              dynamics_list <- foreach(p = 1:nrow(params), 
                                       .inorder = FALSE,
                                       .errorhandling = 'pass',
                                       .packages = c("tidyverse", "data.table", "stats")) %dopar% {
                                         
                                         try(expr = {
                                           
                                           
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
                                           N <- N + 1 # no initial dispersal limitation
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
                                             
                                             # record fitness after environmental effects
                                             fitness_env <- r
                                             fitness_env[is.na(fitness_env)] <- 0
                                             
                                             # who germinates? Binomial distributed
                                             N_germ <- germination(N + D, species_traits)
                                             
                                             # of germinating fraction, grow via BH model
                                             N_hat <- growth(N_germ, species_traits, r, int_mat)
                                             
                                             comp_partitions <- get_comp_effects(N_germ, species_traits, r, int_mat)
                                             
                                             # next extract the actual biotic effect, ie "full" option of previous function
                                             fitness_biotic <- N_hat / N_germ # delta N per capita after accounting for germination numbers in each patch
                                             fitness_biotic[is.na(fitness_biotic)] <- 0
                                             
                                             # of those that didn't germinate, compute seed bank survival via binomial draw
                                             D_hat <- survival((N + D - N_germ), species_traits) 
                                             
                                             N_hat[N_hat < 0] <- 0
                                             
                                             N_hat_deterministic <- N_hat # snapshot before demographic stochasticity
                                             
                                             N_hat <- matrix(rpois(n = species*patches, lambda = N_hat), ncol = species, nrow = patches) # poisson draw on aboveground
                                             
                                             # record number of stochastic extinctions and fitness after accounting for demographic stochasticity
                                             stoch_extinct_demo <- (N_hat == 0) & (N_hat_deterministic != 0)
                                             fitness_stoch_demo <- ((N_hat - N_germ)/N_germ)
                                             fitness_stoch_demo[is.na(fitness_stoch_demo)] <- 0
                                             
                                             # determine emigrants from aboveground community
                                             E <- matrix(nrow = patches, ncol = species)
                                             disp_rates <- species_traits$dispersal_rate
                                             for(s in 1:species){
                                               E[,s] <- rbinom(n = patches, size = (N_hat[,s]), prob = disp_rates[s])
                                             }
                                             
                                             dispSP <- colSums(E)
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
                                             
                                             rescued_by_dispersal <- (stoch_extinct_demo == TRUE) & (I > 0)
                                             extinct_by_emigration <- ((N_hat - E) <= 0) & (N_hat != 0)
                                             
                                             # measure per capita population growth after dispersal
                                             fitness_dispersal <- ((N_hat - E + I) - N_germ)/N_germ
                                             fitness_dispersal[is.na(fitness_dispersal)] <- 0
                                             
                                             N <- N_hat - E + I
                                             D <- D_hat 
                                             
                                             # next impose disturbance
                                             N_before_disturb <- N
                                             N[rbinom(n = species * patches, size = 1, prob = extirp_prob) > 0] <- 0
                                             
                                             stoch_extinct_env <- (N - N_before_disturb) > 0
                                             fitness_stoch_env <- (N - N_germ) / N_germ
                                             fitness_stoch_env[is.na(fitness_stoch_env)] <- 0
                                             
                                             # At this point, N contains a snapshot of the community at time i
                                             dynamics_i <- data.table(N = c(N),
                                                                      species = rep(1:species, each = patches),
                                                                      patch = 1:patches,
                                                                      
                                                                      # fitness and competition strengths
                                                                      W_env = c(fitness_env), # fitness with env filtering
                                                                      W_bio = c(fitness_biotic), # fitness after env, competition
                                                                      W_stoch_demo = c(fitness_stoch_demo), # fitness after env, comp, demo stoch
                                                                      W_dispersal = c(fitness_dispersal), # fitness after env, comp, demo stoch, dispersal
                                                                      W_stoch_env = c(fitness_stoch_env), # fitness after env, comp, demo stoch, dispersal, env stoch
                                                                      comp_effects_full = c(comp_partitions$full), # effects of intra and inter competition on growth
                                                                      comp_effects_inter = c(comp_partitions$inter), # effects of inter competition on growth
                                                                      comp_effects_intra = c(comp_partitions$intra), # effects of intra competition on growth
                                                                      comp_effects_none = c(comp_partitions$nocomp), # effects of no competition on growth
                                                                      
                                                                      # process inferences
                                                                      extinct_by_emigration = c(extinct_by_emigration), # did sp go extinct by too high emigration rate
                                                                      rescued_by_dispersal = c(rescued_by_dispersal), # was sp rescued by immigration in this timestep
                                                                      stoch_extinct_demo = c(stoch_extinct_demo), # did sp go extinct through demo stoch
                                                                      stoch_extinct_env = c(stoch_extinct_env), # did sp go extinct through env stoch
                                                                      
                                                                      # traits
                                                                      emigration = disp,
                                                                      kernel_exp = kernel_exp,
                                                                      comp = x,
                                                                      rep = rep,
                                                                      
                                                                      # environmental variables
                                                                      time = i-initialization-burn_in, 
                                                                      env = env,
                                                                      extirp_prob = extirp_prob,
                                                                      temp_noise_color = temp_noise_color,
                                                                      temp_noise_sd = temp_noise_sd,
                                                                      spat_heterogeneity = spat_heterogeneity,
                                                                      spat_mean = mean(env),
                                                                      spat_sd = sd(env)
                                             ) 
                                             
                                             # add time step i to the full dynamics "dynamics_out"
                                             dynamics_out <- rbind(dynamics_out, 
                                                                   filter(dynamics_i, time > 0))
                                           } # end loop calculating N responses for a timeseries
                                           
                                           # write dynamics to file
                                           write_csv(x = dynamics_out, col_names = TRUE, 
                                                     file = here(paste0("sim_output/dynamics/disp_kernel_", tstamp ,".csv")))
                                           
                                           # here is where do temporal beta and variability partitioning
                                           
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
                                           
                                           env_temp_cv_mean <- dynamics_out |> 
                                             group_by(patch) |> 
                                             summarize(env_temp_cv = sd(env)/mean(env)) |> 
                                             summarize(env_temp_cv_mean = mean(env_temp_cv))
                                           
                                           env_spat_cv_mean <- dynamics_out |> 
                                             group_by(time) |> 
                                             summarize(env_spat_cv = sd(env)/mean(env)) |> 
                                             summarize(env_spat_cv_mean = mean(env_spat_cv))
                                           
                                           
                                           # make an output table
                                           output_summary <- data.table(
                                             rep = rep,
                                             condition = x,
                                             disp_rate = params[p,1],
                                             kernel_exp = params[p,2],
                                             disturb_rate = params[p,3],
                                             temp_noise_color = temp_noise_color,
                                             temp_noise_sd = temp_noise_sd,
                                             spat_heterogeneity = spat_heterogeneity,
                                             env_temp_cv_mean = env_temp_cv_mean,
                                             env_spat_cv_mean = env_spat_cv_mean,
                                             
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
                                         }, silent = FALSE, outFile = "errors.log")
                                       }
              print(paste("---- completed for condition:", x))
              
              fails <- (as.logical(lapply(dynamics_list, is_try_error)) + as.logical(lapply(dynamics_list, is_simple_error)))
              
              dynamics_total_i <- rbindlist(dynamics_list[!fails]) # only binds the runs without errors
              dynamics_total <- bind_rows(dynamics_total, dynamics_total_i)
              
            }
          }
        }
      }
    }
  }
}
end_sims <- Sys.time()
tstamp <- str_replace_all(end_sims, " ", "_") %>% 
  str_replace_all(":", "")
gc()
write_csv(x = dynamics_total, col_names = TRUE, 
          file = here(paste0("sim_output/variability_partitioning_disp_kernel_", tstamp ,".csv")))
