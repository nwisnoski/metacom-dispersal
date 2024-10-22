library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(readr)
library(purrr)
library(data.table)
library(vegan)
library(foreach)
library(doParallel)
library(here)
library(lubridate)
source(here("simulations/metacom_functions.R"))
source(here("analysis/metacommunity_variability_partitioning.R"))
source(here("analysis/diversity_partitioning.R"))


# define simulation parameters
nreps <- 1
x_dim <- 100
y_dim <- 100
patches <- 100
species <- 40
timesteps <- 100
initialization <- 200
burn_in <- 800
full_grid <- FALSE
write_dynamics <- FALSE

# temporal environmental noise parameters
temp_noise_color_vec = 0 #c(-0.8, 0, 0.8) # temporal autocorrelation, can range from -1 (blue noise) to 0 (white noise), to +1 (red noise)
temp_noise_sd_vec = 0 #c(0, 0.1, 0.5)

# temporal trends (set cycles very high to have no trend)
cycles <- 10 # how long are interannual cycles in years/timesteps?
w <- (2*pi) / cycles # angular frequency, oscillations per time interval in radians per second
A <- 5 # amplitude, peak deviations
phi <- 0 # phase, where in the cycle oscillation is at t=0  

# relative spatial heterogeneity (0 = none)
spat_hetero_vec <- c(0, .1, 1000) # spatial heterogeneity, higher numbers, steeper slopes
env_filter_strength <- 0.1

# competition scenarios
conditions <- c("equal", "stable") # can also be "equal" or "priority"
niche_breadth <- 0.25

# dispersal, dormancy, and disturbance rates
disp_rates <- 10^seq(-5, 0, length.out = 20)
kernel_vals <- 10^seq(-4, -0.5, length.out = 10) # c(0, .1, .25, .5)
germ_fracs <- 1
surv_fracs <- 0
sb_responsive <- c(TRUE, FALSE) # is germination responsive to the environment
sb_sensitivity <- 10 # bigger is steeper transition with env condition
sb_maxgerm <- .9 # with responsive sb, does every individual (1) or say, most (.9) germinate? i.e. is there bet hedging?
disturbance_rates <- c(0, 0.01)

# find all combinations of parameters
params <- expand.grid(disp_rates, kernel_vals, disturbance_rates)
# empty, no seedbank
seedbank_params <- data.frame(germ = c(1),
                              surv = c(0))
seedbank_params <- cbind.data.frame(seedbank_params, 
                                    responsive = rep(sb_responsive, each = nrow(seedbank_params)))
seedbank_params <- seedbank_params[1,] # just focus on no seedbank

# source("simulations/test_params.R")

# initialize simulation run
cl <- parallel::makeCluster(16)
registerDoParallel()
start_sim <- Sys.time()

dynamics_total <- data.table()
temp_dyn_per_patch_total <- data.table()
spat_dyn_over_time_total <- data.table()

# set.seed(9468751)
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
          labs(title = paste0("rep = ",rep,", noise color = ", temp_noise_color, ", noise sd = ", temp_noise_sd,", spat het = ", spat_heterogeneity, ", \nA = ", A, ", w = ", w, ", phi = ", phi)) +
          theme_minimal() +
          scale_x_continuous(limits = c(800, 900))
        ggsave(plot = env_plot, filename = paste0("figures/env_plots/env_plot_", rep,"_", temp_noise_color,"_", temp_noise_sd,"_", spat_heterogeneity,"_", as.character(format(Sys.time(), "%X")), ".pdf"), width = 6, height = 4)
        
        
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
          
          for(j in 1:nrow(seedbank_params)){
            germ = seedbank_params[j,1]
            surv = seedbank_params[j,2]
            responsive = seedbank_params[j,3]
              
              print(paste("germ =", germ, "; surv =", surv, "; spat_het =", spat_heterogeneity))
              
              # up until this point, parameters are getting set up for this run 
              
              dynamics_list <- foreach(p = 1:nrow(params), 
                                       .inorder = FALSE,
                                       .errorhandling = 'pass',
                                       .packages = c("tidyr", "dplyr", "ggplot2", "readr", "stringr", "purrr", "data.table", "stats")) %dopar% {
                                         
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
                                                                          env_niche_breadth = niche_breadth, 
                                                                          env_niche_optima = "even",
                                                                          responsive = responsive,
                                                                          sb_sensitivity = sb_sensitivity,
                                                                          sb_maxgerm = sb_maxgerm)
                                           
                                           disp_array <- generate_dispersal_matrices(landscape, species, patches, species_traits, torus = FALSE)
                                           # int_mat <- species_int_mat(species = species, intra = intra,
                                           #                            min_inter = min_inter, max_inter = max_inter,
                                           #                            comp_scaler = comp_scaler, plot = TRUE)
                                           
                                           
                                           N <- init_community(initialization = initialization, species = species, patches = patches)
                                           N <- N + 1 # no initial dispersal limitation
                                           D <- N*0
                                           
                                           for(i in 1:(initialization + burn_in + timesteps)){
                                             if(i <= initialization){
                                               if(i %in% seq(10, 200, by = 10)){
                                                 N <- N + matrix(rpois(n = species*patches, lambda = 0.5), nrow = patches, ncol = species)
                                                 D <- D + matrix(rpois(n = species*patches, lambda = 0.5), nrow = patches, ncol = species)
                                               }
                                               env <- env_df$env[env_df$time == 1]
                                             } else {
                                               env <- env_df$env[env_df$time == (i - initialization)]
                                             }
                                             
                                             # compute r; -1 allows species to have negative fitness in patches with big mismatch
                                             r <- compute_r_xt(species_traits, env = env, species = species) - env_filter_strength
                                             
                                             
                                             # who germinates? Binomial distributed, fixed or responsive
                                             N_germ <- germination(N + D, species_traits, r)
                                             
                                             # record fitness component after environmental effects
                                             # several components will be recorded
                                             fitness_max <- N_germ %*% diag(species_traits$max_r)
                                             fitness_env <- r * N_germ
                                             fitness_env[is.na(fitness_env)] <- 0
                                             
                                             
                                             # of germinating fraction, grow via BH model
                                             N_hat <- growth(N_germ, species_traits, r, int_mat)
                                             
                                             comp_partitions <- get_comp_effects(N_germ, species_traits, r, int_mat)
                                             
                                             # next extract the actual biotic effect, ie "full" option of previous function
                                             fitness_biotic <- N_hat # N after accounting for growth and competition in each patch
                                             fitness_biotic[is.na(fitness_biotic)] <- 0
                                             
                                             # of those that didn't germinate, compute seed bank survival via binomial draw
                                             D_hat <- survival((N + D - N_germ), species_traits) 
                                             
                                             N_hat[N_hat < 0] <- 0
                                             
                                             N_hat_deterministic <- N_hat # snapshot before demographic stochasticity
                                             
                                             N_hat <- matrix(rpois(n = species*patches, lambda = N_hat), ncol = species, nrow = patches) # poisson draw on aboveground
                                             
                                             # record number of stochastic extinctions and fitness after accounting for demographic stochasticity
                                             stoch_extinct_demo <- (N_hat == 0) & (N_hat_deterministic != 0)
                                             fitness_stoch_demo <- N_hat # N after demo stoch
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
                                             
                                             rescued_by_dispersal <- (N_germ != 0) & (N_hat == 0) & (I > 0)
                                             extinct_by_emigration <- ((N_hat - E) <= 0) & (N_hat != 0)
                                             
                                             # measure per capita population growth after dispersal
                                             fitness_dispersal <- (N_hat - E + I) # N after emigration and immigration
                                             fitness_dispersal[is.na(fitness_dispersal)] <- 0
                                             colonizations <- (N_germ == 0) & ((N_hat - E + I) > 0)
                                             
                                             N <- N_hat - E + I
                                             D <- D_hat 
                                             
                                             # next impose disturbance
                                             N_before_disturb <- N
                                             N[rbinom(n = species * patches, size = 1, prob = extirp_prob) > 0] <- 0
                                             
                                             stoch_extinct_env <- (N - N_before_disturb) < 0
                                             fitness_stoch_env <- N # N after any disturbances
                                             fitness_stoch_env[is.na(fitness_stoch_env)] <- 0
                                             
                                             # At this point, N contains a snapshot of the community at time i
                                             dynamics_i <- data.table(N = c(N),
                                                                      species = rep(1:species, each = patches),
                                                                      patch = 1:patches,
                                                                      
                                                                      # fitness and competition strengths
                                                                      W_max = c(fitness_max), # max fitness without any limits
                                                                      W_env = c(fitness_env), # fitness with env filtering
                                                                      W_bio = c(fitness_biotic), # fitness after env, competition
                                                                      #W_bio_full = c(comp_partitions$full), # effects of intra and inter competition on growth, same as W_bio
                                                                      W_bio_inter = c(comp_partitions$inter), # effects of inter competition on growth
                                                                      W_bio_intra = c(comp_partitions$intra), # effects of intra competition on growth
                                                                      #W_bio_none = c(comp_partitions$nocomp), # effects of no competition on growth, same as W_env
                                                                      W_stoch_demo = c(fitness_stoch_demo), # fitness after env, comp, demo stoch
                                                                      W_dispersal = c(fitness_dispersal), # fitness after env, comp, demo stoch, dispersal
                                                                      W_stoch_env = c(fitness_stoch_env), # fitness after env, comp, demo stoch, dispersal, env stoch
                                                                      
                                                                      
                                                                      # process inferences
                                                                      extinct_by_emigration = c(extinct_by_emigration), # did sp go extinct by too high emigration rate
                                                                      rescued_by_dispersal = c(rescued_by_dispersal), # was sp rescued by immigration in this timestep
                                                                      colonization = c(colonizations), # was there a colonization event
                                                                      stoch_extinct_demo = c(stoch_extinct_demo), # did sp go extinct through demo stoch
                                                                      stoch_extinct_env = c(stoch_extinct_env), # did sp go extinct through env stoch
                                                                      
                                                                      # traits
                                                                      emigration = disp,
                                                                      kernel_exp = kernel_exp,
                                                                      germ_rate = germ,
                                                                      surv_rate = surv,
                                                                      responsive = responsive,
                                                                      
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
                                           if(write_dynamics){
                                             current_time <- Sys.time()
                                             tstamp <- str_replace_all(current_time, " ", "_") %>% 
                                               str_replace_all(":", "")
                                             write_csv(x = dynamics_out, col_names = TRUE, 
                                                       file = here(paste0("sim_output/dynamics/disp_kernel_", tstamp ,".csv")))
                                           }
                                           
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
                                           
                                           # summarize local pop dynamics
                                           dynamics_out <- dynamics_out |> 
                                             mutate(sink_env = ((N > 0) & (W_env < 0) & (W_dispersal > 0)),
                                                    sink_biotic = ((N > 0) & (W_env > 0) & (W_bio <= 0) & (W_dispersal > 0)),
                                                    sink_stoch_demo = ((N > 0) & (W_env > 0) & (W_bio > 0) & (W_stoch_demo <= 0) & (W_dispersal > 0)), 
                                             )
                                           
                                           # over time, per patch
                                           pop_dyn_temporal_per_patch <- dynamics_out |> 
                                             group_by(rep, species, patch, emigration, kernel_exp, germ_rate, surv_rate, responsive, extirp_prob, comp, spat_heterogeneity) |> 
                                             left_join(select(species_traits, species, 
                                                              env_niche_optima, env_niche_breadth), by = "species") |> 
                                             summarize(abundance_mean = mean(N),
                                                       abundance_var = var(N),
                                                       occupancy = sum(N>0),
                                                       
                                                       # N after each process
                                                       W_max = mean(W_max), # N if no limits on growth
                                                       W_env_mean = mean(W_env), # N after env filtering
                                                       W_bio_mean = mean(W_bio), # N after env, competition
                                                       W_bio_inter_mean = mean(W_bio_inter), # N after env, and inter comp
                                                       W_bio_intra_mean = mean(W_bio_intra), # N after env, and intra comp
                                                       W_stoch_demo_mean = mean(W_stoch_demo), # N after env, comp, demo stoch
                                                       W_dispersal_mean = mean(W_dispersal), # N after env, comp, demo stoch, dispersal
                                                       W_stoch_env_mean = mean(W_stoch_env), # N after env, comp, demo stoch, dispersal, and disturbance
                                                       
                                                       # delta N due to each process
                                                       delta_env = mean(W_env - W_max),
                                                       delta_bio = mean(W_bio - W_env),
                                                       delta_bio_inter = mean(W_bio_inter - W_env),
                                                       delta_bio_intra = mean(W_bio_intra - W_env),
                                                       delta_stoch_demo = mean(W_stoch_demo - W_bio),
                                                       delta_dispersal = mean(W_dispersal - W_stoch_demo),
                                                       delta_stoch_env = mean(W_stoch_env - W_dispersal),
                                                       
                                                       # env properties
                                                       env_mean = mean(env),
                                                       env_mismatch = mean(abs(env - env_niche_optima)), # how far is this species from its optimum
                                                       extinct_by_emigration = sum(extinct_by_emigration),
                                                       rescued_by_dispersal = sum(rescued_by_dispersal),
                                                       colonization = sum(colonization),
                                                       stoch_extinct_demo = sum(stoch_extinct_demo),
                                                       stoch_extinct_env = sum(stoch_extinct_env),
                                                       
                                                       # sink classification
                                                       sink_env = sum(sink_env),
                                                       sink_biotic = sum(sink_biotic),
                                                       sink_stoch_demo = sum(sink_stoch_demo)) |> 
                                             filter(abundance_mean != 0)
                                           
                                           # across space, per time
                                           pop_dyn_spatial_per_time <- dynamics_out |> 
                                             group_by(rep, species, time, emigration, kernel_exp, germ_rate, surv_rate, responsive, extirp_prob, comp, spat_heterogeneity) |> 
                                             left_join(select(species_traits, species, 
                                                              env_niche_optima, env_niche_breadth), by = "species") |> 
                                             summarize(abundance_mean = mean(N),
                                                       abundance_var = var(N),
                                                       occupancy = sum(N>0),
                                                       
                                                       # N after each process
                                                       W_max = mean(W_max), # N if no limits on growth
                                                       W_env_mean = mean(W_env), # N after env filtering
                                                       W_bio_mean = mean(W_bio), # N after env, competition
                                                       W_bio_inter_mean = mean(W_bio_inter), # N after env, and inter comp
                                                       W_bio_intra_mean = mean(W_bio_intra), # N after env, and intra comp
                                                       W_stoch_demo_mean = mean(W_stoch_demo), # N after env, comp, demo stoch
                                                       W_dispersal_mean = mean(W_dispersal), # N after env, comp, demo stoch, dispersal
                                                       W_stoch_env_mean = mean(W_stoch_env), # N after env, comp, demo stoch, dispersal, and disturbance
                                                       
                                                       # delta N due to each process
                                                       delta_env = mean(W_env - W_max),
                                                       delta_bio = mean(W_bio - W_env),
                                                       delta_bio_inter = mean(W_bio_inter - W_env),
                                                       delta_bio_intra = mean(W_bio_intra - W_env),
                                                       delta_stoch_demo = mean(W_stoch_demo - W_bio),
                                                       delta_dispersal = mean(W_dispersal - W_stoch_demo),
                                                       delta_stoch_env = mean(W_stoch_env - W_dispersal),
                                                       
                                                       # env properties
                                                       env_mean = mean(env),
                                                       env_mismatch = mean(abs(env - env_niche_optima)), # how far is this species from its optimum
                                                       extinct_by_emigration = sum(extinct_by_emigration),
                                                       rescued_by_dispersal = sum(rescued_by_dispersal),
                                                       colonization = sum(colonization),
                                                       stoch_extinct_demo = sum(stoch_extinct_demo),
                                                       stoch_extinct_env = sum(stoch_extinct_env),
                                                       
                                                       # sink classification
                                                       sink_env = sum(sink_env),
                                                       sink_biotic = sum(sink_biotic),
                                                       sink_stoch_demo = sum(sink_stoch_demo)) |> 
                                             filter(abundance_mean != 0)
                                           
                                           
                                           
                                           
                                           
                                           # make an output table
                                           output_summary <- data.table(
                                             rep = rep,
                                             condition = x,
                                             disp_rate = params[p,1],
                                             kernel_exp = params[p,2],
                                             germ_rate = germ,
                                             surv_rate = surv,
                                             responsive = responsive,
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
                                           
                                           # return(output_summary)
                                           return(list(summary = output_summary,
                                                       spatial_per_time = pop_dyn_spatial_per_time,
                                                       time_per_patch = pop_dyn_temporal_per_patch)) # this return means this is final information taken into dynamics_list
                                         }, silent = FALSE, outFile = "errors.log")
                                       }
              
              summary_list <- dynamics_list |> map_depth(1, "summary")
              fails <- (as.logical(lapply(summary_list, is_try_error)) + as.logical(lapply(summary_list, is_simple_error)))
              
              dynamics_total_i <- rbindlist(summary_list[!fails]) # only binds the runs without errors
              dynamics_total <- bind_rows(dynamics_total, dynamics_total_i)
              
              temp_dyn_per_patch_list <- dynamics_list |> map_depth(1, "time_per_patch")
              temp_dyn_per_patch_i <- rbindlist(temp_dyn_per_patch_list[!fails])
              temp_dyn_per_patch_total <- bind_rows(temp_dyn_per_patch_total, temp_dyn_per_patch_i)
              
              spat_dyn_over_time_list <- dynamics_list |> map_depth(1, "spatial_per_time") 
              spat_dyn_over_time_i <- rbindlist(spat_dyn_over_time_list[!fails])
              spat_dyn_over_time_total <- bind_rows(spat_dyn_over_time_total, spat_dyn_over_time_i)
              
              print(paste("---- completed for condition:", x))
              
              
              
            }
        }
      }
    }
  }
}
end_sims <- Sys.time()
tstamp <- str_replace_all(end_sims, " ", "_") %>% 
  str_replace_all(":", "") |> 
  str_replace("[.]", "_")
dir.create(paste0("sim_output/", lubridate::date(tstamp)))
write_csv(x = dynamics_total, col_names = TRUE, 
          file = here(paste0("sim_output/", lubridate::date(tstamp), "/variability_partitioning_disp_kernel_", tstamp ,"_summary.csv")))
write_csv(x = spat_dyn_over_time_total, col_names = TRUE, 
          file = here(paste0("sim_output/", lubridate::date(tstamp), "/variability_partitioning_disp_kernel_", tstamp ,"_spat_per_time.csv")))
write_csv(x = temp_dyn_per_patch_total, col_names = TRUE, 
          file = here(paste0("sim_output/", lubridate::date(tstamp), "/variability_partitioning_disp_kernel_", tstamp ,"_temp_per_patch.csv")))
