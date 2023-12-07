# dist.torus function from the packages som.nn
dist.torus <- function (coors) 
{
  x <- coors[, 1]
  dx <- stats::dist(x, diag = TRUE, upper = TRUE)
  y <- coors[, 2]
  dy <- stats::dist(y, diag = TRUE, upper = TRUE)
  max.x <- max(dx) + 1
  max.y <- max(dy) + 1
  mdx <- max.x - dx
  mdy <- max.y - dy
  dx <- pmin(dx, mdx)
  dy <- pmin(dy, mdy)
  d <- sqrt(dx^2 + dy^2)
  return(d)
}


# initialize_landscape
init_landscape <- function(patches, x_dim = 100, y_dim = 100){
  repeat{
    landscape <- data.frame(x = sample(1:x_dim, size = patches, replace = T),
                            y = sample(1:y_dim, size = patches, replace = T))
    if(dim(unique(landscape))[1] == patches) {break}
  }
  
  # plot(landscape)
  
  return(landscape)
}


init_species <- function(species = 10, 
                         env_niche_optima = "random",
                         env_niche_breadth = 0.5,
                         kernel_exp = 0.1,
                         max_r = 5,
                         dispersal_rate = 0.1,
                         survival = 0.8,
                         germ = 0.5){
  
  # generate niche optima

  if(env_niche_optima == "random"){
    optima <- runif(n = species, min = 0, max = 1)
  } else if(env_niche_optima == "even"){
    optima <- seq(from = 0, to = 1, length = species)
  } else if(env_niche_optima == "identical"){
    optima <- rep(runif(1), species)
  } else stop("Enter a env_niche_optima of 'random', 'even' or 'identical'.")
  
  # generate niche breadths
  if(length(env_niche_breadth) != species & length(env_niche_breadth) != 1){
    stop("Enter a number or a vector of length 'species' for env_niche_breadth.")
  }
  
  # generate dispersal rates
  if(length(kernel_exp) != species & length(kernel_exp) != 1){
    stop("Enter a number or a vector of length 'species' for kernel_exp.")
  }
  if(length(dispersal_rate) != species & length(dispersal_rate) != 1){
    stop("Enter a number or a vector of length 'species' for dispersal_rate.")
  }
  
  # generate max growth rates
  if(length(max_r) != species & length(max_r) != 1){
    stop("Enter a number or a vector of length 'species' for max_r.")
  }
  
  species_traits <- data.frame(
    species = 1:species,
    max_r = max_r,
    env_niche_optima = optima,
    env_niche_breadth = env_niche_breadth,
    kernel_exp = kernel_exp,
    dispersal_rate = dispersal_rate,
    survival = survival,
    germ = germ
  )
  
  # matplot(sapply(X = 1:species, FUN = function(x) {
  #     exp(-((species_traits$env_niche_optima[x]-seq(0, 1, length = 30))/(2*species_traits$env_niche_breadth[x]))^2)
  #   })*rep(max_r,each = 30), 
  #   type = "l", lty = 1, ylab = "r", xlab = "environment", ylim = c(0,max(max_r)))
  
  return(species_traits)

}

init_community <- function(initialization = 200, species = 10, patches = 100){
  
  N <- matrix(rpois(n = species * patches, lambda = 0.5), nrow = patches, ncol = species)
  
  return(N)
}


generate_dispersal_matrices <- function(landscape, species, 
                                        patches = patches, 
                                        species_traits, torus = TRUE){
  
  if(torus == TRUE){
    dist_mat <- as.matrix(dist.torus(coors = landscape))
  } else {
    dist_mat <- as.matrix(dist(landscape))
  }
  
  disp_array <- array(dim = c(patches, patches, species))
  for(k in 1:species){
    spec_dist_mat <- exp(-species_traits[k,"kernel_exp"] * dist_mat)
    # next, make all cols sum to 1
    disp_array[,,k] <- apply(spec_dist_mat, 1, function(x) x / sum(x))
    if (sum(colSums(disp_array[,,k]) > 1.001) > 0) warning (
      "dispersal from a patch to all others exceeds 100%. 
      Make sure the rowSums(disp_mat) <= 1")
    if (sum(colSums(disp_array[,,k]) < 0.999) > 0) warning (
      "dispersal from a patch to all others is less than 100%. 
      Some dispersing individuals will be lost from the metacommunity")
  }
  
  return(disp_array)
  
}

compute_r_xt <- function(species_traits = species_traits, env = env, species = species){
  
  # get env matrix at time t
  env_mat <- matrix(rep(env, each = species), nrow = species, ncol = patches)
  
  env_mismatch <- exp(-((species_traits$env_niche_optima - env_mat) / (2*species_traits$env_niche_breadth))^2)
  
  r_ixt <- species_traits$max_r*t(env_mismatch)
  return(r_ixt)
}

survival <- function(N, species_traits){
  s_prop <- species_traits$survival * (1-species_traits$germ)
  if(nrow(species_traits) != ncol(N)) {stop("Dimensions off")}
  
  N_surv <- N * 0
  for(i in 1:ncol(N)){
    #N_surv[,i] <- N[,i] * s_prop[i]
    N_surv[,i] <- rbinom(n = nrow(N), size = N[,i], prob = s_prop[i])
  }
  return(N_surv)
}

germination <- function(N, species_traits){
  g_prop <- species_traits$germ
  if(nrow(species_traits) != ncol(N)) {stop("Dimensions off")}
  
  N_germ <- N * 0
  for(i in 1:ncol(N)){
    #N_surv[,i] <- N[,i] * s_prop[i]
    N_germ[,i] <- rbinom(n = nrow(N), size = N[,i], prob = g_prop[i])
  }
  return(N_germ)
}

growth <- function(N, species_traits, r, int_mat){
  N_growth <- N*0
  #germ <- matrix(species_traits$germ, nrow = 1)
  
  for(i in 1:ncol(N)){
    # compute number of seeds that germinate
    #N_germ <- rbinom(n = nrow(N), size = N[,i], prob = germ[i])
    # of germinated seeds, grow
    N_growth[,i] <- N[,i] * r[,i]
  }
  N_next <- N_growth / (1 + N%*%int_mat)
  return(N_next)
}


get_comp_effects <- function(N, species_traits, r, int_mat){
  N_growth <- N*0
  
  for(i in 1:ncol(N)){
    N_growth[,i] <- N[,i] * r[,i]
  }
  
  # find inter comp effects on growth
  int_mat_comp <- int_mat
  diag(int_mat_comp) <- 0
  #image(int_mat_comp)
  N_next_inter_only <- N_growth / (1 + N%*%int_mat_comp)
  #N_next_inter_only
  
  # find intra comp effects on growth
  int_mat_comp <- int_mat
  int_mat_comp[upper.tri(int_mat_comp, diag = FALSE)] <- 0
  int_mat_comp[lower.tri(int_mat_comp, diag = FALSE)] <- 0
  #image(int_mat_comp)
  N_next_intra_only <- N_growth / (1 + N%*%int_mat_comp)
  #N_next_intra_only
  
  # find full comp effects on growth
  int_mat_comp <- int_mat
  #image(int_mat_comp)
  N_next_full <- N_growth / (1 + N%*%int_mat_comp)
  #N_next_full
  
  # find effects of no comp on growth
  int_mat_comp <- int_mat * 0
  #image(int_mat_comp)
  N_next_nocomp <- N_growth / (1 + N%*%int_mat_comp)
  #N_next_nocomp
  
  comp_effects <- list()
  comp_effects$full <- N_next_full
  comp_effects$intra <- N_next_intra_only
  comp_effects$inter <- N_next_inter_only
  comp_effects$nocomp <- N_next_nocomp
  
  return(comp_effects)
}

# # THompson model
# env_generate <- function(landscape, env.df, env1Scale = 500, timesteps = 1000, plot = TRUE){
#   if (missing(env.df)){
#     repeat {
#       model <- RandomFields::RMexp(var=0.5, scale=env1Scale) + # with variance 4 and scale 10
#         RandomFields::RMnugget(var=0) + # nugget
#         RandomFields::RMtrend(mean=0.05) # and mean
#       
#       RF <- RandomFields::RFsimulate(model = model,x = landscape$x*10, y = landscape$y*10, T = 1:timesteps, spConform=FALSE)
#       env.df <- data.frame(env1 = vegan::decostand(RF,method = "range"), patch = 1:nrow(landscape), time = rep(1:timesteps, each = nrow(landscape)))
#       env.initial <- env.df[env.df$time == 1,]
#       if((max(env.initial$env1)-min(env.initial$env1)) > 0.6) {break}
#     }
#   } else {
#     if(all.equal(names(env.df), c("env1", "patch", "time")) != TRUE) stop("env.df must be a dataframe with columns: env1, patch, time")
#   }
#   
#   if(plot == TRUE){
#     g<-ggplot2::ggplot(env.df, aes(x = time, y = env1, group = patch, color = factor(patch)))+
#       ggplot2::geom_line()+
#       scale_color_viridis_d(guide=F)
#     print(g)
#   }
#   
#   return(env.df)
# }

species_int_mat <- function(species, intra = 1, min_inter = 0, max_inter = 1.5, int_matrix, comp_scaler = 0.05, plot = TRUE){
  if (missing(int_matrix)){
    int_mat <- matrix(runif(n = species*species, min = min_inter, max = max_inter), nrow = species, ncol = species)
    diag(int_mat) <- intra
    int_mat <- int_mat * comp_scaler
  } else {
    if (is.matrix(int_matrix) == FALSE) stop("int_matrix must be a matrix")
    if (dim(int_matrix) != c(species,species)) stop("int_matrix must be a matrix with a row and column for each species")
    if (is.numeric(int_matrix) == FALSE) stop("int_matrix must be numeric")
  }
  
  if (plot == TRUE){
    colnames(int_mat)<- 1:species
    g <- as.data.frame(int_mat) %>%
      dplyr::mutate(i = 1:species) %>%
      tidyr::gather(key = j, value = competition, -i) %>%
      dplyr::mutate(i = as.numeric(as.character(i)),
                    j = as.numeric(as.character(j))) %>%
      ggplot2::ggplot(ggplot2::aes(x = i, y = j, fill = competition))+
      ggplot2::geom_tile()+
      scale_fill_viridis_c(option = "E")
    
    print(g)
  }
  return(int_mat)
}


# new environmental functions:
generate_noise_ts <- function(a, length, sd = 1){
  sd = sd
  sig_vec <- NULL
  sig_vec[1] <- 0
  a <- a
  b <- (1-a^2)^0.5
  for(t in 2:length){
    sig_vec[t] <- a*sig_vec[t-1] + b * rnorm(n = 1, mean = 0, sd = sd)
  }
  return(sig_vec)
}


# simpler env simulator
env_generate <- function(landscape, 
                         spat_heterogeneity = 0.5, 
                         temp_noise_color = 0,
                         temp_noise_sd = 0,
                         timesteps = 1000, A, w, phi){
  
  env_mat <- matrix(0, nrow = timesteps, ncol = nrow(landscape))
  
  # spatial heterogeneity is from linear increase in x and y directions
  gen_spatial_var <- function(x, y) return(spat_heterogeneity*x + spat_heterogeneity*y)
  env_mat_init <- gen_spatial_var(landscape$x, landscape$y)
  for(patch in 1:ncol(env_mat)){
    env_mat[,patch] = env_mat[,patch] + env_mat_init[patch]
  }
  
  # sine wave params
  A <- A # amplitude, peak deviation
  w <- w # angular frequency, oscillations per time interval in radians per second
  phi <- phi # phase, in radians where the cycle is at t = 0
  
  # add one timeseries on top of starting conditions
  sim_ts_trend <- A*sin(w*(1:timesteps) + phi)
  sim_ts_trend_noise <- sim_ts_trend + generate_noise_ts(a = temp_noise_color, length = nrow(env_mat), sd = temp_noise_sd)
  
  
  env_mat <- env_mat + sim_ts_trend_noise
  
  # spec.mtm(generate_noise_ts(a = 0, length = 10000, sd = temp_noise_sd), dtUnits = "year")
  # spectrum(env_ts)
  # 
  
  env_mat <- (env_mat - min(env_mat)) / (max(env_mat) - min(env_mat))
  # matplot(env_mat, type = "l")
  
  
  env_df <- tidyr::pivot_longer(cbind.data.frame(time = 1:timesteps, env_mat), cols = (1:ncol(env_mat)+1), names_to = "patch", values_to = "env")
  # env_df %>% 
  #   ggplot(aes(x = time, y= env, color = patch)) + geom_line(alpha = 0.2, show.legend = F) + theme_bw() + scale_color_viridis_d()
  
  return(env_df)
}


# env_generate <- function(landscape, x_dim, y_dim, spat_auto = 0.5, temp_auto = 0, timesteps = 1000, A, k, w, phi){
#   grid <- list(x = seq(0, 100, length.out = 100), y = seq(0, 100, length.out = 100)) 
#   
#   repeat{
#     obj <-fields::Exp.image.cov(grid = grid, theta=spat_auto, setup=TRUE)
#     look <- fields::sim.rf(obj)
#     image.plot( grid$x, grid$y, look) 
#     points(landscape)
#     # title("simulated gaussian field")
#     
#     sig_mat <- matrix(NA, nrow = timesteps, ncol = nrow(landscape))
#     sig_mat[1,] <- look[cbind(landscape$x, landscape$y)]
#     
#     # sine wave params
#     A <- A
#     k <- k # wave number, wavelengths per unit distance
#     w <- w # angular frequency, relates to speed of propagation, v = w/k
#     phi <- phi # phase, where in the cycle oscillation is at t=0  
#     
#     for(patch in 1:nrow(landscape)){
#       mean_cond <- sig_mat[1,patch]
#       with_trend <- mean_cond + A*sin(k*landscape$x[patch] - w*(1:timesteps)) + A*sin(k*landscape$y[patch] - w*(1:timesteps)) + phi
#       with_noise <- with_trend + generate_noise_ts(a = temp_auto, length = nrow(sig_mat), sd = 0.1)
#       sig_mat[,patch] <- with_noise
#       #plot(with_noise, type = 'l')
#       #spec.mtm(with_noise)
#       
#     }
#     sig_mat <- (sig_mat - min(sig_mat)) / (max(sig_mat) - min(sig_mat))
#     if(max(sig_mat[1,]) - min(sig_mat[1,]) > 0.6){break}
#   }
#   sig_df <- tidyr::pivot_longer(cbind.data.frame(time = 1:timesteps, sig_mat), cols = (1:ncol(sig_mat)+1), names_to = "patch", values_to = "env")
#   sig_df %>% 
#     filter(as.numeric(patch) < 10) %>% 
#     ggplot(aes(x = time, y= env, color = patch)) + geom_line(alpha = 0.8, show.legend = F) + theme_bw()
#   
#   return(sig_df)
# }

# anim <- sig_df %>% left_join(rownames_to_column(landscape, var = "patch")) %>% 
#   ggplot() +
#   geom_tile(mapping = aes(x = x, y = y, fill = env)) +
#   transition_states(states = time) +
#   theme_minimal() +
#   scale_fill_viridis()
# animate(anim, duration = 10)  

# some error handling
is_simple_error <- function(x) inherits(x, "simpleError")
is_try_error <- function(x) inherits(x, "try-error")


# geometric mean 
gm_mean = function(x, na.rm=TRUE, zero.propagate = FALSE){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
}

