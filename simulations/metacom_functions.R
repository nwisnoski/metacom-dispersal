# Set of functions related to metacommunity simulations run in 
# the "metacom_dispersal-kernel.R" simulations


# dist.torus function from the package som.nn
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


# initialize landscape to form non-overlapping patches
init_landscape <- function(patches, x_dim = 100, y_dim = 100){
  repeat{
    landscape <- data.frame(x = sample(1:x_dim, size = patches, replace = T),
                            y = sample(1:y_dim, size = patches, replace = T))
    if(dim(unique(landscape))[1] == patches) {break}
  }
  
  # plot(landscape)
  
  return(landscape)
}


# create species with given trait distributions according to parameters
# specified in the main simulation script 
# modified from Thompson et al. 2020

init_species <- function(species = 10, 
                         env_niche_optima = "random",
                         env_niche_breadth = 0.5,
                         kernel_exp = 0.1,
                         max_r = 5,
                         dispersal_rate = 0.1,
                         survival = 0,
                         germ = 1){
  
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
  
  # generate dormancy rates
  if(length(germ) != species & length(germ) != 1){
    stop("Enter a number or a vector of length 'species' for germ.")
  }
  if(length(survival) != species & length(survival) != 1){
    stop("Enter a number or a vector of length 'species' for survival.")
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
  
  return(species_traits)

}

# initialize community with poisson-distributed colonization process
init_community <- function(initialization = 200, species = 10, patches = 100){
  
  N <- matrix(rpois(n = species * patches, lambda = 0.5), nrow = patches, ncol = species)
  
  return(N)
}

# compute dispersal matrices for each species
# stored as an array with different slices for each species
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

# local growth, gaussian fit to environment
compute_r_xt <- function(species_traits = species_traits, env = env, species = species){
  
  # get env matrix at time t
  env_mat <- matrix(rep(env, each = species), nrow = species, ncol = patches)
  
  env_mismatch <- exp(-((species_traits$env_niche_optima - env_mat) / (2*species_traits$env_niche_breadth))^2)
  
  r_ixt <- species_traits$max_r*t(env_mismatch)
  return(r_ixt)
}



# competition function
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


# calculations for competition fitness effects
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

# compute interaction strengths for competition
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


# environmental noise process
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


# simulate 2D linear environment with temporal fluctuations
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
  
  
  
  env_mat <- (env_mat - min(env_mat)) / (max(env_mat) - min(env_mat))
  
  
  env_df <- tidyr::pivot_longer(cbind.data.frame(time = 1:timesteps, env_mat), cols = (1:ncol(env_mat)+1), names_to = "patch", values_to = "env")
  
  return(env_df)
}



# some error handling
is_simple_error <- function(x) inherits(x, "simpleError")
is_try_error <- function(x) inherits(x, "try-error")




