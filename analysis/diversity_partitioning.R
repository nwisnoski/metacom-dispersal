# partition alpha, beta, gamma div
diversity.partition <- function(metacomm_tsdata){
  ts_dims <- dim(metacomm_tsdata)
  names(ts_dims) <- c("species", "time", "sites")
  
  div_dynamics <- matrix(ncol = 3, nrow = ts_dims["time"])
  
  for(time in 1:ts_dims["time"]){
    snapshot <- metacomm_tsdata[,time,]
    
    occupancy <- rowSums((snapshot > 0) * 1)
    gamma_div <- sum(occupancy > 0)
    
    alpha_div <- mean(colSums(snapshot > 0))
    
    beta_div <- gamma_div / alpha_div
    
    div_dynamics[time,] <- c(alpha_div, beta_div, gamma_div)
  }
  
  return(colMeans(div_dynamics))
}



# compute spatial and temporal beta
beta.div.calc <- function(metacomm_tsdata){
  ts_dims <- dim(metacomm_tsdata)
  names(ts_dims) <- c("species", "time", "sites")
  
  # how does spatial beta diversity change over time
  spat_beta_dynamics <- numeric(length = ts_dims["time"])
  
  for(time in 1:ts_dims["time"]){
    snapshot <- metacomm_tsdata[,time,]
    sbs <- t(snapshot)
    sbs <- vegan::decostand(sbs, "hellinger")
    tot_ss <- sum((scale(sbs, center = TRUE, scale = FALSE))^2)
    beta_div <- tot_ss / (nrow(sbs)-1)
    
    spat_beta_dynamics[time] <- beta_div
  }
  
  # how does temporal beta diversity vary across space
  temp_beta_distribution <- numeric(length = ts_dims["sites"])
  for(site in 1:ts_dims["sites"]){
    site_ts <- metacomm_tsdata[,,site]
    tbs <- t(site_ts)
    tbs <- vegan::decostand(tbs, "hellinger")
    tot_ss <- sum((scale(tbs, center = TRUE, scale = FALSE))^2)
    beta_div <- tot_ss / (nrow(tbs)-1)
    
    temp_beta_distribution[site] <- beta_div
  }
  
  beta_div <- data.frame(mean_beta_spatial = mean(spat_beta_dynamics),
                         mean_beta_temporal = mean(temp_beta_distribution))
  return(beta_div)
}



# compute patch-level diversity and stability
div.stab.comp <- function(metacomm_tsdata){
  ts_dims <- dim(metacomm_tsdata)
  names(ts_dims) <- c("species", "time", "sites")
  
  
  # how does temporal beta diversity vary across space
  temp_alpha_var_distribution <- matrix(ncol = 2, nrow = ts_dims["sites"])
  
  for(site in 1:ts_dims["sites"]){
    site_ts <- t(metacomm_tsdata[,,site])
    
    site_alpha <- mean(rowSums((site_ts > 0) * 1))
    site_cv <- sd(rowSums(site_ts)) / mean(rowSums(site_ts))
    
    temp_alpha_var_distribution[site,1] <- site_alpha
    temp_alpha_var_distribution[site,2] <- site_cv
  }
  
  div_stab_df <- as.data.frame(temp_alpha_var_distribution)
  names(div_stab_df) <- c("alpha", "cv")
  
  return(data.frame(local_mean_richness = mean(div_stab_df$alpha), 
             local_mean_cv = mean(div_stab_df$cv)))
}

