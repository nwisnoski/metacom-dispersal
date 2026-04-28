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

# define a zero-abundance-safe hellinger function for sites/times with total abund = 0
safe_hellinger <- function(mat) {
  rs <- rowSums(mat, na.rm = TRUE)
  out <- mat
  out[rs > 0, ] <- sqrt(out[rs > 0, , drop = FALSE] / rs[rs > 0])
  out[rs == 0, ] <- 0
  out
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
    sbs <- safe_hellinger(sbs)
    tot_ss <- sum((scale(sbs, center = TRUE, scale = FALSE))^2)
    beta_div <- tot_ss / (nrow(sbs)-1)
    
    spat_beta_dynamics[time] <- beta_div
  }
  
  # how does temporal beta diversity vary across space
  temp_beta_distribution <- numeric(length = ts_dims["sites"])
  for(site in 1:ts_dims["sites"]){
    site_ts <- metacomm_tsdata[,,site]
    tbs <- t(site_ts)
    tbs <- safe_hellinger(tbs)
    tot_ss <- sum((scale(tbs, center = TRUE, scale = FALSE))^2)
    beta_div <- tot_ss / (nrow(tbs)-1)
    
    temp_beta_distribution[site] <- beta_div
  }
  
  beta_div <- data.frame(mean_beta_spatial = mean(spat_beta_dynamics),
                         mean_beta_temporal = mean(temp_beta_distribution))
  return(beta_div)
}



