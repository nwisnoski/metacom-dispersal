
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


