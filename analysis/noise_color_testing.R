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



generate_noise_ts <- function(a, length, ...){
  sig_vec <- NULL
  sig_vec[1] <- 0
  a <- a
  b <- (1-a^2)^0.5
  for(t in 2:length){
    sig_vec[t] <- a*sig_vec[t-1] + b * rnorm(n = 1, mean = 0, sd = 1)
  }
  return(sig_vec)
}
# sig_vec <- generate_noise_ts(a=0.2, length = 1000)
# plot(sig_vec,type='l')
# spectrum(as.ts(sig_vec), method = "pgram")
# sig_spec <- multitaper::spec.mtm(as.ts(sig_vec))
# tibble(freq = sig_spec$freq, spec = sig_spec$spec) %>% 
#   ggplot(aes(x = freq, y = spec)) + 
#   geom_line() +
#   geom_smooth(method = "lm", se = F) +
#   scale_y_log10(labels = scales::label_number()) + 
#   scale_x_log10(labels = scales::label_number()) +
#   labs(x = "Frequency (/time)", y = "Spectral density")


#Simulate a Gaussian random field with an exponential covariance function,  
#range parameter = 2.0 and the domain is  [0,5]X [0,5] evaluating the 
#field at a 100X100 grid.  

env_generate <- function(landscape, x_dim, y_dim, theta = 0.5, temp_auto = 0, timesteps = 1000){
  grid <- list(x = seq(0, 5,, x_dim), y = seq(0, 5,, y_dim)) 
  obj <-fields::Exp.image.cov(grid = grid, theta=theta, setup=TRUE)
  look <- fields::sim.rf(obj)
  # image.plot( grid$x, grid$y, look) 
  # title("simulated gaussian field")
  
  sig_mat <- matrix(NA, nrow = timesteps, ncol = nrow(landscape))
  sig_mat[1,] <- look[cbind(landscape$x, landscape$y)]
  for(patch in 1:nrow(landscape)){
    mean_cond <- sig_mat[1,patch]
    with_noise <- mean_cond + generate_noise_ts(a = temp_auto, length = nrow(sig_mat))
    sig_mat[,patch] <- with_noise
    #plot(with_noise, type = 'l')
    #spec.mtm(with_noise)
    
  }
  sig_mat <- (sig_mat - min(sig_mat)) / (max(sig_mat) - min(sig_mat))
  return(sig_mat)
}

#matplot(env_df[1:100, 50:100], type = 'l')
