library(tidyverse)
source("simulations/metacom_functions.R")

kernel <- function(x, k){
  return(10^((-k) * x))
}

landscape <- init_landscape(patches = 100, x_dim = 100, y_dim = 100)

hist(dist(landscape))

dist_matrix <- as.matrix(dist(landscape))
kernel_vals <- c(0, 10^seq(-4, 0, length.out = 9))
for (k in kernel_vals){
  
  disp_matrix <- exp(-k * dist_matrix)
  disp_matrix <- apply(disp_matrix, 1, function(x) x / sum(x))
  
  dispersers <- data.frame(disperers = sample(dist_matrix, size = 10000, replace = T, prob = disp_matrix))
  available <- data.frame(landscape = as.numeric(as.matrix(dist(landscape))))
  disperser_mean <- signif(mean(dispersers$disperers), 3)
  landscape_mean <- signif(mean(available$landscape), 3)
  
  disp_landscape_fig <- 
    cbind.data.frame(dispersers, available) |> 
    pivot_longer(cols = everything(), names_to = "subset", values_to = "distance") |> 
    ggplot(aes(x = distance, color = subset)) + 
    geom_density( alpha = 0.25, linewidth = 1) +
    geom_density(alpha = 0.25, linewidth = 1) +
    theme_minimal() +
    labs(x = "Distance", y = "Density", color = "",
         title = paste("Dispersal kernel exponent, k =", signif(k,3)),
         subtitle = paste0("Mean inter-patch distance = ", landscape_mean,"; mean dispersal distance = ",disperser_mean))
  
  ggsave(filename = paste0("figures/disp_landscape_mismatch_k",k,".png"), width = 6, height = 6, dpi = 500)
  
}



