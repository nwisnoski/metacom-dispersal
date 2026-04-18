library(igraph)
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
  
  connectivity <- graph_from_adjacency_matrix(disp_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)
  
  png(paste0("figures/connectivity_k",k,".png"), width = 6, height = 6, units = "in", res = 500)
  plot.igraph(connectivity, 
              layout = as.matrix(landscape), 
              rescale = FALSE,
              xlim = c(0,100), ylim = c(0,100),
              axes = TRUE,
              vertex.label = NA,
              vertex.size = 400, 
              vertex.label.cex = 0.5,
              vertex.color = NA,
              edge.color = "gray40",
              edge.weight = E(connectivity)$weight,
              edge.width = 10*E(connectivity)$weight,
              alpha = .5
  )
  dev.off()
  
}


