library(igraph)
source("simulations/metacom_functions.R")
set.seed(10294)
kernel <- function(x, k){
  return(10^((-k) * x))
}

# Figure S2: Visualize connectivity for a landscape under
# differing dispersal kernel exponents

landscape <- init_landscape(patches = 100, x_dim = 100, y_dim = 100)

dist_matrix <- as.matrix(dist(landscape))
kernel_vals <- c(0, 10^seq(-4, 0, length.out = 9))
for (k in kernel_vals){
  
  disp_matrix <- exp(-k * dist_matrix)
  disp_matrix <- apply(disp_matrix, 1, function(x) x / sum(x))
  
  connectivity <- graph_from_adjacency_matrix(disp_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)
  
  png(paste0("figures/FigS2/connectivity_Li",k,".png"), width = 6, height = 6, units = "in", res = 500)
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
              edge.width = 5*E(connectivity)$weight,
              alpha = .5
  )
  dev.off()
  
}


# make 4 and then assemble
plot_connectivity <- function(k){
  disp_matrix <- exp(-k * dist_matrix)
  disp_matrix <- apply(disp_matrix, 1, function(x) x / sum(x))
  
  connectivity <- graph_from_adjacency_matrix(disp_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)
  
  fig = plot.igraph(connectivity, 
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
              edge.width = 5*E(connectivity)$weight,
              alpha = .5
  )
  return(fig)
  
}


png("figures/FigS2/connectivity_Li_multi.png", width = 6, height = 6, units = "in", res = 500)

par(mfrow = c(2, 2),
    mar = c(4, 4, 3, 1))  # give a bit of top margin for labels

idx <- c(1, 4, 6, 8)
labels <- c("A", "B", "C", "D")

for (j in seq_along(idx)) {
  i <- idx[j]
  
  plot_connectivity(k = kernel_vals[i])
  
  mtext(
    paste0(labels[j], ")  Li = ", kernel_vals[i]),
    side = 3,      # top
    line = 1,      # distance from plot
    adj = 0,       # left-align
    cex = 0.9,
    font = 2       # bold
  )
}

dev.off()
