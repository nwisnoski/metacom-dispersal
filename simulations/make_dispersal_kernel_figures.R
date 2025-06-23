kernel <- function(x, k){
  return(10^((-k) * x))
}

pdf("figures/dispersal_kernel_family.pdf", width = 4, height = 4)
curve(expr = kernel(x, k=1), from = 0, to = 10, n = 100, type = "l",xlab = "Distance", ylab = "Probability of Dispersal")
for(k in 0:10/10){
  curve(expr = kernel(x, k=k), from = 0, to = 10, n = 100, type = "l", add = TRUE)
  
}
dev.off()

pdf("figures/dispersal_kernel_0.9.pdf", width = 4, height = 4)
curve(expr = kernel(x, k=0.9), from = 0, to = 100, n = 100, type = "l",xlab = "Distance", ylab = "Probability of Dispersal")
dev.off()

pdf("figures/dispersal_kernel_0.1.pdf", width = 4, height = 4)
curve(expr = kernel(x, k=0.1), 
      from = 0, to = 100, n = 100, 
      type = "l", xlab = "Distance", 
      ylab = "Probability of Dispersal", ylim = c(0,1))
dev.off()


k = 0.25
pdf(paste0("figures/dispersal_kernel_",k,".pdf"), width = 4, height = 4)
curve(expr = kernel(x, k=k), 
      from = 0, to = 100, n = 100, 
      type = "l", xlab = "Distance", 
      ylab = "Probability of Dispersal", ylim = c(0,1))
dev.off()
