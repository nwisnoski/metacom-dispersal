library(tidyverse)
library(here)

equal <- read_csv(here("sim_output/equal_2021-11-03_162607.csv"))
stable <- read_csv(here("sim_output/stable_2021-11-03_170641.csv"))
priority <- read_csv(here("sim_output/priority_2021-11-03_174448.csv"))


compute_alpha <- function(sim_dat){
  alpha_div <- sim_dat %>% 
    group_by(comp, rep, dispersal, kernel_exp, extirp_prob, patch) %>% 
    summarize(alpha = sum(1*(N>0))) %>% 
    ungroup() %>% 
    group_by(comp,rep,dispersal,kernel_exp,extirp_prob) %>% 
    summarize(mean_alpha = mean(alpha))
  return(alpha_div)
  
}



compute_gamma <- function(sim_dat){
  gamma_div <- sim_dat %>% 
    group_by(comp, rep, dispersal, kernel_exp, extirp_prob) %>% 
    filter(N > 0) %>% 
    summarize(gamma = length(unique(species)))
  return(gamma_div)
}


#### 
make_alpha_fig <- function(alpha_div){
  alpha_fig <- alpha_div %>% 
    ggplot(aes(x = dispersal, y = kernel_exp, z = mean_alpha)) + 
    geom_contour_filled() + 
    facet_wrap(~extirp_prob) +
    theme_bw()
  return(alpha_fig)
}

make_gamma_fig <- function(gamma_div){
  gamma_fig <- gamma_div %>% 
    ggplot(aes(x = dispersal, y = kernel_exp, z = gamma)) + 
    geom_contour_filled() + 
    facet_wrap(~extirp_prob) +
    theme_bw()
  return(gamma_fig)
}

make_beta_fig <- function(alpha_div, gamma_div){
  beta_fig <- alpha_div %>% 
    full_join(gamma_div) %>% 
    mutate(beta = gamma / mean_alpha) %>% 
    ggplot(aes(x = dispersal, y = kernel_exp, z = beta)) + 
    geom_contour_filled() +
    facet_wrap(~extirp_prob) + 
    theme_bw()
  return(beta_fig)
}


# apply functions to data
equal_alpha <- compute_alpha(equal)
stable_alpha <- compute_alpha(stable)
priority_alpha <- compute_alpha(priority)

equal_gamma <- compute_gamma(equal)
stable_gamma <- compute_gamma(stable)
priority_gamma <- compute_gamma(priority)

fig_equal_alpha <- make_alpha_fig(equal_alpha)
fig_equal_beta <- make_beta_fig(equal_alpha, equal_gamma)
fig_equal_gamma <- make_gamma_fig(equal_gamma)

fig_stable_alpha <- make_alpha_fig(stable_alpha)
fig_stable_beta <- make_beta_fig(stable_alpha, stable_gamma)
fig_stable_gamma <- make_gamma_fig(stable_gamma)

fig_priority_alpha <- make_alpha_fig(priority_alpha)
fig_priority_beta <- make_beta_fig(priority_alpha, priority_gamma)
fig_priority_gamma <- make_gamma_fig(priority_gamma)

# save figures
ggsave(filename = here("figures/equal_alpha.pdf"), plot = fig_equal_alpha, width = 8, height = 6)
ggsave(filename = here("figures/equal_beta.pdf"), plot = fig_equal_beta, width = 8, height = 6)
ggsave(filename = here("figures/equal_gamma.pdf"), plot = fig_equal_gamma, width = 8, height = 6)

ggsave(filename = here("figures/stable_alpha.pdf"), plot = fig_stable_alpha, width = 8, height = 6)
ggsave(filename = here("figures/stable_beta.pdf"), plot = fig_stable_beta, width = 8, height = 6)
ggsave(filename = here("figures/stable_gamma.pdf"), plot = fig_stable_gamma, width = 8, height = 6)

ggsave(filename = here("figures/priority_alpha.pdf"), plot = fig_priority_alpha, width = 8, height = 6)
ggsave(filename = here("figures/priority_beta.pdf"), plot = fig_priority_beta, width = 8, height = 6)
ggsave(filename = here("figures/priority_gamma.pdf"), plot = fig_priority_gamma, width = 8, height = 6)
