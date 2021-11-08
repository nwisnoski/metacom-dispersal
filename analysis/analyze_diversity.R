library(tidyverse)
library(here)


# subset by dates
files <- dir(here("sim_output/")) %>% 
  tibble() %>% rename("file" = ".") %>% 
  filter(!str_detect(string = file, pattern = "2021-11-03"))

# subset by competition treatments
files_equal <- files %>% 
  filter(str_detect(string = file, pattern = "equal_"))
files_stable <- files %>% 
  filter(str_detect(string = file, pattern = "stable_"))
files_priority <- files %>% 
  filter(str_detect(string = file, pattern = "priority_"))

equal <- NULL
stable <- NULL
priority <- NULL

for (f in files_equal$file) {
  equal.f <- read_csv(here(paste0("sim_output/",f)))
  equal <- bind_rows(equal, equal.f)
  
}

for (f in files_stable$file) {
  stable.f <- read_csv(here(paste0("sim_output/",f)))
  stable <- bind_rows(stable, stable.f)
  
}

for (f in files_priority$file) {
  priority.f <- read_csv(here(paste0("sim_output/",f)))
  priority <- bind_rows(priority, priority.f)
  
}
# equal <- read_csv(here("sim_output/equal_2021-11-05"))
# stable <- read_csv(here("sim_output/stable_2021-11-03_170641.csv"))
# priority <- read_csv(here("sim_output/priority_2021-11-03_174448.csv"))


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
    group_by(comp, rep, dispersal, kernel_exp, extirp_prob, species) %>% 
    summarize(total_abund = sum(N)) %>% 
    group_by(comp, rep, dispersal, kernel_exp, extirp_prob) %>% 
    summarize(gamma = sum(1*(total_abund>0)))
  return(gamma_div)
}

compute_beta <- function(alpha_div, gamma_div){
  beta_div <- full_join(alpha_div, gamma_div) %>% 
    mutate(beta = gamma / mean_alpha,
           beta = ifelse(is.nan(beta), 0, beta))
  
  return(beta_div)
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

make_beta_fig <- function(beta_div){
  beta_fig <- beta_div %>% 
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

equal_beta <- compute_beta(equal_alpha, equal_gamma)
stable_beta <- compute_beta(stable_alpha, stable_gamma)
priority_beta <- compute_beta(priority_alpha, priority_gamma)

rm(list = c("equal", "stable", "priority"))

fig_equal_alpha <- make_alpha_fig(equal_alpha)
fig_equal_beta <- make_beta_fig(equal_beta)
fig_equal_gamma <- make_gamma_fig(equal_gamma)

fig_stable_alpha <- make_alpha_fig(stable_alpha)
fig_stable_beta <- make_beta_fig(stable_beta)
fig_stable_gamma <- make_gamma_fig(stable_gamma)

fig_priority_alpha <- make_alpha_fig(priority_alpha)
fig_priority_beta <- make_beta_fig(priority_beta)
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
