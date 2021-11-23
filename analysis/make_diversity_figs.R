require(tidyverse)
require(here)

diversity_full <- NULL
for(f in c("equal", "priority", "stable")){
  alpha <- read_csv(file = here("data",paste0(f, "_alpha.csv")))
  beta <- read_csv(file = here("data",paste0(f, "_beta.csv")))
  gamma <- read_csv(file = here("data",paste0(f, "_gamma.csv")))
  tot <- left_join(alpha, beta) %>% left_join(gamma)
  diversity_full <- bind_rows(diversity_full, tot)
}


theme_set(theme_minimal())

diversity <- diversity_full %>% mutate(beta = (gamma - mean_alpha)/gamma) %>% 
  mutate(beta = ifelse(is.nan(beta), 0, beta)) %>% 
  filter(extirp_prob %in% c(0, 0.05))

div_equal_fig <- diversity %>% 
  filter(comp == "equal") %>% 
  pivot_longer(cols = c("mean_alpha", "beta", "gamma"), names_to = "div_level", values_to = "diversity") %>% 
  mutate(div_level = factor(div_level, levels = c("mean_alpha", "beta", "gamma"))) %>% 
  
  ggplot(aes(x = dispersal, y = diversity, group = interaction(kernel_exp), color = kernel_exp)) + 
  facet_grid(div_level~extirp_prob, scales = "free") + 
  geom_point(alpha = 0.01) + 
  geom_line(stat = "smooth", method = "loess", span = 0.75, alpha = 0.75, size = 1, se = FALSE) +
  scale_x_log10() + 
  scale_color_viridis_c() +
  labs(title = "Equal competition")



div_stable_fig <- diversity %>% 
  filter(comp == "stable") %>% 
  pivot_longer(cols = c("mean_alpha", "beta", "gamma"), names_to = "div_level", values_to = "diversity") %>% 
  mutate(div_level = factor(div_level, levels = c("mean_alpha", "beta", "gamma"))) %>% 
  
  ggplot(aes(x = dispersal, y = diversity, group = interaction(kernel_exp), color = kernel_exp)) + 
  facet_grid(div_level~extirp_prob, scales = "free") + 
  geom_point(alpha = 0.01) + 
  geom_line(stat = "smooth", method = "loess", span = 0.75, alpha = 0.75, size = 1, se = FALSE) +
  scale_x_log10() + 
  scale_color_viridis_c() +
  labs(title = "Stabilizing competition")


div_prior_fig <- diversity %>% 
  filter(comp == "priority") %>% 
  pivot_longer(cols = c("mean_alpha", "beta", "gamma"), names_to = "div_level", values_to = "diversity") %>% 
  mutate(div_level = factor(div_level, levels = c("mean_alpha", "beta", "gamma"))) %>% 
  
  ggplot(aes(x = dispersal, y = diversity, group = interaction(kernel_exp), color = kernel_exp)) + 
  facet_grid(div_level~extirp_prob, scales = "free") + 
  geom_point(alpha = 0.01) + 
  geom_line(stat = "smooth", method = "loess", span = 0.75, alpha = 0.75, size = 1, se = FALSE) +
  scale_x_log10() + 
  scale_color_viridis_c() +
  labs(title = "Priority effects")

# view
ggsave(filename = "figures/diversity_stable.pdf", plot = div_stable_fig, width = 10, height = 8)
ggsave(filename = "figures/diversity_equal.pdf", plot = div_equal_fig, width = 10, height = 8)
ggsave(filename = "figures/diversity_priority.pdf", plot = div_prior_fig, width = 10, height = 8)

