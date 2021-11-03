library(tidyverse)
library(here)

equal <- read_csv(here("sim_output/equal_2021-11-03_162607.csv"))
stable <- read_csv(here("sim_output/stable_2021-11-03_143719.csv"))
priority <- read_csv(here("sim_output/priority_2021-11-03_151300.csv"))

equal_sub <- equal[1:100000,]

alpha_div <- equal %>% 
  group_by(comp, rep, dispersal, kernel_exp, extirp_prob, patch) %>% 
  summarize(alpha = sum(1*(N>0))) %>% 
  ungroup() %>% 
  group_by(comp,rep,dispersal,kernel_exp,extirp_prob) %>% 
  summarize(mean_alpha = mean(alpha))

gamma_div <- equal %>% 
  group_by(comp, rep, dispersal, kernel_exp, extirp_prob) %>% 
  filter(N > 0) %>% 
  summarize(gamma = length(unique(species)))

#### 
alpha_fig <- alpha_div %>% 
  ggplot(aes(x = dispersal, y = kernel_exp, z = mean_alpha)) + 
  geom_contour_filled() + 
  facet_wrap(~extirp_prob) +
  theme_bw()

gamma_fig <- gamma_div %>% 
  ggplot(aes(x = dispersal, y = kernel_exp, z = gamma)) + 
  geom_contour_filled() + 
  facet_wrap(~extirp_prob) +
  theme_bw()

beta_fig <- alpha_div %>% 
  full_join(gamma_div) %>% 
  mutate(beta = gamma / mean_alpha) %>% 
  ggplot(aes(x = dispersal, y = kernel_exp, z = beta)) + 
  geom_contour_filled() +
  facet_wrap(~extirp_prob) + 
  theme_bw()

ggsave(filename = here("figures/alpha.pdf"), plot = alpha_fig, width = 8, height = 6)
ggsave(filename = here("figures/beta.pdf"), plot = beta_fig, width = 8, height = 6)
ggsave(filename = here("figures/gamma.pdf"), plot = gamma_fig, width = 8, height = 6)
