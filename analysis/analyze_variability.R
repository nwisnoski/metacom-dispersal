library(tidyverse)
library(here)
library(patchwork)

theme_set(theme_bw())

variability <- read_csv(here("data/variability_partitioning_2022-11-14_112210.csv"))

kernel_exps <- unique(variability$kernel_exp)
disp_rates <- unique(variability$disp_rate)
disturb_rates <- unique(variability$disturb_rate)

variability <- variability %>% 
  mutate(local_dsr_cv = replace_na(local_dsr_cv, 0),
         CV_S_L = replace_na(CV_S_L, 0),
         CV_C_L = replace_na(CV_C_L, 0), 
         CV_S_R = replace_na(CV_S_R, 0),
         CV_C_R = replace_na(CV_C_R, 0), 
         phi_S_L2R = replace_na(phi_S_L2R, 0),
         phi_C_L2R = replace_na(phi_C_L2R, 0),
         phi_S2C_L = replace_na(phi_S2C_L, 0),
         phi_S2C_R = replace_na(phi_S2C_R, 0),
         beta_div = replace_na(beta_div, 0))


# Diversity patterns
div_long_stable <- variability %>% 
  #filter(disturb_rate < 0.1) %>% 
  filter(condition == "stable") %>% 
  # pivot_longer(cols = c(CV_S_L, CV_C_L, CV_S_R, CV_C_R), 
  #              names_to = "variability", values_to = "CV") %>% 
  # pivot_longer(cols = c(phi_S_L2R, phi_C_L2R, phi_S2C_L, phi_S2C_R), 
  #              names_to = "synchrony", values_to = "phi") %>% 
  pivot_longer(cols = c(alpha_div, beta_div, gamma_div, beta_spatial, beta_temporal),
               names_to = "div_type", values_to = "diversity") %>% 
  mutate(#variability = factor(variability, levels = c("CV_S_L", "CV_C_L", "CV_S_R", "CV_C_R")),
         #synchrony = factor(synchrony, levels = c("phi_S_L2R", "phi_S2C_L", "phi_S2C_R", "phi_C_L2R")),
         div_type = factor(div_type, levels = c("alpha_div", "beta_div", "gamma_div", "beta_spatial", "beta_temporal")))

div_long_stable %>% 
  filter(disturb_rate < 0.1) %>% 
  ggplot(aes(x = disp_rate, y = diversity, color = as.factor(round(kernel_exp,2)))) + 
  geom_point(alpha = 0.2) + 
  #geom_smooth(method = "loess", span = .5, se = FALSE) + 
  geom_smooth(method = "gam", formula = y ~ s(x, k = 20, bs = "cs"), se = F) +
  scale_x_log10() +
  scale_color_viridis_d() +
  facet_grid(div_type ~ round(disturb_rate,2), scales = "free_y")


# variability
var_long_stable <- variability %>% 
  filter(disturb_rate < 0.1) %>% 
  filter(condition == "stable") %>% 
  pivot_longer(cols = c(CV_S_L, CV_C_L, CV_S_R, CV_C_R), 
               names_to = "variability", values_to = "CV") %>% 
  pivot_longer(cols = c(phi_S_L2R, phi_C_L2R, phi_S2C_L, phi_S2C_R), 
               names_to = "synchrony", values_to = "phi") %>% 
  mutate(variability = factor(variability, levels = c("CV_S_L", "CV_C_L", "CV_S_R", "CV_C_R")),
         synchrony = factor(synchrony, levels = c("phi_S_L2R", "phi_S2C_L", "phi_S2C_R", "phi_C_L2R")))



# Plot figures

# First, we'll just look at variability at different scales


panel_stable_var_disturbed <- var_long_stable %>% 
  filter(disturb_rate > 0) %>% 
  ggplot(aes(x = disp_rate, y = CV, color = as.factor(round(kernel_exp, 2)))) + 
  geom_point(alpha = 0.2) + 
  #geom_smooth(method = "loess", se = FALSE) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 20, bs = "cs"), se = FALSE) + 
  facet_grid(variability~round(disturb_rate,2), scales = "free_y") +
  scale_x_log10() +
  scale_color_viridis_d(option = "A") +
  labs(color = "Dispersal kernel \nexponent",
       x = "Emigration rate")

panel_stable_var_undisturbed <- var_long_stable %>% 
  filter(disturb_rate == 0) %>% 
  ggplot(aes(x = disp_rate, y = CV, color = as.factor(round(kernel_exp, 2)))) + 
  geom_point(alpha = 0.2) + 
  #geom_smooth(method = "loess", se = FALSE) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 20, bs = "cs"), se = FALSE) + 
  facet_grid(variability~round(disturb_rate,2), scales = "free_y") +
  scale_x_log10() +
  scale_color_viridis_d(option = "A") +
  labs(color = "Dispersal kernel \nexponent",
       x = "Emigration rate")

fig_stable_disturb_v_undisturb_variability <- panel_stable_var_undisturbed + panel_stable_var_disturbed +
  plot_layout(ncol = 2, guides = "collect") + 
  plot_annotation(tag_levels = "A", tag_suffix = ")")
fig_stable_disturb_v_undisturb_variability
ggsave(filename = "figures/variability_stable_disturb_vs_undisturb.pdf",
       plot = fig_stable_disturb_v_undisturb_variability,
       height = 8, width = 10)





# now look at the phi values, the scaling coefficients

panel_stable_phi_disturbed <- var_long_stable %>% 
  filter(disturb_rate > 0) %>% 
  ggplot(aes(x = disp_rate, y = phi, color = as.factor(round(kernel_exp, 2)))) + 
  geom_point(alpha = 0.2) + 
  geom_smooth(method = "gam", formula = y ~ s(x, k = 20, bs = "cs"), se = FALSE) + 
  facet_grid(synchrony~round(disturb_rate,2), scales = "free_y") +
  scale_x_log10() +
  scale_color_viridis_d(option = "A") +
  labs(color = "Dispersal kernel \nexponent",
       x = "Emigration rate", y = "Synchrony (phi)") 
panel_stable_phi_undisturbed <- var_long_stable %>% 
  filter(disturb_rate == 0) %>% 
  ggplot(aes(x = disp_rate, y = phi, color = as.factor(round(kernel_exp, 2)))) + 
  geom_point(alpha = 0.2) + 
  geom_smooth(method = "gam", formula = y ~ s(x, k = 20, bs = "cs"), se = FALSE) + 
  facet_grid(synchrony~round(disturb_rate,2), scales = "free_y") +
  scale_x_log10() +
  scale_color_viridis_d(option = "A") +
  labs(color = "Dispersal kernel \nexponent",
       x = "Emigration rate", y = "Synchrony (phi)") 
fig_stable_disturb_v_undisturb_phi <- panel_stable_phi_undisturbed + panel_stable_phi_disturbed +
  plot_layout(ncol = 2, guides = "collect") + 
  plot_annotation(tag_levels = "A", tag_suffix = ")")
fig_stable_disturb_v_undisturb_phi
ggsave(filename = "figures/variability_phi_stable_disturb_vs_undisturb.pdf",
       plot = fig_stable_disturb_v_undisturb_phi,
       height = 8, width = 10)




# now diversity-stability relationships

variability %>% 
  filter(temp_auto == 0, spat_auto == 10) %>% 
  filter(condition == "stable", disturb_rate < 0.1) %>% 
  filter(disturb_rate == 0, alpha_div > 0) %>% 
  ggplot(aes(x = alpha_div, y = CV_C_L)) + 
  geom_point(alpha = 0.3) + 
  geom_smooth(method = "lm", se=FALSE) +
  facet_grid(as.factor(kernel_exp) ~ as.factor(disp_rate), scales = "free")

variability %>% 
  filter(temp_auto == 0, spat_auto == 10) %>% 
  filter(condition == "stable", disturb_rate < 0.1) %>% 
  filter(disturb_rate == 0, alpha_div > 0) %>% 
  ggplot(aes(x = beta_div, y = phi_C_L2R)) + 
  geom_point(alpha = 0.3) + 
  geom_smooth(method = "lm", se=FALSE) +
  facet_grid(as.factor(kernel_exp) ~ as.factor(disp_rate), scales = "free")

variability %>% 
  filter(temp_auto == 0, spat_auto == 10) %>% 
  filter(condition == "stable", disturb_rate < 0.1) %>% 
  filter(disturb_rate == 0, alpha_div > 0) %>% 
  filter(kernel_exp %in% c(0,1)) %>% 
  ggplot(aes(x = alpha_div, y = CV_C_L, color = as.factor(disp_rate))) + 
  geom_point(alpha = 0.3) + 
  geom_smooth(method = "lm", se=FALSE) +
  facet_wrap(~as.factor(kernel_exp), scales = "free") + 
  scale_color_viridis_d() +
  theme_minimal() +
  theme(legend.position = "none")


# local dsrs
variability %>% 
  filter(condition == "stable") %>% 
  filter(disturb_rate == 0) %>% 
  ggplot(aes(x = local_dsr_richness, y = local_dsr_cv, color = as.factor(kernel_exp))) + 
  geom_point(alpha = .3) + 
  geom_smooth(method = "lm", se = FALSE) + 
  facet_wrap(~as.factor(disp_rate), scales = "free") + 
  scale_color_viridis_d() + 
  theme_minimal()


variability %>% 
  filter(condition == "stable") %>% 
  filter(disturb_rate == 0) %>% 
  ggplot(aes(x = local_dsr_richness, y = local_dsr_cv, color = as.factor(disp_rate))) + 
  geom_point(alpha = .3) + 
  geom_smooth(method = "lm", se = FALSE) + 
  scale_color_viridis_d() + 
  theme_minimal()
