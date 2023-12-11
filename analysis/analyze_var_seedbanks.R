library(tidyverse)
library(here)
library(patchwork)

theme_set(theme_bw())

variability <- read_csv(here("sim_output/variability_partitioning_disp_kernel_seedbank_2023-12-08_140313.336936.csv")) |> 
  bind_rows(read_csv(here("sim_output/variability_partitioning_disp_kernel_seedbank_2023-12-08_153919.129852.csv"))) |> 
  bind_rows(read_csv(here("sim_output/variability_partitioning_disp_kernel_seedbank_2023-12-08_173926.212731.csv"))) |> 
  bind_rows(read_csv(here("sim_output/variability_partitioning_disp_kernel_seedbank_2023-12-11_140319.539096.csv")))
  
variability <- variability |> mutate(beta_div = (gamma_div - alpha_div))
  

kernel_exps <- sort(unique(variability$kernel_exp))
disp_rates <- sort(unique(variability$disp_rate))
disturb_rates <- sort(unique(variability$disturb_rate))

#### Stable


# Diversity patterns
div_long_stable <- variability %>% 
  filter(disturb_rate == 0) %>% 
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

diversity_plot_no_sb <- div_long_stable %>% 
  filter(germ_rate == 1) %>% 
  ggplot(aes(x = disp_rate, y = diversity, color = as.factor(round(kernel_exp,2)))) + 
  geom_point(alpha = 0.2) + 
  #geom_smooth(method = "loess", span = .5, se = FALSE) + 
  geom_smooth(method = "gam", formula = y ~ s(x, k = 20, bs = "cs"), se = F) +
  scale_x_log10() +
  scale_color_viridis_d() +
  facet_grid(div_type ~ round(germ_rate,2), scales = "free_y") +
  labs(x = "Emigration rate", y = "Diversity", color = "Kernel exponent")

diversity_plot_with_sb <- div_long_stable %>%
  filter(germ_rate < 1) %>%
  ggplot(aes(x = disp_rate, y = diversity, color = as.factor(round(kernel_exp,2)))) +
  geom_point(alpha = 0.2) +
  #geom_smooth(method = "loess", span = .5, se = FALSE) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 20, bs = "cs"), se = F) +
  scale_x_log10() +
  scale_color_viridis_d() +
  facet_grid(div_type ~ round(germ_rate,2), scales = "free_y") +
  labs(x = "Emigration rate", y = "Diversity", color = "Kernel exponent")

figure_diversity <- diversity_plot_no_sb + diversity_plot_with_sb +
  plot_layout(ncol = 2, guides = "collect") + 
  plot_annotation(tag_levels = "A", tag_suffix = ")")
figure_diversity
# ggsave(plot = figure_diversity, filename = "figures/fig_diversity.pdf", width = 10, height = 8)


# Diversity patterns with disturbance
div_long_stable <- variability %>% 
  filter(disturb_rate == 0.01) %>% 
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

diversity_plot_no_sb <- div_long_stable %>% 
  filter(germ_rate == 1) %>% 
  ggplot(aes(x = disp_rate, y = diversity, color = as.factor(round(kernel_exp,2)))) + 
  geom_point(alpha = 0.2) + 
  #geom_smooth(method = "loess", span = .5, se = FALSE) + 
  geom_smooth(method = "gam", formula = y ~ s(x, k = 20, bs = "cs"), se = F) +
  scale_x_log10() +
  scale_color_viridis_d() +
  facet_grid(div_type ~ round(germ_rate,2), scales = "free_y") +
  labs(x = "Emigration rate", y = "Diversity", color = "Kernel exponent")

diversity_plot_with_sb <- div_long_stable %>%
  filter(germ_rate < 1) %>%
  ggplot(aes(x = disp_rate, y = diversity, color = as.factor(round(kernel_exp,2)))) +
  geom_point(alpha = 0.2) +
  #geom_smooth(method = "loess", span = .5, se = FALSE) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 20, bs = "cs"), se = F) +
  scale_x_log10() +
  scale_color_viridis_d() +
  facet_grid(div_type ~ round(germ_rate,2), scales = "free_y") +
  labs(x = "Emigration rate", y = "Diversity", color = "Kernel exponent")

figure_diversity <- diversity_plot_no_sb + diversity_plot_with_sb +
  plot_layout(ncol = 2, guides = "collect") + 
  plot_annotation(tag_levels = "A", tag_suffix = ")")
figure_diversity
# ggsave(plot = figure_diversity, filename = "figures/fig_diversity.pdf", width = 10, height = 8)



#### 


# variability no disturbance
var_long_stable <- variability %>% 
  filter(disturb_rate == 0) %>% 
  filter(condition == "stable") %>% 
  pivot_longer(cols = c(CV_S_L, CV_C_L, CV_S_R, CV_C_R), 
               names_to = "variability", values_to = "CV") %>% 
  pivot_longer(cols = c(phi_S_L2R, phi_C_L2R, phi_S2C_L, phi_S2C_R), 
               names_to = "synchrony", values_to = "phi") %>% 
  mutate(variability = factor(variability, levels = c("CV_S_L", "CV_C_L", "CV_S_R", "CV_C_R")),
         synchrony = factor(synchrony, levels = c("phi_S_L2R", "phi_S2C_L", "phi_S2C_R", "phi_C_L2R")))

# Plot figures

# First, we'll just look at variability at different scales

panel_stable_var_no_sb <- var_long_stable %>% 
  filter(germ_rate == 1, gamma_div > 0) %>% 
  ggplot(aes(x = disp_rate, y = CV, color = as.factor(round(kernel_exp, 2)))) + 
  geom_point(alpha = 0.2) + 
  #geom_smooth(method = "loess", se = FALSE) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5, bs = "cs"), se = FALSE) + 
  facet_grid(variability~round(germ_rate,2), scales = "free_y") +
  scale_x_log10() +
  scale_color_viridis_d(option = "A") +
  labs(color = "Dispersal kernel \nexponent",
       x = "Emigration rate")

panel_stable_var_with_sb <- var_long_stable %>% 
  filter(germ_rate < 1, gamma_div > 0) %>% 
  ggplot(aes(x = disp_rate, y = CV, color = as.factor(round(kernel_exp, 2)))) + 
  geom_point(alpha = 0.2) + 
  #geom_smooth(method = "loess", se = FALSE) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5, bs = "cs"), se = FALSE) + 
  facet_grid(variability~round(germ_rate,2), scales = "free_y") +
  scale_x_log10() +
  scale_color_viridis_d(option = "A") +
  labs(color = "Dispersal kernel \nexponent",
       x = "Emigration rate")


fig_stable_variability <- panel_stable_var_no_sb + panel_stable_var_with_sb +
  plot_layout(ncol = 2, guides = "collect") + 
  plot_annotation(tag_levels = "A", tag_suffix = ")")
fig_stable_variability

# ggsave(filename = "figures/variability_stable_disturb_vs_undisturb.pdf",
#        plot = fig_stable_disturb_v_undisturb_variability,
#        height = 8, width = 10)

var_long_stable %>% 
  filter(gamma_div > 0) %>% 
  ggplot(aes(x = disp_rate, y = CV, color = as.factor(round(kernel_exp, 2)))) + 
  geom_point(alpha = 0.2) + 
  #geom_smooth(method = "loess", se = FALSE) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5, bs = "cs"), se = FALSE) + 
  facet_grid(variability~round(germ_rate,2), scales = "free_y") +
  scale_x_log10() +
  scale_color_viridis_d(option = "A") +
  labs(color = "Dispersal kernel \nexponent",
       x = "Emigration rate")



# variability with disturbance
var_long_stable <- variability %>% 
  filter(disturb_rate == .01) %>% 
  filter(condition == "stable") %>% 
  pivot_longer(cols = c(CV_S_L, CV_C_L, CV_S_R, CV_C_R), 
               names_to = "variability", values_to = "CV") %>% 
  pivot_longer(cols = c(phi_S_L2R, phi_C_L2R, phi_S2C_L, phi_S2C_R), 
               names_to = "synchrony", values_to = "phi") %>% 
  mutate(variability = factor(variability, levels = c("CV_S_L", "CV_C_L", "CV_S_R", "CV_C_R")),
         synchrony = factor(synchrony, levels = c("phi_S_L2R", "phi_S2C_L", "phi_S2C_R", "phi_C_L2R")))

# Plot figures

# First, we'll just look at variability at different scales

panel_stable_var_no_sb <- var_long_stable %>% 
  filter(germ_rate == 1, gamma_div > 0) %>% 
  ggplot(aes(x = disp_rate, y = CV, color = as.factor(round(kernel_exp, 2)))) + 
  geom_point(alpha = 0.2) + 
  #geom_smooth(method = "loess", se = FALSE) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5, bs = "cs"), se = FALSE) + 
  facet_grid(variability~round(germ_rate,2), scales = "free_y") +
  scale_x_log10() +
  #scale_y_continuous(limits = c(0, 2)) +
  scale_color_viridis_d(option = "A") +
  labs(color = "Dispersal kernel \nexponent",
       x = "Emigration rate")

panel_stable_var_with_sb <- var_long_stable %>% 
  filter(germ_rate < 1, gamma_div > 0) %>% 
  ggplot(aes(x = disp_rate, y = CV, color = as.factor(round(kernel_exp, 2)))) + 
  geom_point(alpha = 0.2) + 
  #geom_smooth(method = "loess", se = FALSE) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5, bs = "cs"), se = FALSE) + 
  facet_grid(variability~round(germ_rate,2), scales = "free_y") +
  scale_x_log10() + 
  #scale_y_continuous(limits = c(0,2)) +
  scale_color_viridis_d(option = "A") +
  labs(color = "Dispersal kernel \nexponent",
       x = "Emigration rate")


fig_stable_variability <- panel_stable_var_no_sb + panel_stable_var_with_sb +
  plot_layout(ncol = 2, guides = "collect") + 
  plot_annotation(tag_levels = "A", tag_suffix = ")")
fig_stable_variability

var_long_stable %>% 
  filter(gamma_div > 0) %>% 
  ggplot(aes(x = disp_rate, y = CV, color = as.factor(round(kernel_exp, 2)))) + 
  geom_point(alpha = 0.2) + 
  #geom_smooth(method = "loess", se = FALSE) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5, bs = "cs"), se = FALSE) + 
  facet_grid(variability~round(germ_rate,2), scales = "free_y") +
  scale_x_log10() +
  scale_color_viridis_d(option = "A") +
  labs(color = "Dispersal kernel \nexponent",
       x = "Emigration rate")


# now look at the phi values, the scaling coefficients
var_long_stable <- variability %>% 
  #filter(disturb_rate == .01) %>% 
  filter(condition == "stable") %>% 
  pivot_longer(cols = c(CV_S_L, CV_C_L, CV_S_R, CV_C_R), 
               names_to = "variability", values_to = "CV") %>% 
  pivot_longer(cols = c(phi_S_L2R, phi_C_L2R, phi_S2C_L, phi_S2C_R), 
               names_to = "synchrony", values_to = "phi") %>% 
  mutate(variability = factor(variability, levels = c("CV_S_L", "CV_C_L", "CV_S_R", "CV_C_R")),
         synchrony = factor(synchrony, levels = c("phi_S_L2R", "phi_S2C_L", "phi_S2C_R", "phi_C_L2R")))

# no disturbance
panel_stable_phi_disturbed <- var_long_stable %>% 
  filter(disturb_rate > 0, gamma_div > 0) %>% 
  ggplot(aes(x = disp_rate, y = phi, color = as.factor(round(kernel_exp, 2)))) + 
  geom_point(alpha = 0.2) + 
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5, bs = "cs"), se = FALSE) + 
  facet_grid(synchrony~round(germ_rate,2), scales = "free_y") +
  scale_x_log10() +
  scale_color_viridis_d(option = "A") +
  labs(color = "Dispersal kernel \nexponent",
       x = "Emigration rate", y = "Synchrony (phi)") 
panel_stable_phi_disturbed

panel_stable_phi_undisturbed <- var_long_stable %>% 
  filter(disturb_rate == 0, gamma_div > 0) %>% 
  ggplot(aes(x = disp_rate, y = phi, color = as.factor(round(kernel_exp, 2)))) + 
  geom_point(alpha = 0.2) + 
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5, bs = "cs"), se = FALSE) + 
  facet_grid(synchrony~round(germ_rate,2), scales = "free_y") +
  scale_x_log10() +
  scale_color_viridis_d(option = "A") +
  labs(color = "Dispersal kernel \nexponent",
       x = "Emigration rate", y = "Synchrony (phi)") 
panel_stable_phi_undisturbed

# fig_stable_disturb_v_undisturb_phi <- panel_stable_phi_undisturbed + panel_stable_phi_disturbed +
#   plot_layout(ncol = 2, guides = "collect") + 
#   plot_annotation(tag_levels = "A", tag_suffix = ")")
# fig_stable_disturb_v_undisturb_phi
# ggsave(filename = "figures/variability_phi_stable_disturb_vs_undisturb.pdf",
#        plot = fig_stable_disturb_v_undisturb_phi,
#        height = 8, width = 10)
# 
# var_long_stable %>% 
#   filter(gamma_div > 0) %>% 
#   ggplot(aes(x = disp_rate, y = phi, color = as.factor(round(kernel_exp, 2)))) + 
#   geom_point(alpha = 0.2) + 
#   geom_smooth(method = "gam", formula = y ~ s(x, k = 20, bs = "cs"), se = FALSE) + 
#   facet_grid(synchrony~round(germ_rate,2), scales = "free_y") +
#   scale_x_log10() +
#   scale_color_viridis_d(option = "A") +
#   labs(color = "Dispersal kernel \nexponent",
#        x = "Emigration rate", y = "Synchrony (phi)") 
# 

# now diversity-stability relationships


variability %>% 
  filter(condition == "stable", disturb_rate ==0) %>% 
  filter(germ_rate == 1, alpha_div > 0) %>% 
  ggplot(aes(x = beta_div, y = phi_C_L2R)) + 
  geom_point(alpha = 0.3) + 
  geom_smooth(se=FALSE) 

variability %>% 
  filter(condition == "stable", disturb_rate ==0) %>% 
  filter(germ_rate < 1, alpha_div > 0) %>% 
  ggplot(aes(x = beta_div, y = phi_C_L2R)) + 
  geom_point(alpha = 0.3) + 
  geom_smooth(se=FALSE) 

variability %>% 
  filter(condition == "stable", disturb_rate ==.01) %>% 
  filter(germ_rate == 1, alpha_div > 0) %>% 
  ggplot(aes(x = beta_div, y = phi_C_L2R)) + 
  geom_point(alpha = 0.3) + 
  geom_smooth(se=FALSE) 

variability %>% 
  filter(condition == "stable", disturb_rate ==0.01) %>% 
  filter(germ_rate < 1, alpha_div > 0) %>% 
  ggplot(aes(x = beta_div, y = phi_C_L2R)) + 
  geom_point(alpha = 0.3) + 
  geom_smooth(se=FALSE) 


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
  ggplot(aes(x = local_dsr_richness, y = local_dsr_cv, color = as.factor(kernel_exp))) + 
  geom_point(alpha = .3) + 
  geom_smooth(method = "lm", se = FALSE) + 
  facet_grid(as.factor(round(germ_rate,2))~as.factor(disp_rate), scales = "free") + 
  scale_color_viridis_d() + 
  theme_minimal()


variability %>% 
  filter(condition == "stable") %>% 
  filter(disturb_rate == 0) %>% 
  ggplot(aes(x = local_dsr_richness, y = local_dsr_cv, color = as.factor(kernel_exp))) + 
  geom_point(alpha = .3) + 
  geom_smooth(method = "lm", se = FALSE) + 
  scale_color_viridis_d() + 
  theme_minimal() +
  facet_grid(~round(germ_rate,2), scales = "free_y")









###################### Equal

#### Stable


# Diversity patterns
div_long_equal <- variability %>% 
  filter(disturb_rate == 0) %>% 
  filter(condition == "equal") %>% 
  # pivot_longer(cols = c(CV_S_L, CV_C_L, CV_S_R, CV_C_R), 
  #              names_to = "variability", values_to = "CV") %>% 
  # pivot_longer(cols = c(phi_S_L2R, phi_C_L2R, phi_S2C_L, phi_S2C_R), 
  #              names_to = "synchrony", values_to = "phi") %>% 
  pivot_longer(cols = c(alpha_div, beta_div, gamma_div, beta_spatial, beta_temporal),
               names_to = "div_type", values_to = "diversity") %>% 
  mutate(#variability = factor(variability, levels = c("CV_S_L", "CV_C_L", "CV_S_R", "CV_C_R")),
    #synchrony = factor(synchrony, levels = c("phi_S_L2R", "phi_S2C_L", "phi_S2C_R", "phi_C_L2R")),
    div_type = factor(div_type, levels = c("alpha_div", "beta_div", "gamma_div", "beta_spatial", "beta_temporal")))

diversity_plot_no_sb <- div_long_equal %>% 
  filter(germ_rate == 1) %>% 
  ggplot(aes(x = disp_rate, y = diversity, color = as.factor(round(kernel_exp,2)))) + 
  geom_point(alpha = 0.2) + 
  #geom_smooth(method = "loess", span = .5, se = FALSE) + 
  geom_smooth(method = "gam", formula = y ~ s(x, k = 20, bs = "cs"), se = F) +
  scale_x_log10() +
  scale_color_viridis_d() +
  facet_grid(div_type ~ round(germ_rate,2), scales = "free_y") +
  labs(x = "Emigration rate", y = "Diversity", color = "Kernel exponent")

diversity_plot_with_sb <- div_long_equal %>%
  filter(germ_rate < 1) %>%
  ggplot(aes(x = disp_rate, y = diversity, color = as.factor(round(kernel_exp,2)))) +
  geom_point(alpha = 0.2) +
  #geom_smooth(method = "loess", span = .5, se = FALSE) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 20, bs = "cs"), se = F) +
  scale_x_log10() +
  scale_color_viridis_d() +
  facet_grid(div_type ~ round(germ_rate,2), scales = "free_y") +
  labs(x = "Emigration rate", y = "Diversity", color = "Kernel exponent")

figure_diversity <- diversity_plot_no_sb + diversity_plot_with_sb +
  plot_layout(ncol = 2, guides = "collect") + 
  plot_annotation(tag_levels = "A", tag_suffix = ")")
figure_diversity
# ggsave(plot = figure_diversity, filename = "figures/fig_diversity.pdf", width = 10, height = 8)


# Diversity patterns with disturbance
div_long_equal <- variability %>% 
  filter(disturb_rate == 0.01) %>% 
  filter(condition == "equal") %>% 
  # pivot_longer(cols = c(CV_S_L, CV_C_L, CV_S_R, CV_C_R), 
  #              names_to = "variability", values_to = "CV") %>% 
  # pivot_longer(cols = c(phi_S_L2R, phi_C_L2R, phi_S2C_L, phi_S2C_R), 
  #              names_to = "synchrony", values_to = "phi") %>% 
  pivot_longer(cols = c(alpha_div, beta_div, gamma_div, beta_spatial, beta_temporal),
               names_to = "div_type", values_to = "diversity") %>% 
  mutate(#variability = factor(variability, levels = c("CV_S_L", "CV_C_L", "CV_S_R", "CV_C_R")),
    #synchrony = factor(synchrony, levels = c("phi_S_L2R", "phi_S2C_L", "phi_S2C_R", "phi_C_L2R")),
    div_type = factor(div_type, levels = c("alpha_div", "beta_div", "gamma_div", "beta_spatial", "beta_temporal")))

diversity_plot_no_sb <- div_long_equal %>% 
  filter(germ_rate == 1) %>% 
  ggplot(aes(x = disp_rate, y = diversity, color = as.factor(round(kernel_exp,2)))) + 
  geom_point(alpha = 0.2) + 
  #geom_smooth(method = "loess", span = .5, se = FALSE) + 
  geom_smooth(method = "gam", formula = y ~ s(x, k = 20, bs = "cs"), se = F) +
  scale_x_log10() +
  scale_color_viridis_d() +
  facet_grid(div_type ~ round(germ_rate,2), scales = "free_y") +
  labs(x = "Emigration rate", y = "Diversity", color = "Kernel exponent")

diversity_plot_with_sb <- div_long_equal %>%
  filter(germ_rate < 1) %>%
  ggplot(aes(x = disp_rate, y = diversity, color = as.factor(round(kernel_exp,2)))) +
  geom_point(alpha = 0.2) +
  #geom_smooth(method = "loess", span = .5, se = FALSE) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 20, bs = "cs"), se = F) +
  scale_x_log10() +
  scale_color_viridis_d() +
  facet_grid(div_type ~ round(germ_rate,2), scales = "free_y") +
  labs(x = "Emigration rate", y = "Diversity", color = "Kernel exponent")

figure_diversity <- diversity_plot_no_sb + diversity_plot_with_sb +
  plot_layout(ncol = 2, guides = "collect") + 
  plot_annotation(tag_levels = "A", tag_suffix = ")")
figure_diversity
# ggsave(plot = figure_diversity, filename = "figures/fig_diversity.pdf", width = 10, height = 8)



#### 


# variability no disturbance
var_long_equal <- variability %>% 
  filter(disturb_rate == 0) %>% 
  filter(condition == "equal") %>% 
  pivot_longer(cols = c(CV_S_L, CV_C_L, CV_S_R, CV_C_R), 
               names_to = "variability", values_to = "CV") %>% 
  pivot_longer(cols = c(phi_S_L2R, phi_C_L2R, phi_S2C_L, phi_S2C_R), 
               names_to = "synchrony", values_to = "phi") %>% 
  mutate(variability = factor(variability, levels = c("CV_S_L", "CV_C_L", "CV_S_R", "CV_C_R")),
         synchrony = factor(synchrony, levels = c("phi_S_L2R", "phi_S2C_L", "phi_S2C_R", "phi_C_L2R")))

# Plot figures

# First, we'll just look at variability at different scales

panel_equal_var_no_sb <- var_long_equal %>% 
  filter(germ_rate == 1, gamma_div > 0) %>% 
  ggplot(aes(x = disp_rate, y = CV, color = as.factor(round(kernel_exp, 2)))) + 
  geom_point(alpha = 0.2) + 
  #geom_smooth(method = "loess", se = FALSE) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5, bs = "cs"), se = FALSE) + 
  facet_grid(variability~round(germ_rate,2), scales = "free_y") +
  scale_x_log10() +
  scale_color_viridis_d(option = "A") +
  labs(color = "Dispersal kernel \nexponent",
       x = "Emigration rate")

panel_equal_var_with_sb <- var_long_equal %>% 
  filter(germ_rate < 1, gamma_div > 0) %>% 
  ggplot(aes(x = disp_rate, y = CV, color = as.factor(round(kernel_exp, 2)))) + 
  geom_point(alpha = 0.2) + 
  #geom_smooth(method = "loess", se = FALSE) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 20, bs = "cs"), se = FALSE) + 
  facet_grid(variability~round(germ_rate,2), scales = "free_y") +
  scale_x_log10() +
  scale_color_viridis_d(option = "A") +
  labs(color = "Dispersal kernel \nexponent",
       x = "Emigration rate")


fig_equal_variability <- panel_equal_var_no_sb + panel_equal_var_with_sb +
  plot_layout(ncol = 2, guides = "collect") + 
  plot_annotation(tag_levels = "A", tag_suffix = ")")
fig_equal_variability

ggsave(filename = "figures/variability_equal_disturb_vs_undisturb.pdf",
       plot = fig_equal_disturb_v_undisturb_variability,
       height = 8, width = 10)

var_long_equal %>% 
  filter(gamma_div > 0) %>% 
  ggplot(aes(x = disp_rate, y = CV, color = as.factor(round(kernel_exp, 2)))) + 
  geom_point(alpha = 0.2) + 
  #geom_smooth(method = "loess", se = FALSE) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 20, bs = "cs"), se = FALSE) + 
  facet_grid(variability~round(germ_rate,2), scales = "free_y") +
  scale_x_log10() +
  scale_color_viridis_d(option = "A") +
  labs(color = "Dispersal kernel \nexponent",
       x = "Emigration rate")



# variability with disturbance
var_long_equal <- variability %>% 
  filter(disturb_rate == .01) %>% 
  filter(condition == "equal") %>% 
  pivot_longer(cols = c(CV_S_L, CV_C_L, CV_S_R, CV_C_R), 
               names_to = "variability", values_to = "CV") %>% 
  pivot_longer(cols = c(phi_S_L2R, phi_C_L2R, phi_S2C_L, phi_S2C_R), 
               names_to = "synchrony", values_to = "phi") %>% 
  mutate(variability = factor(variability, levels = c("CV_S_L", "CV_C_L", "CV_S_R", "CV_C_R")),
         synchrony = factor(synchrony, levels = c("phi_S_L2R", "phi_S2C_L", "phi_S2C_R", "phi_C_L2R")))

# Plot figures

# First, we'll just look at variability at different scales

panel_equal_var_no_sb <- var_long_equal %>% 
  filter(germ_rate == 1, gamma_div > 0, disp_rate > 1e-4) %>% 
  ggplot(aes(x = disp_rate, y = CV, color = as.factor(round(kernel_exp, 2)))) + 
  geom_point(alpha = 0.2) + 
  #geom_smooth(method = "loess", se = FALSE) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5, bs = "cs"), se = FALSE) + 
  facet_grid(variability~round(germ_rate,2), scales = "free_y") +
  scale_x_log10() +
  #scale_y_continuous(limits = c(0, 2)) +
  scale_color_viridis_d(option = "A") +
  labs(color = "Dispersal kernel \nexponent",
       x = "Emigration rate")

panel_equal_var_with_sb <- var_long_equal %>% 
  filter(germ_rate < 1, gamma_div > 0, disp_rate > 1e-4) %>% 
  ggplot(aes(x = disp_rate, y = CV, color = as.factor(round(kernel_exp, 2)))) + 
  geom_point(alpha = 0.2) + 
  #geom_smooth(method = "loess", se = FALSE) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 20, bs = "cs"), se = FALSE) + 
  facet_grid(variability~round(germ_rate,2), scales = "free_y") +
  scale_x_log10() + 
  #scale_y_continuous(limits = c(0,2)) +
  scale_color_viridis_d(option = "A") +
  labs(color = "Dispersal kernel \nexponent",
       x = "Emigration rate")


fig_equal_variability <- panel_equal_var_no_sb + panel_equal_var_with_sb +
  plot_layout(ncol = 2, guides = "collect") + 
  plot_annotation(tag_levels = "A", tag_suffix = ")")
fig_equal_variability

var_long_equal %>% 
  filter(gamma_div > 0) %>% 
  ggplot(aes(x = disp_rate, y = CV, color = as.factor(round(kernel_exp, 2)))) + 
  geom_point(alpha = 0.2) + 
  #geom_smooth(method = "loess", se = FALSE) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 20, bs = "cs"), se = FALSE) + 
  facet_grid(variability~round(germ_rate,2), scales = "free_y") +
  scale_x_log10() +
  scale_color_viridis_d(option = "A") +
  labs(color = "Dispersal kernel \nexponent",
       x = "Emigration rate")


# now look at the phi values, the scaling coefficients
var_long_equal <- variability %>% 
  #filter(disturb_rate == .01) %>% 
  filter(condition == "equal") %>% 
  pivot_longer(cols = c(CV_S_L, CV_C_L, CV_S_R, CV_C_R), 
               names_to = "variability", values_to = "CV") %>% 
  pivot_longer(cols = c(phi_S_L2R, phi_C_L2R, phi_S2C_L, phi_S2C_R), 
               names_to = "synchrony", values_to = "phi") %>% 
  mutate(variability = factor(variability, levels = c("CV_S_L", "CV_C_L", "CV_S_R", "CV_C_R")),
         synchrony = factor(synchrony, levels = c("phi_S_L2R", "phi_S2C_L", "phi_S2C_R", "phi_C_L2R")))

# no disturbance
panel_equal_phi_disturbed <- var_long_equal %>% 
  filter(disturb_rate > 0, gamma_div > 0) %>% 
  ggplot(aes(x = disp_rate, y = phi, color = as.factor(round(kernel_exp, 2)))) + 
  geom_point(alpha = 0.2) + 
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5, bs = "cs"), se = FALSE) + 
  facet_grid(synchrony~round(germ_rate,2), scales = "free_y") +
  scale_x_log10() +
  scale_color_viridis_d(option = "A") +
  labs(color = "Dispersal kernel \nexponent",
       x = "Emigration rate", y = "Synchrony (phi)") 
panel_equal_phi_disturbed

panel_equal_phi_undisturbed <- var_long_equal %>% 
  filter(disturb_rate == 0, gamma_div > 0) %>% 
  ggplot(aes(x = disp_rate, y = phi, color = as.factor(round(kernel_exp, 2)))) + 
  geom_point(alpha = 0.2) + 
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5, bs = "cs"), se = FALSE) + 
  facet_grid(synchrony~round(germ_rate,2), scales = "free_y") +
  scale_x_log10() +
  scale_color_viridis_d(option = "A") +
  labs(color = "Dispersal kernel \nexponent",
       x = "Emigration rate", y = "Synchrony (phi)") 
panel_equal_phi_undisturbed

# fig_equal_disturb_v_undisturb_phi <- panel_equal_phi_undisturbed + panel_equal_phi_disturbed +
#   plot_layout(ncol = 2, guides = "collect") + 
#   plot_annotation(tag_levels = "A", tag_suffix = ")")
# fig_equal_disturb_v_undisturb_phi
# ggsave(filename = "figures/variability_phi_equal_disturb_vs_undisturb.pdf",
#        plot = fig_equal_disturb_v_undisturb_phi,
#        height = 8, width = 10)
# 
# var_long_equal %>% 
#   filter(gamma_div > 0) %>% 
#   ggplot(aes(x = disp_rate, y = phi, color = as.factor(round(kernel_exp, 2)))) + 
#   geom_point(alpha = 0.2) + 
#   geom_smooth(method = "gam", formula = y ~ s(x, k = 20, bs = "cs"), se = FALSE) + 
#   facet_grid(synchrony~round(germ_rate,2), scales = "free_y") +
#   scale_x_log10() +
#   scale_color_viridis_d(option = "A") +
#   labs(color = "Dispersal kernel \nexponent",
#        x = "Emigration rate", y = "Synchrony (phi)") 
# 

# now diversity-stability relationships


variability %>% 
  filter(condition == "equal", disturb_rate ==0) %>% 
  filter(germ_rate == 1, alpha_div > 0) %>% 
  ggplot(aes(x = beta_div, y = phi_C_L2R)) + 
  geom_point(alpha = 0.3) + 
  geom_smooth(se=FALSE) 

variability %>% 
  filter(condition == "equal", disturb_rate ==0) %>% 
  filter(germ_rate < 1, alpha_div > 0) %>% 
  ggplot(aes(x = beta_div, y = phi_C_L2R)) + 
  geom_point(alpha = 0.3) + 
  geom_smooth(se=FALSE) 

variability %>% 
  filter(condition == "equal", disturb_rate ==.01) %>% 
  filter(germ_rate == 1, alpha_div > 0) %>% 
  ggplot(aes(x = beta_div, y = phi_C_L2R)) + 
  geom_point(alpha = 0.3) + 
  geom_smooth(se=FALSE) 

variability %>% 
  filter(condition == "equal", disturb_rate ==0.01) %>% 
  filter(germ_rate < 1, alpha_div > 0) %>% 
  ggplot(aes(x = beta_div, y = phi_C_L2R)) + 
  geom_point(alpha = 0.3) + 
  geom_smooth(se=FALSE) 


variability %>% 
  filter(temp_auto == 0, spat_auto == 10) %>% 
  filter(condition == "equal", disturb_rate < 0.1) %>% 
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
  filter(condition == "equal") %>% 
  ggplot(aes(x = local_dsr_richness, y = local_dsr_cv, color = as.factor(kernel_exp))) + 
  geom_point(alpha = .3) + 
  geom_smooth(method = "lm", se = FALSE) + 
  facet_grid(as.factor(round(germ_rate,2))~as.factor(disp_rate), scales = "free") + 
  scale_color_viridis_d() + 
  theme_minimal()


variability %>% 
  filter(condition == "equal") %>% 
  filter(disturb_rate == 0) %>% 
  ggplot(aes(x = local_dsr_richness, y = local_dsr_cv, color = as.factor(kernel_exp))) + 
  geom_point(alpha = .3) + 
  geom_smooth(method = "lm", se = FALSE) + 
  scale_color_viridis_d() + 
  theme_minimal() +
  facet_grid(~round(germ_rate,2), scales = "free_y")