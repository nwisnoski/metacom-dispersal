library(tidyverse)
library(here)
library(patchwork)

theme_set(theme_bw())

variability <- read_csv(here("data/variability_partitioning_2022-10-26_163959.csv"))


var_long_stable <- variability %>% 
  na.omit() %>% 
  filter(disturb_rate < 0.1) %>% 
  filter(condition == "stable") %>% 
  pivot_longer(cols = c(CV_S_L, CV_C_L, CV_S_R, CV_C_R), 
               names_to = "variability", values_to = "CV") %>% 
  pivot_longer(cols = c(phi_S_L2R, phi_C_L2R, phi_S2C_L, phi_S2C_R), 
               names_to = "synchrony", values_to = "phi") %>% 
  mutate(variability = factor(variability, levels = c("CV_S_L", "CV_C_L", "CV_S_R", "CV_C_R")),
         synchrony = factor(synchrony, levels = c("phi_S_L2R", "phi_S2C_L", "phi_S2C_R", "phi_C_L2R")))

var_long_priority <- variability %>% 
  na.omit() %>% 
  filter(disturb_rate < 0.1) %>% 
  filter(condition == "priority") %>% 
  pivot_longer(cols = c(CV_S_L, CV_C_L, CV_S_R, CV_C_R), 
               names_to = "variability", values_to = "CV") %>% 
  pivot_longer(cols = c(phi_S_L2R, phi_C_L2R, phi_S2C_L, phi_S2C_R), 
               names_to = "synchrony", values_to = "phi") %>% 
  mutate(variability = factor(variability, levels = c("CV_S_L", "CV_C_L", "CV_S_R", "CV_C_R")),
         synchrony = factor(synchrony, levels = c("phi_S_L2R", "phi_S2C_L", "phi_S2C_R", "phi_C_L2R")))


# Plot figures

# First, we'll just look at variability at different scales

# stable conditions
var_long_stable %>% 
  #filter(disturb_rate == 0) %>% 
  ggplot(aes(x = disp_rate, y = CV, color = as.factor(round(kernel_exp, 2)))) + 
  geom_point(alpha = 0.2) + 
  geom_smooth(se = FALSE) + 
  facet_grid(round(disturb_rate,2)~variability, scales = "free_y") +
  scale_x_log10() +
  scale_color_viridis_d() + 
  labs(color = "Dispersal kernel \nexponent",
       x = "Emigration rate")

panel_stable_var_disturbed <- var_long_stable %>% 
  filter(disturb_rate > 0) %>% 
  ggplot(aes(x = disp_rate, y = CV, color = as.factor(round(kernel_exp, 2)))) + 
  geom_point(alpha = 0.2) + 
  geom_smooth(se = FALSE) + 
  facet_grid(variability~round(disturb_rate,2), scales = "free_y") +
  scale_x_log10() +
  scale_color_viridis_d(option = "A") +
  labs(color = "Dispersal kernel \nexponent",
       x = "Emigration rate")

panel_stable_var_undisturbed <- var_long_stable %>% 
  filter(disturb_rate == 0) %>% 
  ggplot(aes(x = disp_rate, y = CV, color = as.factor(round(kernel_exp, 2)))) + 
  geom_point(alpha = 0.2) + 
  geom_smooth(se = FALSE) + 
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




# priority effects

panel_priority_var_disturbed <- var_long_priority %>% 
  filter(disturb_rate > 0) %>% 
  ggplot(aes(x = disp_rate, y = CV, color = as.factor(round(kernel_exp, 2)))) + 
  geom_point(alpha = 0.2) + 
  geom_smooth(se = FALSE) + 
  facet_grid(variability~round(disturb_rate,2), scales = "free_y") +
  scale_x_log10() +
  scale_color_viridis_d(option = "A") +
  labs(color = "Dispersal kernel \nexponent",
       x = "Emigration rate")

panel_priority_var_undisturbed <- var_long_priority %>% 
  filter(disturb_rate == 0) %>% 
  ggplot(aes(x = disp_rate, y = CV, color = as.factor(round(kernel_exp, 2)))) + 
  geom_point(alpha = 0.2) + 
  geom_smooth(se = FALSE) + 
  facet_grid(variability~round(disturb_rate,2), scales = "free_y") +
  scale_x_log10() +
  scale_color_viridis_d(option = "A") +
  labs(color = "Dispersal kernel \nexponent",
       x = "Emigration rate")

fig_priority_disturb_v_undisturb_variability <- panel_priority_var_undisturbed + panel_priority_var_disturbed +
  plot_layout(ncol = 2, guides = "collect") + 
  plot_annotation(tag_levels = "A", tag_suffix = ")")
fig_priority_disturb_v_undisturb_variability
ggsave(filename = "figures/variability_priority_disturb_vs_undisturb.pdf",
       plot = fig_priority_disturb_v_undisturb_variability,
       height = 8, width = 10)

var_long_priority %>% 
  ggplot(aes(x = disp_rate, y = CV, color = as.factor(round(kernel_exp, 2)))) + 
  geom_point(alpha = 0.2) + 
  geom_smooth(se = FALSE) + 
  facet_grid(variability~round(disturb_rate,2), scales = "free_y") +
  scale_x_log10() +
  labs(color = "Dispersal kernel \nexponent",
       x = "Emigration rate")


# now look at the phi values, the scaling coefficients

panel_stable_phi_disturbed <- var_long_stable %>% 
  filter(disturb_rate > 0) %>% 
  ggplot(aes(x = disp_rate, y = phi, color = as.factor(round(kernel_exp, 2)))) + 
  geom_point(alpha = 0.2) + 
  geom_smooth(se = FALSE) + 
  facet_grid(synchrony~round(disturb_rate,2), scales = "free_y") +
  scale_x_log10() +
  scale_color_viridis_d(option = "A") +
  labs(color = "Dispersal kernel \nexponent",
       x = "Emigration rate") 
panel_stable_phi_undisturbed <- var_long_stable %>% 
  filter(disturb_rate == 0) %>% 
  ggplot(aes(x = disp_rate, y = phi, color = as.factor(round(kernel_exp, 2)))) + 
  geom_point(alpha = 0.2) + 
  geom_smooth(se = FALSE) + 
  facet_grid(synchrony~round(disturb_rate,2), scales = "free_y") +
  scale_x_log10() +
  scale_color_viridis_d(option = "A") +
  labs(color = "Dispersal kernel \nexponent",
       x = "Emigration rate") 
fig_stable_disturb_v_undisturb_phi <- panel_stable_phi_undisturbed + panel_stable_phi_disturbed +
  plot_layout(ncol = 2, guides = "collect") + 
  plot_annotation(tag_levels = "A", tag_suffix = ")")
fig_stable_disturb_v_undisturb_phi
ggsave(filename = "figures/variability_phi_stable_disturb_vs_undisturb.pdf",
       plot = fig_stable_disturb_v_undisturb_phi,
       height = 8, width = 10)


panel_priority_phi_disturbed <- var_long_priority %>% 
  filter(disturb_rate > 0) %>% 
  ggplot(aes(x = disp_rate, y = phi, color = as.factor(round(kernel_exp, 2)))) + 
  geom_point(alpha = 0.2) + 
  geom_smooth(se = FALSE) + 
  facet_grid(synchrony~round(disturb_rate,2), scales = "free_y") +
  scale_x_log10() +
  scale_color_viridis_d(option = "A") +
  labs(color = "Dispersal kernel \nexponent",
       x = "Emigration rate") 
panel_priority_phi_undisturbed <- var_long_priority %>% 
  filter(disturb_rate == 0) %>% 
  ggplot(aes(x = disp_rate, y = phi, color = as.factor(round(kernel_exp, 2)))) + 
  geom_point(alpha = 0.2) + 
  geom_smooth(se = FALSE) + 
  facet_grid(synchrony~round(disturb_rate,2), scales = "free_y") +
  scale_x_log10() +
  scale_color_viridis_d(option = "A") +
  labs(color = "Dispersal kernel \nexponent",
       x = "Emigration rate") 
fig_priority_disturb_v_undisturb_phi <- panel_priority_phi_undisturbed + panel_priority_phi_disturbed +
  plot_layout(ncol = 2, guides = "collect") + 
  plot_annotation(tag_levels = "A", tag_suffix = ")")
fig_priority_disturb_v_undisturb_phi
ggsave(filename = "figures/variability_phi_priority_disturb_vs_undisturb.pdf",
       plot = fig_priority_disturb_v_undisturb_phi,
       height = 8, width = 10)






variability %>% 
  filter(condition == "stable") %>% 
  ggplot(aes(x = alpha_div, y = CV_S_L, color = as.factor(disp_rate))) + 
  geom_point(alpha = 0.3) + 
  #geom_smooth(se=FALSE) +
  facet_grid(as.factor(kernel_exp) ~ disturb_rate, scales = "free")

variability %>% 
  filter(condition == "stable") %>% 
  ggplot(aes(x = alpha_div, y = CV_C_L, color = as.factor(disp_rate))) + 
  geom_point(alpha = 0.3) + 
  #geom_smooth(se=FALSE) +
  facet_grid(as.factor(kernel_exp) ~ disturb_rate, scales = "free")

variability %>% 
  filter(condition == "stable") %>% 
  ggplot(aes(x = gamma_div, y = CV_S_R, color = as.factor(disp_rate))) + 
  geom_point(alpha = 0.3) + 
  #geom_smooth(se=FALSE) +
  facet_grid(as.factor(kernel_exp) ~ disturb_rate, scales = "free")

variability %>% 
  filter(condition == "stable") %>% 
  ggplot(aes(x = gamma_div, y = CV_C_R, color = as.factor(disp_rate))) + 
  geom_point(alpha = 0.3) + 
  #geom_smooth(se=FALSE) +
  facet_grid(as.factor(kernel_exp) ~ disturb_rate, scales = "free")

variability %>% 
  filter(condition == "stable") %>% 
  ggplot(aes(x = beta_div, y = phi_C_L2R, color = as.factor(disp_rate))) + 
  geom_point(alpha = 0.3) + 
  #geom_smooth(se=FALSE) +
  facet_grid(as.factor(kernel_exp) ~ disturb_rate, scales = "free")

variability %>% 
  filter(condition == "stable") %>% 
  filter(disturb_rate == 0) %>% 
  ggplot(aes(x = alpha_div, y = CV_C_L, color = as.factor(kernel_exp))) + 
  geom_point(alpha = 0.3) + 
  geom_smooth(method = "lm", se=FALSE) +
  facet_wrap(~as.factor(disp_rate), scales = "free") + 
  scale_color_viridis_d() +
  theme_minimal()

variability %>% 
  filter(condition == "stable") %>% 
  filter(disturb_rate == 0) %>% 
  ggplot(aes(x = beta_div, y = phi_C_L2R, color = as.factor(kernel_exp))) + 
  geom_point(alpha = 0.3) + 
  geom_smooth(method = "lm", se=FALSE) +
  facet_wrap(~as.factor(disp_rate), scales = "free") + 
  scale_color_viridis_d() +
  theme_minimal()

variability %>% 
  filter(condition == "stable") %>% 
  filter(disturb_rate == 0) %>% 
  ggplot(aes(x = gamma_div, y = CV_C_R, color = as.factor(kernel_exp))) + 
  geom_point(alpha = 0.3) + 
  geom_smooth(method = "lm", se=FALSE) +
  facet_wrap(~as.factor(disp_rate), scales = "free") + 
  scale_color_viridis_d() +
  theme_minimal()



# variability %>% 
#   filter(condition == "stable") %>% 
#   filter(disturb_rate == 0) %>% 
#   ggplot(aes(x = alpha_div, y = CV_C_L, color = as.factor(disp_rate))) + 
#   geom_point(alpha = 0.3) + 
#   geom_smooth(method = "lm", se=FALSE) +
#   facet_wrap(~as.factor(round(kernel_exp,2)), scales = "free") + 
#   scale_color_viridis_d() +
#   theme_minimal() +
#   labs(color = "Emigration rate", 
#        x = "Mean alpha-diversity", 
#        y = "Local community variability")



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
