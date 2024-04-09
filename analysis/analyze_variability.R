library(tidyverse)
library(here)
library(patchwork)
library(data.table)
library(lme4)
library(lmerTest)

theme_set(theme_bw())

subfolder <- "sim_output/2024-02-09_sims/"
comp_scenario <- "equal"
# comp_scenario <- "stable"
# comp_scenario <- "priority"

# load and combine reps
variability <- data.table()
i = 1
for(file in list.files(path = subfolder, pattern = "summary\\.csv$")){
  this_run <- fread(paste(subfolder, file, sep = "/"))
  this_run$rep <- i
  variability <- bind_rows(variability, this_run)
  i = i + 1
} 

variability <- variability |> 
  rename(env_temp_cv_mean = env_temp_cv_mean.env_temp_cv_mean, 
         env_spat_cv_mean = env_spat_cv_mean.env_spat_cv_mean) |> 
  mutate(local_dsr_cv = replace_na(local_dsr_cv, 0),
         CV_S_L = replace_na(CV_S_L, 0),
         CV_C_L = replace_na(CV_C_L, 0), 
         CV_S_R = replace_na(CV_S_R, 0),
         CV_C_R = replace_na(CV_C_R, 0), 
         phi_S_L2R = replace_na(phi_S_L2R, 0),
         phi_C_L2R = replace_na(phi_C_L2R, 0),
         phi_S2C_L = replace_na(phi_S2C_L, 0),
         phi_S2C_R = replace_na(phi_S2C_R, 0),
         beta_div = replace_na(beta_div, 0),
         env_spat_cv_mean = replace_na(env_spat_cv_mean, 0),
         env_temp_cv_mean = replace_na(env_temp_cv_mean, 0)) |> 
  filter(spat_heterogeneity != 1) 
  

kernel_exps <- sort(unique(variability$kernel_exp))
disp_rates <- sort(unique(variability$disp_rate))
disturb_rates <- sort(unique(variability$disturb_rate))

# Diversity patterns
div_long <- variability |> 
  filter(condition == comp_scenario) |>  
  # pivot_longer(cols = c(CV_S_L, CV_C_L, CV_S_R, CV_C_R), 
  #              names_to = "variability", values_to = "CV")  |>  
  # pivot_longer(cols = c(phi_S_L2R, phi_C_L2R, phi_S2C_L, phi_S2C_R), 
  #              names_to = "synchrony", values_to = "phi") |> 
  pivot_longer(cols = c(alpha_div, gamma_div, beta_spatial, beta_temporal),
               names_to = "div_type", values_to = "diversity") |> 
  mutate(#variability = factor(variability, levels = c("CV_S_L", "CV_C_L", "CV_S_R", "CV_C_R")),
         #synchrony = factor(synchrony, levels = c("phi_S_L2R", "phi_S2C_L", "phi_S2C_R", "phi_C_L2R")),
         div_type = factor(div_type, 
                           levels = c("alpha_div", "gamma_div", "beta_spatial", "beta_temporal"),
                           labels = c("alpha diversity", "gamma diversity", str_wrap("spatial beta diversity", width = 15), str_wrap("temporal beta diversity", width = 15)))) |> 
  mutate(kernel_exp = as.factor(kernel_exp),
         spat_heterogeneity = factor(spat_heterogeneity, 
                                     levels = c(0, 0.1, 1, 1000), 
                                     labels = c("temporal\nvariation", "spatiotemporal\nvariation", "spatiotemporal\nvariation (more spatial)", "spatial\nvariation")))


# div_long |> 
#   filter(disturb_rate == 0.0) |> 
#   ggplot(aes(x = disp_rate, y = diversity, color = kernel_exp)) + 
#   geom_point(alpha = 0.2) + 
#   geom_smooth(method = "gam", formula = y ~ s(x, k = 5, bs = "cs"), se = F) +
#   scale_x_log10() +
#   scale_color_viridis_d(option = "B", end = .9) +
#   facet_grid(div_type ~ spat_heterogeneity, scales = "free_y") +
#   labs(x = "Emigration rate", y = "Diversity", color = "Kernel exponent")

fig_diversity_nodisturb <- div_long |> 
  filter(disturb_rate == 0.0) |> 
  ggplot(aes(x = disp_rate, y = diversity, color = kernel_exp)) + 
  geom_point(alpha = 0.2) + 
  geom_line(data = div_long |> 
              filter(disturb_rate == 0.0) |> 
              group_by(disp_rate, kernel_exp, spat_heterogeneity, div_type) |> 
              summarize(diversity = mean(diversity))) +
  scale_x_log10() +
  scale_color_viridis_d(option = "B", end = .9) +
  facet_grid(div_type ~ spat_heterogeneity, scales = "free_y") +
  labs(x = "Emigration rate", y = "Diversity", color = "Kernel exponent")


# div_long |> 
#   filter(disturb_rate == 0.01) |> 
#   ggplot(aes(x = disp_rate, y = diversity, color = kernel_exp)) + 
#   geom_point(alpha = 0.2) + 
#   geom_smooth(method = "gam", formula = y ~ s(x, k = 5, bs = "cs"), se = F) +
#   scale_x_log10() +
#   scale_color_viridis_d(option = "B", end = .9) +
#   facet_grid(div_type ~ spat_heterogeneity, scales = "free_y") +
#   labs(x = "Emigration rate", y = "Diversity", color = "Kernel exponent")

fig_diversity_disturb <- div_long |> 
  filter(disturb_rate == 0.01) |> 
  ggplot(aes(x = disp_rate, y = diversity, color = kernel_exp)) + 
  geom_point(alpha = 0.2) + 
  geom_line(data = div_long |> 
              filter(disturb_rate == 0.01) |> 
              group_by(disp_rate, kernel_exp, spat_heterogeneity, div_type) |> 
              summarize(diversity = mean(diversity))) +
  scale_x_log10() +
  scale_color_viridis_d(option = "B", end = .9) +
  facet_grid(div_type ~ spat_heterogeneity, scales = "free_y") +
  labs(x = "Emigration rate", y = "Diversity", color = "Kernel exponent")


ggsave(paste0("figures/",comp_scenario,"_diversity_nodisturb.pdf"), plot = fig_diversity_nodisturb, width = 8, height = 6)
ggsave(paste0("figures/",comp_scenario,"_diversity_disturb.pdf"), plot = fig_diversity_disturb, width = 8, height = 6)
ggsave(paste0("figures/",comp_scenario,"_diversity_nodisturb.png"), plot = fig_diversity_nodisturb, width = 8, height = 6, dpi = 500, bg = "white")
ggsave(paste0("figures/",comp_scenario,"_diversity_disturb.png"), plot = fig_diversity_disturb, width = 8, height = 6, dpi = 500, bg = "white")



# variability
var_long <- variability |> 
  filter(disturb_rate <= 0.01) |> 
  filter(condition == comp_scenario) |> 
  pivot_longer(cols = c(CV_S_L, CV_C_L, CV_S_R, CV_C_R), 
               names_to = "variability", values_to = "CV") |> 
  pivot_longer(cols = c(phi_S_L2R, phi_C_L2R, phi_S2C_L, phi_S2C_R), 
               names_to = "synchrony", values_to = "phi") |> 
  mutate(variability = factor(variability, 
                              levels = c("CV_S_L", "CV_C_L", "CV_S_R", "CV_C_R"),
                              labels = c("population", "community", "metapopulation", "metacommunity")),
         synchrony = factor(synchrony, 
                            levels = c("phi_S2C_L", "phi_S2C_R", "phi_S_L2R", "phi_C_L2R"),
                            labels = c(str_wrap("local species synchrony", width = 15),
                                       str_wrap("regional species synchrony", width = 15),
                                       str_wrap("spatial species synchrony", width = 15),
                                       str_wrap("spatial community synchrony", width = 15)))) |> 
  mutate(kernel_exp = as.factor(kernel_exp),
         spat_heterogeneity = factor(spat_heterogeneity, 
                                     levels = c(0, 0.1, 1, 1000), 
                                     labels = c("temporal\nvariation", "spatiotemporal\nvariation", "spatiotemporal\nvariation (more spatial)", "spatial\nvariation")))




# Plot figures

# First, we'll just look at variability at different scales


# var_long |> 
#   filter(disturb_rate == 0.01, gamma_div > 0) |> 
#   ggplot(aes(x = disp_rate, y = 1/CV, color = kernel_exp)) + 
#   geom_point(alpha = 0.2) + 
#   #geom_smooth(method = "loess", se = FALSE) +
#   geom_smooth(method = "gam", formula = y ~ s(x, k = 5, bs = "cs"), se = FALSE) + 
#   facet_grid(variability~spat_heterogeneity, scales = "free_y") +
#   scale_x_log10() +
#   scale_y_log10() +
#   scale_color_viridis_d(option = "B", end = .9) +
#   labs(color = "Dispersal kernel \nexponent",
#        x = "Emigration rate", y = "Stability (1/CV)")

fig_var_disturb <- var_long |> 
  filter(disturb_rate == 0.01, gamma_div > 0) |> 
  ggplot(aes(x = disp_rate, y = 1/CV, color = kernel_exp)) + 
  geom_point(alpha = 0.2) + 
  geom_line(data = var_long |> 
              filter(disturb_rate == 0.01, gamma_div > 0) |> 
              group_by(disp_rate, kernel_exp, spat_heterogeneity, variability) |> 
              summarize(CV = mean(CV))) +
  facet_grid(variability~spat_heterogeneity, scales = "free_y") +
  scale_x_log10() +
  #scale_y_log10() +
  scale_color_viridis_d(option = "B", end = .9) +
  labs(color = "Dispersal kernel \nexponent",
       x = "Emigration rate", y = "Stability (1/CV)")

# var_long |> 
#   filter(disturb_rate == 0, gamma_div > 0) |> 
#   ggplot(aes(x = disp_rate, y = 1/CV, color = kernel_exp)) + 
#   geom_point(alpha = 0.2) + 
#   #geom_smooth(method = "loess", se = FALSE) +
#   geom_smooth(method = "gam", formula = y ~ s(x, k = 5, bs = "cs"), se = FALSE) + 
#   facet_grid(variability~spat_heterogeneity, scales = "free_y") +
#   scale_x_log10() +
#   scale_y_log10() +
#   scale_color_viridis_d(option = "B", end = .9) +
#   labs(color = "Dispersal kernel \nexponent",
#        x = "Emigration rate", y = "Stability (1/CV)")
fig_var_nodisturb <- var_long |> 
  filter(disturb_rate == 0, gamma_div > 0) |> 
  ggplot(aes(x = disp_rate, y = 1/CV, color = kernel_exp)) + 
  geom_point(alpha = 0.2) + 
  geom_line(data = var_long |> 
              filter(disturb_rate == 0.0, gamma_div > 0) |> 
              group_by(disp_rate, kernel_exp, spat_heterogeneity, variability) |> 
              summarize(CV = mean(CV))) +
  facet_grid(variability~spat_heterogeneity, scales = "free_y") +
  scale_x_log10() +
  #scale_y_log10() +
  scale_color_viridis_d(option = "B", end = .9) +
  labs(color = "Dispersal kernel \nexponent",
       x = "Emigration rate", y = "Stability (1/CV)")

ggsave(paste0("figures/",comp_scenario,"_variability_nodisturb.pdf"), plot = fig_var_nodisturb, width = 8, height = 6)
ggsave(paste0("figures/",comp_scenario,"_variability_disturb.pdf"), plot = fig_var_disturb, width = 8, height = 6)
ggsave(paste0("figures/",comp_scenario,"_variability_nodisturb.png"), plot = fig_var_nodisturb, width = 8, height = 6, dpi = 500, bg = "white")
ggsave(paste0("figures/",comp_scenario,"_variability_disturb.png"), plot = fig_var_disturb, width = 8, height = 6, dpi = 500, bg = "white")



# variability |> 
#   filter(disturb_rate == 0) |> 
#   ggplot(aes(x = disp_rate, y = 1/CV_C_R, color = as.factor(kernel_exp))) + 
#   geom_point(alpha = 0.3) + 
#   geom_smooth(method = "lm", se=FALSE) +
#   geom_smooth(aes(color = NULL), method = "lm", se=FALSE, color = "black") +
#   facet_grid(spat_heterogeneity~condition, scales = "free") +
#   scale_x_log10() +
#   labs(x = "Emigration rate", y = "Metacommunity stability (1/CV)")
# metacommunities overall are quite stable when there's spatial heterogeneity
# but this is enhanced by stable species interactions too

# variability |> 
#   filter(disturb_rate == 0) |> 
#   ggplot(aes(x = disp_rate, y = 1/CV_C_L, color = as.factor(kernel_exp))) + 
#   geom_point(alpha = 0.3) + 
#   geom_smooth(method = "lm", se=FALSE) +
#   geom_smooth(aes(color = NULL), method = "lm", se=FALSE, color = "black") +
#   facet_grid(spat_heterogeneity~condition, scales = "free") +
#   scale_x_log10() +
#   labs(x = "Emigration rate", y = "Local stability (1/CV)")
# local stability enhanced by emigration rate and flat kernels in non-spatial environments
# but reduced by emigration and steep kernels in spatial environments

# now look at the phi values, the scaling coefficients


# var_long |> 
#   filter(disturb_rate == 0.01, gamma_div > 0) |> 
#   ggplot(aes(x = disp_rate, y = phi, color = kernel_exp)) + 
#   geom_point(alpha = 0.2) + 
#   geom_smooth(method = "gam", formula = y ~ s(x, k = 5, bs = "cs"), se = FALSE) + 
#   facet_grid(synchrony~spat_heterogeneity, scales = "free_y") +
#   scale_x_log10() +
#   scale_y_continuous(limits = c(0,1)) +
#   scale_color_viridis_d(option = "B", end = .9) +
#   labs(color = "Dispersal kernel \nexponent",
#        x = "Emigration rate", y = "Synchrony (phi)") 

fig_phi_disturb <- var_long |> 
  filter(disturb_rate == 0.01, gamma_div > 0) |> 
  ggplot(aes(x = disp_rate, y = phi, color = kernel_exp)) + 
  geom_point(alpha = 0.2) + 
  geom_line(data = var_long |> 
              filter(disturb_rate == 0.01, gamma_div > 0) |> 
              group_by(disp_rate, kernel_exp, spat_heterogeneity, synchrony) |> 
              summarize(phi = mean(phi))) +
  facet_grid(synchrony~spat_heterogeneity, scales = "free_y") +
  scale_x_log10() +
  scale_y_continuous(limits = c(0,1)) +
  scale_color_viridis_d(option = "B", end = .9) +
  labs(color = "Dispersal kernel \nexponent",
       x = "Emigration rate", y = "Synchrony (phi)") 


# var_long |> 
#   filter(disturb_rate == 0, gamma_div > 0) |> 
#   ggplot(aes(x = disp_rate, y = phi, color = kernel_exp)) + 
#   geom_point(alpha = 0.2) + 
#   geom_smooth(method = "gam", formula = y ~ s(x, k = 5, bs = "cs"), se = FALSE) + 
#   facet_grid(synchrony~spat_heterogeneity, scales = "free_y") +
#   scale_x_log10() +
#   scale_color_viridis_d(option = "B", end = .9) +
#   labs(color = "Dispersal kernel \nexponent",
#        x = "Emigration rate", y = "Synchrony (phi)") 
# 

fig_phi_nodisturb <- var_long |> 
  filter(disturb_rate == 0, gamma_div > 0) |> 
  ggplot(aes(x = disp_rate, y = phi, color = kernel_exp)) + 
  geom_point(alpha = 0.2) + 
  geom_line(data = var_long |> 
              filter(disturb_rate == 0.0, gamma_div > 0) |> 
              group_by(disp_rate, kernel_exp, spat_heterogeneity, synchrony) |> 
              summarize(phi = mean(phi))) +
  facet_grid(synchrony~spat_heterogeneity, scales = "free_y") +
  scale_x_log10() +
  scale_color_viridis_d(option = "B", end = .9) +
  labs(color = "Dispersal kernel \nexponent",
       x = "Emigration rate", y = "Synchrony (phi)") 

ggsave(paste0("figures/",comp_scenario,"_synchrony_nodisturb.pdf"), plot = fig_phi_nodisturb, width = 8, height = 6)
ggsave(paste0("figures/",comp_scenario,"_synchrony_disturb.pdf"), plot = fig_phi_disturb, width = 8, height = 6)
ggsave(paste0("figures/",comp_scenario,"_synchrony_nodisturb.png"), plot = fig_phi_nodisturb, width = 8, height = 6, dpi = 500, bg = "white")
ggsave(paste0("figures/",comp_scenario,"_synchrony_disturb.png"), plot = fig_phi_disturb, width = 8, height = 6, dpi = 500, bg = "white")

# now diversity-stability relationships

# variability |> 
#   filter(condition == comp_scenario, disturb_rate == 0) |> 
#   mutate(spat_heterogeneity = as.factor(spat_heterogeneity)) |> 
#   ggplot(aes(x = alpha_div, y = 1/CV_C_R, color = spat_heterogeneity)) +
#   geom_point(alpha = 0.2) +
#   geom_smooth(method = 'lm') + 
#   scale_y_log10() +
#   labs(x = "gamma diversity", 
#        y = "Metacommunity stability (1/CV)")
# 
# variability |> 
#   filter(condition == comp_scenario, disturb_rate == 0) |> 
#   mutate(spat_heterogeneity = as.factor(spat_heterogeneity)) |> 
#   ggplot(aes(x = beta_div, y = 1/CV_C_R, color = spat_heterogeneity)) +
#   geom_point(alpha = 0.2) +
#   geom_smooth(method = 'lm') + 
#   scale_y_log10() +
#   labs(x = "beta diversity", 
#        y = "Metacommunity stability (1/CV)")
# 
# variability |> 
#   filter(condition == comp_scenario, disturb_rate == 0) |> 
#   mutate(spat_heterogeneity = as.factor(spat_heterogeneity)) |> 
#   ggplot(aes(x = gamma_div, y = 1/CV_C_R, color = spat_heterogeneity)) +
#   geom_point(alpha = 0.2) +
#   geom_smooth(method = 'lm') + 
#   scale_y_log10() +
#   labs(x = "alpha diversity", 
#        y = "Metacommunity stability (1/CV)")


# species-env matches
# load and combine reps
patches_over_time <- data.table()
i = 1
for(file in list.files(path = subfolder, pattern = "temp_per_patch\\.csv$")){
  this_run <- fread(paste(subfolder, file, sep = "/"))
  print(paste("Read file:", file))
  this_run$rep <- i
  patches_over_time <- bind_rows(patches_over_time, this_run)
  i = i + 1
} 

patches_over_time$species <- as.factor(patches_over_time$species)
patches_over_time$rep <- as.factor(patches_over_time$rep)
patches_over_time$patch <- as.factor(patches_over_time$patch)
patches_over_time$kernel_exp <- as.factor(patches_over_time$kernel_exp)
patches_over_time <- patches_over_time |> 
  filter(spat_heterogeneity != 1) |> 
  mutate(kernel_exp = as.factor(kernel_exp),
         spat_heterogeneity = factor(spat_heterogeneity, 
                                     levels = c(0, 0.1, 1, 1000), 
                                     labels = c("temporal\nvariation", "spatiotemporal\nvariation", "spatiotemporal\nvariation (more spatial)", "spatial\nvariation")))


# analyze environmental filtering
patches_over_time |> 
  filter(comp == comp_scenario, extirp_prob == 0) |> 
  group_by(emigration, kernel_exp, species, spat_heterogeneity) |> 
  summarize(env_mismatch_mean = mean(env_mismatch),
            env_mismatch_sd = sd(env_mismatch)) |> 
  ggplot(aes(x = emigration, y = env_mismatch_mean, ymin = env_mismatch_mean - env_mismatch_sd,
             ymax = env_mismatch_mean + env_mismatch_sd, color = species, shape = spat_heterogeneity)) + 
  geom_point(alpha = 0.5) +
  geom_errorbar(alpha = 0.5) +
  geom_line() +
  scale_x_log10() +
  scale_fill_viridis_c() +
  facet_grid(species~kernel_exp, scales = "free_y") +
  labs(x = "Emigration rate", y = "Environmental mismatch (mean ± s.d.)")

fig_env_mismatch <- patches_over_time |> 
  filter(comp == comp_scenario, extirp_prob == 0) |> 
  group_by(emigration, kernel_exp, spat_heterogeneity) |> 
  summarize(env_mismatch_mean = mean(env_mismatch),
            env_mismatch_sd = sd(env_mismatch)) |> 
  ggplot(aes(x = emigration, y = env_mismatch_mean, ymin = env_mismatch_mean - env_mismatch_sd,
             ymax = env_mismatch_mean + env_mismatch_sd, color = kernel_exp)) + 
  geom_point(alpha = 0.5) +
  # geom_errorbar(alpha = 0.5) +
  geom_line() +
  scale_x_log10() +
  scale_color_viridis_d(option = "B", end = .9) +
  facet_grid(. ~ spat_heterogeneity) +
  labs(x = "Emigration rate", 
       y = "Env. mismatch (mean)",
       color = "Dispersal kernel\n exponent")
ggsave(filename = paste0("figures/",comp_scenario,"_env-mismatch_nodisturb.png"),
       plot = fig_env_mismatch, width = 6, height = 3, dpi = 500, bg = "white")
ggsave(filename = paste0("figures/",comp_scenario,"_env-mismatch_nodisturb.pdf"),
       plot = fig_env_mismatch, width = 6, height = 3)
  
fig_env_sinks <-  patches_over_time |> 
  filter(comp == comp_scenario, extirp_prob == 0.0) |> 
  group_by(emigration, kernel_exp, spat_heterogeneity) |> 
  summarize(env_sink_mean = mean(sink_env),
            env_sink_sd = sd(sink_env)) |> 
  ggplot(aes(x = emigration, y = env_sink_mean, ymin = env_sink_mean - env_sink_sd,
             ymax = env_sink_mean + env_sink_sd, color = kernel_exp)) + 
  geom_point(alpha = 0.5) +
  # geom_errorbar(alpha = 0.5) +
  geom_line() +
  scale_x_log10() +
  scale_color_viridis_d(option = "B", end = .9) +
  facet_grid(. ~ spat_heterogeneity) +
  labs(x = "Emigration rate", 
       y = "Env. sinks (mean)",
       color = "Dispersal kernel\n exponent")
ggsave(filename = paste0("figures/",comp_scenario,"_env-sinks_nodisturb.png"),
       plot = fig_env_sinks, width = 6, height = 3, dpi = 500, bg = "white")
ggsave(filename = paste0("figures/",comp_scenario,"_env-sinks_nodisturb.pdf"),
       plot = fig_env_sinks, width = 6, height = 3)

fig_env_costs <- patches_over_time |> 
  filter(comp == comp_scenario, extirp_prob == 0) |> 
  mutate(tot_comp = (delta_env) / abundance_mean) |> 
  group_by(emigration, kernel_exp, spat_heterogeneity) |> 
  summarize(comp_mean = mean(tot_comp)) |> 
  ggplot(aes(x = emigration, y = comp_mean, color = kernel_exp)) + 
  geom_point(alpha = 0.5) +
  geom_line() +
  scale_x_log10() +
  scale_color_viridis_d(option = "B", end = .9) +
  facet_grid(. ~ spat_heterogeneity) +
  labs(x = "Emigration rate", 
       y = "Mean fitness effects \n of env. filtering",
       color = "Dispersal kernel\n exponent")
ggsave(filename = paste0("figures/",comp_scenario,"_env-costs_nodisturb.png"),
       plot = fig_env_costs, width = 6, height = 3, dpi = 500, bg = "white")
ggsave(filename = paste0("figures/",comp_scenario,"_env-costs_nodisturb.pdf"),
       plot = fig_env_costs, width = 6, height = 3)

fig_env_filtering <- fig_env_costs + fig_env_mismatch + fig_env_sinks +
  plot_layout(nrow = 3, guides = "collect")
ggsave(filename = paste0("figures/",comp_scenario,"_env-filtering_nodisturb.png"),
       plot = fig_env_filtering, width = 6, height = 6, dpi = 500, bg = "white")
ggsave(filename = paste0("figures/",comp_scenario,"_env-filtering_nodisturb.pdf"),
       plot = fig_env_filtering, width = 6, height = 6)



# analyze biotic filtering
fig_competition_effects <- patches_over_time |> 
  filter(comp == comp_scenario, extirp_prob == 0) |> 
  mutate(intra_comp = (delta_bio_intra)/abundance_mean,
         inter_comp = (delta_bio_inter)/abundance_mean) |> 
  group_by(emigration, kernel_exp, spat_heterogeneity) |> 
  summarize(intra_comp_mean = mean(intra_comp),
            inter_comp_mean = mean(inter_comp)) |> 
  ggplot(aes(x = emigration, color = kernel_exp)) + 
  geom_point(aes(y = intra_comp_mean), alpha = 0.5) +
  geom_line(aes(y = intra_comp_mean), linetype = "solid") +
  
  geom_point(aes(y = inter_comp_mean), alpha = 0.5) +
  geom_line(aes(y = inter_comp_mean), linetype = "dashed") +
  
  scale_x_log10() +
  scale_color_viridis_d(option = "B", end = .9) +
  facet_grid(. ~ spat_heterogeneity) +
  labs(x = "Emigration rate", 
       y = "Mean fitness effects of \n intra (solid) + inter (dashed) competition",
       color = "Dispersal kernel\n exponent")

fig_comp_ratios <- patches_over_time |> 
  filter(comp == comp_scenario, extirp_prob == 0.00) |> 
  mutate(intra_comp = (delta_bio_intra)/abundance_mean,
         inter_comp = (delta_bio_inter)/abundance_mean) |> 
  mutate(comp_ratio = inter_comp/intra_comp) |> 
  filter(!is.na(comp_ratio)) |> 
  group_by(emigration, kernel_exp, spat_heterogeneity, species) |> 
  summarize(comp_ratio = mean(comp_ratio)) |> 
  ggplot(aes(x = emigration, y = comp_ratio, color = kernel_exp, fill = kernel_exp)) + 
  geom_hline(yintercept = 1) +
  geom_point(alpha = 0.5) +
  geom_smooth(se=F) +
  #geom_line() +
  #scale_y_log10() +
  scale_x_log10() +
  scale_color_viridis_d(option = "B", end = .9) +
  scale_fill_viridis_d(option = "B", end = .9) +
  facet_grid(. ~ spat_heterogeneity) +
  labs(x = "Emigration rate", 
       y = "Interspecific : intraspecific\ncompetition ratio",
       color = "Dispersal kernel\n exponent",
       fill = "Dispersal kernel\n exponent")

ggsave(filename = paste0("figures/",comp_scenario,"_comp-ratio_nodisturb.png"),
       plot = fig_comp_ratios, width = 6, height = 4, dpi = 500, bg = "white")
ggsave(filename = paste0("figures/",comp_scenario,"_comp-ratio_nodisturb.pdf"),
       plot = fig_comp_ratios, width = 6, height = 4)


patches_over_time |> 
  filter(comp == comp_scenario, extirp_prob == 0) |> 
  #mutate(tot_comp = (W_max + delta_bio) / abundance_mean) |> 
  mutate(tot_comp = (delta_bio) / abundance_mean) |> 
  group_by(emigration, kernel_exp, spat_heterogeneity) |> 
  summarize(comp_mean = mean(tot_comp)) |> 
  ggplot(aes(x = emigration, y = comp_mean, color = kernel_exp)) + 
  geom_point(alpha = 0.5) +
  geom_line() +
  scale_x_log10() +
  scale_color_viridis_d(option = "B", end = .9) +
  facet_grid(. ~ spat_heterogeneity) +
  labs(x = "Emigration rate", 
       y = "Mean fitness effects \n of competition",
       color = "Dispersal kernel\n exponent")

# now understand demographic effects
fig_demo_sink_effects <- patches_over_time |> 
  filter(comp == comp_scenario, extirp_prob == 0) |> 
  #mutate(demo_stoch = (W_max + delta_stoch_demo) / abundance_mean) |> 
  mutate(demo_sink = sink_stoch_demo) |> 
  group_by(emigration, kernel_exp, spat_heterogeneity) |> 
  summarize(demo_sink = mean(demo_sink)) |> 
  ggplot(aes(x = emigration, y = demo_sink, color = kernel_exp)) + 
  geom_point(alpha = 0.5) +
  geom_line() +
  scale_x_log10() +
  scale_color_viridis_d(option = "B", end = .9) +
  facet_grid(. ~ spat_heterogeneity) +
  labs(x = "Emigration rate", 
       y = "Sink populations due to \n demographic stochasticity",
       color = "Dispersal kernel\n exponent")

fig_demo_extinctions <- patches_over_time |> 
  filter(comp == comp_scenario, extirp_prob == 0) |> 
  group_by(emigration, kernel_exp, spat_heterogeneity) |> 
  summarize(demo_extinctions = mean(stoch_extinct_demo)) |> 
  ggplot(aes(x = emigration, y = demo_extinctions, color = kernel_exp)) + 
  geom_point(alpha = 0.5) +
  geom_line() +
  scale_x_log10() +
  scale_color_viridis_d(option = "B", end = .9) +
  facet_grid(. ~ spat_heterogeneity) +
  labs(x = "Emigration rate", 
       y = "Local extinctions due to \n demographic stochasticity",
       color = "Dispersal kernel\n exponent")

fig_demog_effects <- fig_demo_extinctions + 
  fig_demo_sink_effects + 
  plot_layout(nrow = 2, guides = "collect")

ggsave(filename = paste0("figures/",comp_scenario,"_demographic_nodisturb.png"),
       plot = fig_demog_effects, width = 6, height = 4, dpi = 500, bg = "white")
ggsave(filename = paste0("figures/",comp_scenario,"_demographic_nodisturb.pdf"),
       plot = fig_demog_effects, width = 6, height = 4)

# focus on dispersal
fig_dispersal_effects <- patches_over_time |> 
  filter(comp == comp_scenario, extirp_prob == 0) |> 
  #mutate(demo_stoch = (W_max + delta_stoch_demo) / abundance_mean) |> 
  mutate(dispersal_effect = (delta_dispersal) / abundance_mean) |> 
  group_by(emigration, kernel_exp, spat_heterogeneity) |> 
  summarize(dispersal_effect_mean = mean(dispersal_effect),
            dispersal_effect_sd = sd(dispersal_effect)) |> 
  ggplot(aes(x = emigration, y = dispersal_effect_mean,
             # ymin = dispersal_effect_mean-dispersal_effect_sd,
             # ymax = dispersal_effect_mean+dispersal_effect_sd,
             color = kernel_exp)) + 
  geom_point(alpha = 0.5) +
  # geom_errorbar(alpha = .5) +
  geom_line() +
  scale_x_log10() +
  scale_color_viridis_d(option = "B", end = .9) +
  facet_grid(.~ spat_heterogeneity) +
  labs(x = "Emigration rate", 
       y = "Mean fitness effects of \n dispersal",
       color = "Dispersal kernel\n exponent")

patches_over_time |> 
  filter(comp == comp_scenario, extirp_prob == 0) |> 
  group_by(emigration, kernel_exp, spat_heterogeneity) |> 
  summarize(dispersal_extinctions = mean(extinct_by_emigration),
            dispersal_rescues = sum(rescued_by_dispersal)) |>
  ggplot(aes(x = emigration, y = dispersal_rescues+1, color = kernel_exp)) + 
  geom_point(alpha = 0.5) +
  geom_line() + 
  scale_x_log10() +
  scale_y_log10() +
  scale_color_viridis_d(option = "B", end = .9) +
  facet_grid(.~ spat_heterogeneity) +
  labs(x = "Emigration rate", 
       y = "Dispersal-induced rescues",
       color = "Dispersal kernel\n exponent")

fig_dispersal_fitness_effects <- patches_over_time |> 
  filter(comp == comp_scenario, extirp_prob == 0) |> 
  #mutate(demo_stoch = (W_max + delta_stoch_demo) / abundance_mean) |> 
  mutate(dispersal_effect = (delta_dispersal) / abundance_mean) |> 
  group_by(emigration, kernel_exp, spat_heterogeneity, extirp_prob) |> 
  summarize(dispersal_effect_mean = mean(dispersal_effect),
            dispersal_effect_sd = sd(dispersal_effect)) |> 
  ggplot(aes(x = emigration, y = dispersal_effect_mean,
             # ymin = dispersal_effect_mean-dispersal_effect_sd,
             # ymax = dispersal_effect_mean+dispersal_effect_sd,
             color = kernel_exp)) + 
  geom_point(alpha = 0.5) +
  # geom_errorbar(alpha = .5) +
  geom_line() +
  scale_x_log10() +
  scale_color_viridis_d(option = "B", end = .9) +
  facet_grid(. ~ spat_heterogeneity) +
  labs(x = "Emigration rate", 
       y = "Mean fitness effects of dispersal",
       color = "Dispersal kernel\n exponent")
ggsave(filename = paste0("figures/",comp_scenario,"_dispersal_effects.png"),
       plot = fig_dispersal_fitness_effects, width = 6, height = 3, dpi = 500, bg = "white")
ggsave(filename = paste0("figures/",comp_scenario,"_dispersal_effects.pdf"),
       plot = fig_dispersal_fitness_effects, width = 6, height = 3)


## spatial
space_by_time <- data.table()
i = 1
for(file in list.files(path = subfolder, pattern = "spat_per_time\\.csv$")){
  this_run <- fread(paste(subfolder, file, sep = "/"))
  print(paste("Read file:", file))
  this_run$rep <- i
  space_by_time <- bind_rows(space_by_time, this_run)
  i = i + 1
} 

space_by_time$species <- as.factor(space_by_time$species)
space_by_time$rep <- as.factor(space_by_time$rep)
space_by_time <- space_by_time |> 
  filter(spat_heterogeneity != 1) |> 
  mutate(kernel_exp = as.factor(kernel_exp),
         spat_heterogeneity = factor(spat_heterogeneity, 
                                     levels = c(0, 0.1, 1, 1000), 
                                     labels = c("temporal\nvariation", "spatiotemporal\nvariation", "spatiotemporal\nvariation (more spatial)", "spatial\nvariation")))


fig_occupancy_var <- space_by_time |> 
  group_by(species, spat_heterogeneity, emigration, kernel_exp) |> 
  summarize(occupancy_sd = sd(occupancy)) |> 
  ggplot(aes(x = emigration, y = occupancy_sd, grouping = species, color = species)) + 
  geom_point(alpha = 0.5) + 
  geom_line() +
  scale_x_log10() +
  scale_color_viridis_d(option = "H") +
  facet_grid(kernel_exp ~ spat_heterogeneity) +
  labs(x = "Emigration rate", 
       y = "Occupancy variability (s.d.)",
       color = "Species")
ggsave(filename = paste0("figures/",comp_scenario,"_occupancy-var.png"),
       plot = fig_occupancy_var, width = 6, height = 6, dpi = 500, bg = "white")
ggsave(filename = paste0("figures/",comp_scenario,"_occupancy-var.pdf"),
       plot = fig_occupancy_var, width = 6, height = 6)


fig_spatial_abund_var <- space_by_time |> 
  filter(comp == comp_scenario) |> 
  group_by(species, spat_heterogeneity, emigration, kernel_exp) |> 
  summarize(abundance_sd_mean = mean(sqrt(abundance_var))) |> 
  ggplot(aes(x = emigration, y = abundance_sd_mean, grouping = species, color = species)) + 
  geom_point(alpha = 0.5) + 
  geom_line() +
  scale_x_log10() +
  scale_color_viridis_d(option = "H") +
  facet_grid(kernel_exp ~ spat_heterogeneity) +
  labs(x = "Emigration rate", 
       y = "Mean spatial s.d. in abundance",
       color = "Species")
ggsave(filename = paste0("figures/",comp_scenario,"_spatial_abund-var.png"),
       plot = fig_spatial_abund_var, width = 6, height = 6, dpi = 500, bg = "white")
ggsave(filename = paste0("figures/",comp_scenario,"_spatial_abund-var.pdf"),
       plot = fig_spatial_abund_var, width = 6, height = 6)

space_by_time |> 
  group_by(species, spat_heterogeneity, emigration, kernel_exp) |> 
  summarize(occupancy = mean(occupancy), occupancy_sd = sd(occupancy)) |> 
  ggplot(aes(x = emigration, y = occupancy, color = as.factor(kernel_exp))) + 
  geom_point() + 
  geom_errorbar(aes(ymin = occupancy - occupancy_sd, ymax = occupancy + occupancy_sd)) + 
  geom_line() +
  scale_x_log10() +
  facet_grid(spat_heterogeneity ~ species) 

space_by_time |> 
  group_by(species, spat_heterogeneity, emigration, kernel_exp) |> 
  summarize(sink_env = mean(sink_env),
            sink_biotic = mean(sink_biotic),
            sink_stoch_demo = mean(sink_stoch_demo)) |> 
  ggplot(aes(x = emigration, y = sink_stoch_demo, color = species)) + 
  geom_point() + 
  #geom_errorbar(aes(ymin = occupancy - occupancy_sd, ymax = occupancy + occupancy_sd)) + 
  geom_line() + 
  scale_x_log10() +
  facet_grid(spat_heterogeneity ~ kernel_exp) 


space_by_time |> 
  group_by(species, time) |> 
  summarize(occupancy = mean(occupancy)) |> 
  ggplot(aes(x = time, y = occupancy, color = species)) + 
  geom_point() + 
  #geom_errorbar(aes(ymin = occupancy - occupancy_sd, ymax = occupancy + occupancy_sd)) + 
  geom_line()
