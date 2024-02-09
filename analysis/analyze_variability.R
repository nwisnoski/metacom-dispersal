library(tidyverse)
library(here)
library(patchwork)
library(data.table)
library(lme4)
library(lmerTest)

theme_set(theme_bw())

subfolder <- "sim_output/2024-02-09_sims/"
comp_scenario <- "equal"

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
         env_temp_cv_mean = replace_na(env_temp_cv_mean, 0))
  

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
         div_type = factor(div_type, levels = c("alpha_div", "gamma_div", "beta_spatial", "beta_temporal")))


div_long |> 
  filter(disturb_rate == 0.0) |> 
  ggplot(aes(x = disp_rate, y = diversity, color = as.factor(round(kernel_exp,2)))) + 
  geom_point(alpha = 0.2) + 
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5, bs = "cs"), se = F) +
  scale_x_log10() +
  scale_color_viridis_d() +
  facet_grid(div_type ~ round(spat_heterogeneity,2), scales = "free_y") +
  labs(x = "Emigration rate", y = "Diversity", color = "Kernel exponent")


div_long |> 
  filter(disturb_rate == 0.01) |> 
  ggplot(aes(x = disp_rate, y = diversity, color = as.factor(round(kernel_exp,2)))) + 
  geom_point(alpha = 0.2) + 
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5, bs = "cs"), se = F) +
  scale_x_log10() +
  scale_color_viridis_d() +
  facet_grid(div_type ~ round(spat_heterogeneity,2), scales = "free_y") +
  labs(x = "Emigration rate", y = "Diversity", color = "Kernel exponent")



# variability
var_long <- variability |> 
  filter(disturb_rate <= 0.01) |> 
  filter(condition == comp_scenario) |> 
  pivot_longer(cols = c(CV_S_L, CV_C_L, CV_S_R, CV_C_R), 
               names_to = "variability", values_to = "CV") |> 
  pivot_longer(cols = c(phi_S_L2R, phi_C_L2R, phi_S2C_L, phi_S2C_R), 
               names_to = "synchrony", values_to = "phi") |> 
  mutate(variability = factor(variability, levels = c("CV_S_L", "CV_C_L", "CV_S_R", "CV_C_R")),
         synchrony = factor(synchrony, levels = c("phi_S_L2R", "phi_S2C_L", "phi_S2C_R", "phi_C_L2R")))



# Plot figures

# First, we'll just look at variability at different scales


panel_var_disturbed <- var_long |> 
  filter(disturb_rate == 0.01, gamma_div > 0) |> 
  ggplot(aes(x = disp_rate, y = 1/CV, color = as.factor(round(kernel_exp, 2)))) + 
  geom_point(alpha = 0.2) + 
  #geom_smooth(method = "loess", se = FALSE) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5, bs = "cs"), se = FALSE) + 
  facet_grid(variability~spat_heterogeneity, scales = "free_y") +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_viridis_d(option = "A") +
  labs(color = "Dispersal kernel \nexponent",
       x = "Emigration rate", y = "Stability (1/CV)")
panel_var_disturbed

panel_var_undisturbed <- var_long |> 
  filter(disturb_rate == 0, gamma_div > 0) |> 
  ggplot(aes(x = disp_rate, y = 1/CV, color = as.factor(round(kernel_exp, 2)))) + 
  geom_point(alpha = 0.2) + 
  #geom_smooth(method = "loess", se = FALSE) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5, bs = "cs"), se = FALSE) + 
  facet_grid(variability~spat_heterogeneity, scales = "free_y") +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_viridis_d(option = "A") +
  labs(color = "Dispersal kernel \nexponent",
       x = "Emigration rate", y = "Stability (1/CV)")
panel_var_undisturbed


variability |> 
  filter(disturb_rate == 0) |> 
  ggplot(aes(x = disp_rate, y = 1/CV_C_R, color = as.factor(kernel_exp))) + 
  geom_point(alpha = 0.3) + 
  geom_smooth(method = "lm", se=FALSE) +
  geom_smooth(aes(color = NULL), method = "lm", se=FALSE, color = "black") +
  facet_grid(spat_heterogeneity~condition, scales = "free") +
  scale_x_log10() +
  labs(x = "Emigration rate", y = "Metacommunity stability (1/CV)")
# metacommunities overall are quite stable when there's spatial heterogeneity
# but this is enhanced by stable species interactions too

variability |> 
  filter(disturb_rate == 0) |> 
  ggplot(aes(x = disp_rate, y = 1/CV_C_L, color = as.factor(kernel_exp))) + 
  geom_point(alpha = 0.3) + 
  geom_smooth(method = "lm", se=FALSE) +
  geom_smooth(aes(color = NULL), method = "lm", se=FALSE, color = "black") +
  facet_grid(spat_heterogeneity~condition, scales = "free") +
  scale_x_log10() +
  labs(x = "Emigration rate", y = "Local stability (1/CV)")
# local stability enhanced by emigration rate and flat kernels in non-spatial environments
# but reduced by emigration and steep kernels in spatial environments

# now look at the phi values, the scaling coefficients


panel_phi_disturbed <- var_long |> 
  filter(disturb_rate == 0.01, gamma_div > 0) |> 
  ggplot(aes(x = disp_rate, y = phi, color = as.factor(round(kernel_exp, 2)))) + 
  geom_point(alpha = 0.2) + 
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5, bs = "cs"), se = FALSE) + 
  facet_grid(synchrony~spat_heterogeneity, scales = "free_y") +
  scale_x_log10() +
  scale_y_continuous(limits = c(0,1)) +
  scale_color_viridis_d(option = "A") +
  labs(color = "Dispersal kernel \nexponent",
       x = "Emigration rate", y = "Synchrony (phi)") 
panel_phi_disturbed

panel_phi_undisturbed <- var_long |> 
  filter(disturb_rate == 0, gamma_div > 0) |> 
  ggplot(aes(x = disp_rate, y = phi, color = as.factor(round(kernel_exp, 2)))) + 
  geom_point(alpha = 0.2) + 
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5, bs = "cs"), se = FALSE) + 
  facet_grid(synchrony~spat_heterogeneity, scales = "free_y") +
  scale_x_log10() +
  scale_color_viridis_d(option = "A") +
  labs(color = "Dispersal kernel \nexponent",
       x = "Emigration rate", y = "Synchrony (phi)") 
panel_phi_undisturbed



# now diversity-stability relationships

variability |> 
  filter(condition == comp_scenario, disturb_rate == 0) |> 
  ggplot(aes(x = alpha_div, y = 1/CV_C_R, color = as.factor(disp_rate))) + 
  geom_point(alpha = 0.3) + 
  geom_smooth(method = "lm", se=FALSE) +
  geom_smooth(aes(color = NULL), method = "lm", se=FALSE, color = "black") +
  facet_grid(as.factor(kernel_exp) ~ spat_heterogeneity, scales = "free") +
  labs(x = "Mean alpha-diversity", y = "Metacommunity stability (1/CV)")

variability |> 
  filter(condition == comp_scenario, disturb_rate == 0) |> 
  ggplot(aes(x = beta_div, y = 1/CV_C_R, color = as.factor(disp_rate))) + 
  geom_point(alpha = 0.3) + 
  geom_smooth(method = "lm", se=FALSE) +
  geom_smooth(aes(color = NULL), method = "lm", se=FALSE, color = "black") +
  facet_grid(as.factor(kernel_exp) ~ spat_heterogeneity, scales = "free") +
  labs(x = "Beta-diversity", y = "Metacommunity stability (1/CV)")

variability |> 
  filter(condition == comp_scenario, disturb_rate == 0) |> 
  ggplot(aes(x = gamma_div, y = 1/CV_C_R, color = as.factor(disp_rate))) + 
  geom_point(alpha = 0.3) + 
  geom_smooth(method = "lm", se=FALSE) +
  geom_smooth(aes(color = NULL), method = "lm", se=FALSE, color = "black") +
  facet_grid(as.factor(kernel_exp) ~ spat_heterogeneity, scales = "free") +
  labs(x = "Gamma-diversity", y = "Metacommunity stability (1/CV)")



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

patches_over_time |> 
  filter(comp == "equal", extirp_prob == 0,
         spat_heterogeneity == 1) |> 
  ggplot(aes(x = emigration, y = env_mismatch, color = as.factor(kernel_exp))) + 
  geom_point(alpha = 0.05) +
  scale_x_log10() +
  geom_smooth(formula = y ~ s(x, k = 5, bs = "cs"), se = F) +
  scale_fill_viridis_c() +
  facet_wrap(~species, scales = "free_y")
  
patches_over_time |> 
  ggplot(aes(x = colonization, y = rescued_by_dispersal, color = comp)) + 
  geom_point(alpha = 0.05) +
  geom_abline(slope = 1, intercept = 0, color = "black") +
  geom_smooth(formula = y ~ s(x, k = 4, bs = "cs"), se = F) +
  facet_grid(extirp_prob~spat_heterogeneity, scales = "free_y")
