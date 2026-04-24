library(tidyverse)
library(here)
library(patchwork)
library(data.table)
library(lme4)
library(lmerTest)

theme_set(theme_bw())

subfolder <- "sim_output/2026-04-17/"
comp_scenario <- "stable"
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

# first for 20 species (runs 1-5)
variability <- variability |> filter(rep <=5) |> 
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
         env_temp_cv_mean = replace_na(env_temp_cv_mean, 0),
         kernel_exp = round(kernel_exp, 4)) |> 
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
  mutate(#kernel_exp = as.factor(kernel_exp),
    spat_heterogeneity = factor(spat_heterogeneity, 
                                levels = c(0, 0.1, 1, 1000), 
                                labels = c("temporal\nvariation", "spatial + temporal\nvariation", "spatial + temporal\nvariation (more spatial)", "spatial\nvariation")))



# heatmaps
heat_alpha_gam_nodisturb <- div_long |> 
  filter(disturb_rate == 0.0, 
         div_type %in% c("alpha diversity", "gamma diversity")) |> 
  group_by(disp_rate, kernel_exp, spat_heterogeneity, div_type) |> 
  summarize(diversity = mean(diversity)) |> 
  ggplot(aes(x = disp_rate, y = kernel_exp, fill = diversity)) + 
  geom_tile(height = .53) + 
  scale_x_log10() +
  scale_y_log10() +
  scale_fill_viridis_c() +
  coord_fixed() +
  facet_grid(div_type ~ spat_heterogeneity) +
  labs(x = "Emigration rate", y = "Dispersal kernel exponent", fill = "Diversity", title = "No disturbance",
       subtitle = "20 species in regional species pool")

heat_alpha_gam_disturb <- div_long |> 
  filter(disturb_rate == 0.01, 
         div_type %in% c("alpha diversity", "gamma diversity")) |> 
  group_by(disp_rate, kernel_exp, spat_heterogeneity, div_type) |> 
  summarize(diversity = mean(diversity)) |> 
  ggplot(aes(x = disp_rate, y = kernel_exp, fill = diversity)) + 
  geom_tile(height = .53) + 
  scale_x_log10() +
  scale_y_log10() +
  scale_fill_viridis_c() +
  coord_fixed() +
  facet_grid(div_type ~ spat_heterogeneity) +
  labs(x = "Emigration rate", y = "Dispersal kernel exponent", fill = "Diversity", title = "Local disturbances",
       subtitle = "20 species in regional species pool")


heat_beta_nodisturb <- div_long |> 
  filter(disturb_rate == 0.0, 
         div_type %in% c("spatial beta\ndiversity", "temporal beta\ndiversity")) |> 
  group_by(disp_rate, kernel_exp, spat_heterogeneity, div_type) |> 
  summarize(diversity = mean(diversity)) |> 
  ggplot(aes(x = disp_rate, y = kernel_exp, fill = diversity)) + 
  geom_tile(height = .53) + 
  scale_x_log10() +
  scale_y_log10() +
  scale_fill_viridis_c() +
  coord_fixed() +
  facet_grid(div_type ~ spat_heterogeneity) +
  labs(x = "Emigration rate", y = "Dispersal kernel exponent", fill = "Diversity", title = "No disturbance",
       subtitle = "20 species in regional species pool")

heat_beta_disturb <- div_long |> 
  filter(disturb_rate == 0.01, 
         div_type %in% c("spatial beta\ndiversity", "temporal beta\ndiversity")) |> 
  group_by(disp_rate, kernel_exp, spat_heterogeneity, div_type) |> 
  summarize(diversity = mean(diversity)) |> 
  ggplot(aes(x = disp_rate, y = kernel_exp, fill = diversity)) + 
  geom_tile(height = .53) + 
  scale_x_log10() +
  scale_y_log10() +
  scale_fill_viridis_c() +
  coord_fixed() +
  facet_grid(div_type ~ spat_heterogeneity) +
  labs(x = "Emigration rate", y = "Dispersal kernel exponent", fill = "Diversity", title = "Local disturbances",
       subtitle = "20 species in regional species pool")

# group by scale

heat_alpha <- div_long |> 
  filter(div_type == "alpha diversity") |> 
  mutate(disturb_rate = ifelse(disturb_rate == 0, "undisturbed", "disturbed"), 
         disturb_rate = factor(disturb_rate, levels = c("undisturbed", "disturbed"))) |> 
  group_by(disp_rate, kernel_exp, spat_heterogeneity, div_type, disturb_rate) |> 
  summarize(diversity = mean(diversity)) |> 
  ggplot(aes(x = disp_rate, y = kernel_exp, fill = diversity)) + 
  geom_tile(height = .53) + 
  scale_x_log10() +
  scale_y_log10() +
  scale_fill_viridis_c() +
  coord_fixed() +
  facet_grid(disturb_rate ~ spat_heterogeneity) +
  labs(x = "Emigration rate", y = "Dispersal kernel exponent", fill = "Diversity", title = "Alpha diversity",
       subtitle = "20 species in regional species pool")

heat_gamma <- div_long |> 
  filter(div_type == "gamma diversity") |> 
  mutate(disturb_rate = ifelse(disturb_rate == 0, "undisturbed", "disturbed"), 
         disturb_rate = factor(disturb_rate, levels = c("undisturbed", "disturbed"))) |> 
  group_by(disp_rate, kernel_exp, spat_heterogeneity, div_type, disturb_rate) |> 
  summarize(diversity = mean(diversity)) |> 
  ggplot(aes(x = disp_rate, y = kernel_exp, fill = diversity)) + 
  geom_tile(height = .53) + 
  scale_x_log10() +
  scale_y_log10() +
  scale_fill_viridis_c() +
  coord_fixed() +
  facet_grid(disturb_rate ~ spat_heterogeneity) +
  labs(x = "Emigration rate", y = "Dispersal kernel exponent", fill = "Diversity", title = "Gamma diversity",
       subtitle = "20 species in regional species pool")


heat_beta_spatial <- div_long |> 
  filter(div_type == "spatial beta\ndiversity") |> 
  mutate(disturb_rate = ifelse(disturb_rate == 0, "undisturbed", "disturbed"), 
         disturb_rate = factor(disturb_rate, levels = c("undisturbed", "disturbed"))) |> 
  group_by(disp_rate, kernel_exp, spat_heterogeneity, div_type, disturb_rate) |> 
  summarize(diversity = mean(diversity)) |> 
  ggplot(aes(x = disp_rate, y = kernel_exp, fill = diversity)) + 
  geom_tile(height = .53) + 
  scale_x_log10() +
  scale_y_log10() +
  scale_fill_viridis_c() +
  coord_fixed() +
  facet_grid(disturb_rate ~ spat_heterogeneity) +
  labs(x = "Emigration rate", y = "Dispersal kernel exponent", fill = "Diversity", title = "Spatial beta diversity",
       subtitle = "20 species in regional species pool")

heat_beta_temporal <- div_long |> 
  filter(div_type == "temporal beta\ndiversity") |> 
  mutate(disturb_rate = ifelse(disturb_rate == 0, "undisturbed", "disturbed"), 
         disturb_rate = factor(disturb_rate, levels = c("undisturbed", "disturbed"))) |> 
  group_by(disp_rate, kernel_exp, spat_heterogeneity, div_type, disturb_rate) |> 
  summarize(diversity = mean(diversity)) |> 
  ggplot(aes(x = disp_rate, y = kernel_exp, fill = diversity)) + 
  geom_tile(height = .53) + 
  scale_x_log10() +
  scale_y_log10() +
  scale_fill_viridis_c() +
  coord_fixed() +
  facet_grid(disturb_rate ~ spat_heterogeneity) +
  labs(x = "Emigration rate", y = "Dispersal kernel exponent", fill = "Diversity", title = "Temporal beta diversity",
       subtitle = "20 species in regional species pool")


heat_alpha_combined <- heat_alpha_gam_nodisturb + heat_alpha_gam_disturb +
  plot_layout(ncol = 1) 
ggsave(paste0("figures/",comp_scenario,"_diversity_alpha_gamma_heatmap_20species.pdf"), plot = heat_alpha_combined, width = 8, height = 8)
ggsave(paste0("figures/",comp_scenario,"_diversity_alpha_gamma_heatmap_20species.png"), plot = heat_alpha_combined, width = 8, height = 8, dpi = 500, bg = "white")

heat_beta_combined <- heat_beta_nodisturb + heat_beta_disturb +
  plot_layout(ncol = 1) 
ggsave(paste0("figures/",comp_scenario,"_diversity_beta_heatmap_20species.pdf"), plot = heat_beta_combined, width = 8, height = 8)
ggsave(paste0("figures/",comp_scenario,"_diversity_beta_heatmap_20species.png"), plot = heat_beta_combined, width = 8, height = 8, dpi = 500, bg = "white")


heat_alpha_gamma_combined <- heat_alpha + heat_gamma +
  plot_layout(ncol = 1) +
  plot_annotation(tag_levels = "A")
ggsave(paste0("figures/",comp_scenario,"_diversity_alpha_gamma_heatmap_20species.pdf"), plot = heat_alpha_gamma_combined, width = 8, height = 8)
ggsave(paste0("figures/",comp_scenario,"_diversity_alpha_gamma_heatmap_20species.png"), plot = heat_alpha_gamma_combined, width = 8, height = 8, dpi = 500, bg = "white")

heat_beta_combined <- heat_beta_spatial + heat_beta_temporal +
  plot_layout(ncol = 1)  +
  plot_annotation(tag_levels = "A")
ggsave(paste0("figures/",comp_scenario,"_diversity_beta_heatmap_20species.pdf"), plot = heat_beta_combined, width = 8, height = 8)
ggsave(paste0("figures/",comp_scenario,"_diversity_beta_heatmap_20species.png"), plot = heat_beta_combined, width = 8, height = 8, dpi = 500, bg = "white")






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
patches_over_time <- filter(patches_over_time, rep <=5)
patches_over_time$species <- as.factor(patches_over_time$species)
patches_over_time$rep <- as.factor(patches_over_time$rep)
patches_over_time$patch <- as.factor(patches_over_time$patch)
patches_over_time$kernel_exp <- factor(signif(patches_over_time$kernel_exp, 1),
                                       levels = c(0, 1e-04, 3e-04, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1),
                                       labels = c("0", "0.0001", "0.0003", "0.001", "0.003", "0.01", "0.03", "0.1", "0.3", "1"))
patches_over_time <- patches_over_time |> 
  filter(spat_heterogeneity != 1) |> 
  mutate(kernel_exp = as.factor(kernel_exp),
         spat_heterogeneity = factor(spat_heterogeneity, 
                                     levels = c(0, 0.1, 1, 1000), 
                                     labels = c("temporal\nvariation", "spatial + temporal\nvariation", "spatial + temporal\nvariation (more spatial)", "spatial\nvariation")))


# analyze environmental filtering
# patches_over_time |> 
#   filter(comp == comp_scenario, extirp_prob == 0) |> 
#   group_by(emigration, kernel_exp, species, spat_heterogeneity) |> 
#   summarize(env_mismatch_mean = mean(env_mismatch),
#             env_mismatch_sd = sd(env_mismatch)) |> 
#   ggplot(aes(x = emigration, y = env_mismatch_mean, ymin = env_mismatch_mean - env_mismatch_sd,
#              ymax = env_mismatch_mean + env_mismatch_sd, color = species, shape = spat_heterogeneity)) + 
#   geom_point(alpha = 0.5) +
#   geom_errorbar(alpha = 0.5) +
#   geom_line() +
#   scale_x_log10() +
#   scale_fill_viridis_c() +
#   facet_grid(species~kernel_exp, scales = "free_y") +
#   labs(x = "Emigration rate", y = "Environmental mismatch (mean ± s.d.)")

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
fig_env_mismatch

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
fig_env_sinks

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
fig_env_costs


# analyze biotic filtering
fig_competition_effects <- patches_over_time |> 
  filter(comp == comp_scenario, extirp_prob == 0) |> 
  mutate(intra_comp = (delta_bio_intra)/abundance_mean,
         inter_comp = (delta_bio_inter)/abundance_mean) |> 
  group_by(emigration, kernel_exp, spat_heterogeneity) |> 
  summarize(intra_comp_mean = mean(intra_comp),
            inter_comp_mean = mean(inter_comp)) |> 
  ggplot(aes(x = emigration, color = kernel_exp)) + 
  geom_point(aes(y = intra_comp_mean), alpha = 0.5, shape = "square") +
  geom_line(aes(y = intra_comp_mean), linetype = "solid") +
  
  geom_point(aes(y = inter_comp_mean), alpha = 0.5) +
  geom_line(aes(y = inter_comp_mean), linetype = "dashed") +
  
  scale_x_log10() +
  scale_color_viridis_d(option = "B", end = .9) +
  facet_grid(. ~ spat_heterogeneity) +
  labs(x = "Emigration rate", 
       y = "Mean fitness effects of intra (solid) + \n inter (dashed) competition",
       color = "Dispersal kernel\n exponent")

fig_competition_intra <- patches_over_time |> 
  filter(comp == comp_scenario, extirp_prob == 0) |> 
  mutate(intra_comp = (delta_bio_intra)/abundance_mean,
         inter_comp = (delta_bio_inter)/abundance_mean) |> 
  group_by(emigration, kernel_exp, spat_heterogeneity) |> 
  summarize(intra_comp_mean = mean(intra_comp),
            inter_comp_mean = mean(inter_comp)) |> 
  ggplot(aes(x = emigration, color = kernel_exp)) + 
  geom_point(aes(y = intra_comp_mean), alpha = 0.5) +
  geom_line(aes(y = intra_comp_mean), linetype = "solid") +
  
  scale_x_log10() +
  scale_color_viridis_d(option = "B", end = .9) +
  facet_grid(. ~ spat_heterogeneity) +
  labs(x = "Emigration rate", 
       y = "Mean fitness effects of \n intraspecific competition",
       color = "Dispersal kernel\n exponent")

fig_competition_inter <- patches_over_time |> 
  filter(comp == comp_scenario, extirp_prob == 0) |> 
  mutate(intra_comp = (delta_bio_intra)/abundance_mean,
         inter_comp = (delta_bio_inter)/abundance_mean) |> 
  group_by(emigration, kernel_exp, spat_heterogeneity) |> 
  summarize(intra_comp_mean = mean(intra_comp),
            inter_comp_mean = mean(inter_comp)) |> 
  ggplot(aes(x = emigration, color = kernel_exp)) + 
  geom_point(aes(y = inter_comp_mean), alpha = 0.5) +
  geom_line(aes(y = inter_comp_mean), linetype = "solid") +
  
  scale_x_log10() +
  scale_color_viridis_d(option = "B", end = .9) +
  facet_grid(. ~ spat_heterogeneity) +
  labs(x = "Emigration rate", 
       y = "Mean fitness effects of \n interspecific competition",
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
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "ad"), se = F) +
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
  summarize(demo_extinctions = mean(stoch_extinct_demo),
            demo_extinctions_se = sd(stoch_extinct_demo)/5) |> 
  ggplot(aes(x = emigration, y = demo_extinctions, 
             ymin = demo_extinctions-demo_extinctions_se,
             ymax = demo_extinctions+demo_extinctions_se,
             color = kernel_exp)) + 
  geom_point(alpha = 0.5) +
  geom_errorbar() +
  geom_line() +
  scale_x_log10() +
  scale_color_viridis_d(option = "B", end = .9) +
  facet_grid(. ~ spat_heterogeneity) +
  labs(x = "Emigration rate", 
       y = "Mean stochastic \nlocal extinctions",
       color = "Dispersal kernel\n exponent",
       subtitle = "20 species in regional species pool")


# combine abiotic and biotic in one graph
fig_env_bio_filter <- 
  fig_env_costs + theme(legend.position = "null") + labs(subtitle = "20 species in regional species pool") +
  fig_competition_intra + theme(legend.position = "null") + 
  fig_competition_inter + 
  plot_layout(nrow = 3, guides = "collect") +
  plot_annotation(tag_levels = "A")
ggsave(filename = paste0("figures/",comp_scenario,"_env_bio_filtering_nodisturb_20species.png"),
       plot = fig_env_bio_filter, width = 7, height = 7, dpi = 500, bg = "white")
ggsave(filename = paste0("figures/",comp_scenario,"_env_bio_filtering_nodisturb_20species.pdf"),
       plot = fig_env_bio_filter, width = 7, height = 7, dpi = 500)

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


fig_dispersal_fitness_effects <- patches_over_time |> 
  filter(comp == comp_scenario, extirp_prob == 0) |> 
  #mutate(demo_stoch = (W_max + delta_stoch_demo) / abundance_mean) |> 
  mutate(dispersal_effect = (delta_dispersal) / abundance_mean) |> 
  group_by(emigration, kernel_exp, spat_heterogeneity, extirp_prob) |> 
  summarize(dispersal_effect_mean = mean(dispersal_effect),
            dispersal_effect_sd = sd(dispersal_effect)/10) |> 
  ggplot(aes(x = emigration, y = dispersal_effect_mean,
             ymin = dispersal_effect_mean-dispersal_effect_sd,
             ymax = dispersal_effect_mean+dispersal_effect_sd,
             color = kernel_exp)) + 
  geom_point(alpha = 0.5) +
  geom_errorbar(alpha = .5) +
  geom_line() +
  scale_x_log10() +
  scale_color_viridis_d(option = "B", end = .9) +
  facet_grid(. ~ spat_heterogeneity) +
  labs(x = "Emigration rate", 
       y = "Mean fitness effects \n of dispersal",
       color = "Dispersal kernel\n exponent", 
       subtitle = "20 species in regional species pool")

fig_demog <- fig_demo_extinctions + fig_dispersal_fitness_effects +
  plot_layout(nrow = 2, guides = "collect") +
  plot_annotation(tag_levels = "A")
ggsave(filename = paste0("figures/",comp_scenario,"_stoch_dispersal_20species.png"),
       plot = fig_demog, width = 8, height = 6, dpi = 500, bg = "white")
ggsave(filename = paste0("figures/",comp_scenario,"_stoch_dispersal_20species.pdf"),
       plot = fig_demog, width = 8, height = 6)

#### now for 75 species

# load and combine reps
variability <- data.table()
i = 1
for(file in list.files(path = subfolder, pattern = "summary\\.csv$")){
  this_run <- fread(paste(subfolder, file, sep = "/"))
  this_run$rep <- i
  variability <- bind_rows(variability, this_run)
  i = i + 1
} 

# (runs 6-10)
variability <- variability |> filter(rep >5) |> 
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
         env_temp_cv_mean = replace_na(env_temp_cv_mean, 0),
         kernel_exp = round(kernel_exp, 4)) |> 
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
  mutate(#kernel_exp = as.factor(kernel_exp),
    spat_heterogeneity = factor(spat_heterogeneity, 
                                levels = c(0, 0.1, 1, 1000), 
                                labels = c("temporal\nvariation", "spatial + temporal\nvariation", "spatial + temporal\nvariation (more spatial)", "spatial\nvariation")))



# heatmaps
heat_alpha_gam_nodisturb <- div_long |> 
  filter(disturb_rate == 0.0, 
         div_type %in% c("alpha diversity", "gamma diversity")) |> 
  group_by(disp_rate, kernel_exp, spat_heterogeneity, div_type) |> 
  summarize(diversity = mean(diversity)) |> 
  ggplot(aes(x = disp_rate, y = kernel_exp, fill = diversity)) + 
  geom_tile(height = .53) + 
  scale_x_log10() +
  scale_y_log10() +
  scale_fill_viridis_c() +
  coord_fixed() +
  facet_grid(div_type ~ spat_heterogeneity) +
  labs(x = "Emigration rate", y = "Dispersal kernel exponent", fill = "Diversity", title = "No disturbance",
       subtitle = "75 species in regional species pool")

heat_alpha_gam_disturb <- div_long |> 
  filter(disturb_rate == 0.01, 
         div_type %in% c("alpha diversity", "gamma diversity")) |> 
  group_by(disp_rate, kernel_exp, spat_heterogeneity, div_type) |> 
  summarize(diversity = mean(diversity)) |> 
  ggplot(aes(x = disp_rate, y = kernel_exp, fill = diversity)) + 
  geom_tile(height = .53) + 
  scale_x_log10() +
  scale_y_log10() +
  scale_fill_viridis_c() +
  coord_fixed() +
  facet_grid(div_type ~ spat_heterogeneity) +
  labs(x = "Emigration rate", y = "Dispersal kernel exponent", fill = "Diversity", title = "Local disturbances",
       subtitle = "75 species in regional species pool")


heat_beta_nodisturb <- div_long |> 
  filter(disturb_rate == 0.0, 
         div_type %in% c("spatial beta\ndiversity", "temporal beta\ndiversity")) |> 
  group_by(disp_rate, kernel_exp, spat_heterogeneity, div_type) |> 
  summarize(diversity = mean(diversity)) |> 
  ggplot(aes(x = disp_rate, y = kernel_exp, fill = diversity)) + 
  geom_tile(height = .53) + 
  scale_x_log10() +
  scale_y_log10() +
  scale_fill_viridis_c() +
  coord_fixed() +
  facet_grid(div_type ~ spat_heterogeneity) +
  labs(x = "Emigration rate", y = "Dispersal kernel exponent", fill = "Diversity", title = "No disturbance",
       subtitle = "75 species in regional species pool")

heat_beta_disturb <- div_long |> 
  filter(disturb_rate == 0.01, 
         div_type %in% c("spatial beta\ndiversity", "temporal beta\ndiversity")) |> 
  group_by(disp_rate, kernel_exp, spat_heterogeneity, div_type) |> 
  summarize(diversity = mean(diversity)) |> 
  ggplot(aes(x = disp_rate, y = kernel_exp, fill = diversity)) + 
  geom_tile(height = .53) + 
  scale_x_log10() +
  scale_y_log10() +
  scale_fill_viridis_c() +
  coord_fixed() +
  facet_grid(div_type ~ spat_heterogeneity) +
  labs(x = "Emigration rate", y = "Dispersal kernel exponent", fill = "Diversity", title = "Local disturbances",
       subtitle = "75 species in regional species pool")

# group by scale

heat_alpha <- div_long |> 
  filter(div_type == "alpha diversity") |> 
  mutate(disturb_rate = ifelse(disturb_rate == 0, "undisturbed", "disturbed"), 
         disturb_rate = factor(disturb_rate, levels = c("undisturbed", "disturbed"))) |> 
  group_by(disp_rate, kernel_exp, spat_heterogeneity, div_type, disturb_rate) |> 
  summarize(diversity = mean(diversity)) |> 
  ggplot(aes(x = disp_rate, y = kernel_exp, fill = diversity)) + 
  geom_tile(height = .53) + 
  scale_x_log10() +
  scale_y_log10() +
  scale_fill_viridis_c() +
  coord_fixed() +
  facet_grid(disturb_rate ~ spat_heterogeneity) +
  labs(x = "Emigration rate", y = "Dispersal kernel exponent", fill = "Diversity", title = "Alpha diversity",
       subtitle = "75 species in regional species pool")

heat_gamma <- div_long |> 
  filter(div_type == "gamma diversity") |> 
  mutate(disturb_rate = ifelse(disturb_rate == 0, "undisturbed", "disturbed"), 
         disturb_rate = factor(disturb_rate, levels = c("undisturbed", "disturbed"))) |> 
  group_by(disp_rate, kernel_exp, spat_heterogeneity, div_type, disturb_rate) |> 
  summarize(diversity = mean(diversity)) |> 
  ggplot(aes(x = disp_rate, y = kernel_exp, fill = diversity)) + 
  geom_tile(height = .53) + 
  scale_x_log10() +
  scale_y_log10() +
  scale_fill_viridis_c() +
  coord_fixed() +
  facet_grid(disturb_rate ~ spat_heterogeneity) +
  labs(x = "Emigration rate", y = "Dispersal kernel exponent", fill = "Diversity", title = "Gamma diversity",
       subtitle = "75 species in regional species pool")


heat_beta_spatial <- div_long |> 
  filter(div_type == "spatial beta\ndiversity") |> 
  mutate(disturb_rate = ifelse(disturb_rate == 0, "undisturbed", "disturbed"), 
         disturb_rate = factor(disturb_rate, levels = c("undisturbed", "disturbed"))) |> 
  group_by(disp_rate, kernel_exp, spat_heterogeneity, div_type, disturb_rate) |> 
  summarize(diversity = mean(diversity)) |> 
  ggplot(aes(x = disp_rate, y = kernel_exp, fill = diversity)) + 
  geom_tile(height = .53) + 
  scale_x_log10() +
  scale_y_log10() +
  scale_fill_viridis_c() +
  coord_fixed() +
  facet_grid(disturb_rate ~ spat_heterogeneity) +
  labs(x = "Emigration rate", y = "Dispersal kernel exponent", fill = "Diversity", title = "Spatial beta diversity",
       subtitle = "75 species in regional species pool")

heat_beta_temporal <- div_long |> 
  filter(div_type == "temporal beta\ndiversity") |> 
  mutate(disturb_rate = ifelse(disturb_rate == 0, "undisturbed", "disturbed"), 
         disturb_rate = factor(disturb_rate, levels = c("undisturbed", "disturbed"))) |> 
  group_by(disp_rate, kernel_exp, spat_heterogeneity, div_type, disturb_rate) |> 
  summarize(diversity = mean(diversity)) |> 
  ggplot(aes(x = disp_rate, y = kernel_exp, fill = diversity)) + 
  geom_tile(height = .53) + 
  scale_x_log10() +
  scale_y_log10() +
  scale_fill_viridis_c() +
  coord_fixed() +
  facet_grid(disturb_rate ~ spat_heterogeneity) +
  labs(x = "Emigration rate", y = "Dispersal kernel exponent", fill = "Diversity", title = "Temporal beta diversity",
       subtitle = "75 species in regional species pool")


heat_alpha_combined <- heat_alpha_gam_nodisturb + heat_alpha_gam_disturb +
  plot_layout(ncol = 1) 
ggsave(paste0("figures/",comp_scenario,"_diversity_alpha_gamma_heatmap_75species.pdf"), plot = heat_alpha_combined, width = 8, height = 8)
ggsave(paste0("figures/",comp_scenario,"_diversity_alpha_gamma_heatmap_75species.png"), plot = heat_alpha_combined, width = 8, height = 8, dpi = 500, bg = "white")

heat_beta_combined <- heat_beta_nodisturb + heat_beta_disturb +
  plot_layout(ncol = 1) 
ggsave(paste0("figures/",comp_scenario,"_diversity_beta_heatmap_75species.pdf"), plot = heat_beta_combined, width = 8, height = 8)
ggsave(paste0("figures/",comp_scenario,"_diversity_beta_heatmap_75species.png"), plot = heat_beta_combined, width = 8, height = 8, dpi = 500, bg = "white")


heat_alpha_gamma_combined <- heat_alpha + heat_gamma +
  plot_layout(ncol = 1) +
  plot_annotation(tag_levels = "A")
ggsave(paste0("figures/",comp_scenario,"_diversity_alpha_gamma_heatmap_75species.pdf"), plot = heat_alpha_gamma_combined, width = 8, height = 8)
ggsave(paste0("figures/",comp_scenario,"_diversity_alpha_gamma_heatmap_75species.png"), plot = heat_alpha_gamma_combined, width = 8, height = 8, dpi = 500, bg = "white")

heat_beta_combined <- heat_beta_spatial + heat_beta_temporal +
  plot_layout(ncol = 1)  +
  plot_annotation(tag_levels = "A")
ggsave(paste0("figures/",comp_scenario,"_diversity_beta_heatmap_75species.pdf"), plot = heat_beta_combined, width = 8, height = 8)
ggsave(paste0("figures/",comp_scenario,"_diversity_beta_heatmap_75species.png"), plot = heat_beta_combined, width = 8, height = 8, dpi = 500, bg = "white")






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
patches_over_time <- filter(patches_over_time, rep >5)
patches_over_time$species <- as.factor(patches_over_time$species)
patches_over_time$rep <- as.factor(patches_over_time$rep)
patches_over_time$patch <- as.factor(patches_over_time$patch)
patches_over_time$kernel_exp <- factor(signif(patches_over_time$kernel_exp, 1),
                                       levels = c(0, 1e-04, 3e-04, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1),
                                       labels = c("0", "0.0001", "0.0003", "0.001", "0.003", "0.01", "0.03", "0.1", "0.3", "1"))
patches_over_time <- patches_over_time |> 
  filter(spat_heterogeneity != 1) |> 
  mutate(kernel_exp = as.factor(kernel_exp),
         spat_heterogeneity = factor(spat_heterogeneity, 
                                     levels = c(0, 0.1, 1, 1000), 
                                     labels = c("temporal\nvariation", "spatial + temporal\nvariation", "spatial + temporal\nvariation (more spatial)", "spatial\nvariation")))


# analyze environmental filtering
# patches_over_time |> 
#   filter(comp == comp_scenario, extirp_prob == 0) |> 
#   group_by(emigration, kernel_exp, species, spat_heterogeneity) |> 
#   summarize(env_mismatch_mean = mean(env_mismatch),
#             env_mismatch_sd = sd(env_mismatch)) |> 
#   ggplot(aes(x = emigration, y = env_mismatch_mean, ymin = env_mismatch_mean - env_mismatch_sd,
#              ymax = env_mismatch_mean + env_mismatch_sd, color = species, shape = spat_heterogeneity)) + 
#   geom_point(alpha = 0.5) +
#   geom_errorbar(alpha = 0.5) +
#   geom_line() +
#   scale_x_log10() +
#   scale_fill_viridis_c() +
#   facet_grid(species~kernel_exp, scales = "free_y") +
#   labs(x = "Emigration rate", y = "Environmental mismatch (mean ± s.d.)")

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
fig_env_mismatch

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
fig_env_sinks

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
fig_env_costs


# analyze biotic filtering
fig_competition_effects <- patches_over_time |> 
  filter(comp == comp_scenario, extirp_prob == 0) |> 
  mutate(intra_comp = (delta_bio_intra)/abundance_mean,
         inter_comp = (delta_bio_inter)/abundance_mean) |> 
  group_by(emigration, kernel_exp, spat_heterogeneity) |> 
  summarize(intra_comp_mean = mean(intra_comp),
            inter_comp_mean = mean(inter_comp)) |> 
  ggplot(aes(x = emigration, color = kernel_exp)) + 
  geom_point(aes(y = intra_comp_mean), alpha = 0.5, shape = "square") +
  geom_line(aes(y = intra_comp_mean), linetype = "solid") +
  
  geom_point(aes(y = inter_comp_mean), alpha = 0.5) +
  geom_line(aes(y = inter_comp_mean), linetype = "dashed") +
  
  scale_x_log10() +
  scale_color_viridis_d(option = "B", end = .9) +
  facet_grid(. ~ spat_heterogeneity) +
  labs(x = "Emigration rate", 
       y = "Mean fitness effects of intra (solid) + \n inter (dashed) competition",
       color = "Dispersal kernel\n exponent")

fig_competition_intra <- patches_over_time |> 
  filter(comp == comp_scenario, extirp_prob == 0) |> 
  mutate(intra_comp = (delta_bio_intra)/abundance_mean,
         inter_comp = (delta_bio_inter)/abundance_mean) |> 
  group_by(emigration, kernel_exp, spat_heterogeneity) |> 
  summarize(intra_comp_mean = mean(intra_comp),
            inter_comp_mean = mean(inter_comp)) |> 
  ggplot(aes(x = emigration, color = kernel_exp)) + 
  geom_point(aes(y = intra_comp_mean), alpha = 0.5) +
  geom_line(aes(y = intra_comp_mean), linetype = "solid") +
  
  scale_x_log10() +
  scale_color_viridis_d(option = "B", end = .9) +
  facet_grid(. ~ spat_heterogeneity) +
  labs(x = "Emigration rate", 
       y = "Mean fitness effects of \n intraspecific competition",
       color = "Dispersal kernel\n exponent")

fig_competition_inter <- patches_over_time |> 
  filter(comp == comp_scenario, extirp_prob == 0) |> 
  mutate(intra_comp = (delta_bio_intra)/abundance_mean,
         inter_comp = (delta_bio_inter)/abundance_mean) |> 
  group_by(emigration, kernel_exp, spat_heterogeneity) |> 
  summarize(intra_comp_mean = mean(intra_comp),
            inter_comp_mean = mean(inter_comp)) |> 
  ggplot(aes(x = emigration, color = kernel_exp)) + 
  geom_point(aes(y = inter_comp_mean), alpha = 0.5) +
  geom_line(aes(y = inter_comp_mean), linetype = "solid") +
  
  scale_x_log10() +
  scale_color_viridis_d(option = "B", end = .9) +
  facet_grid(. ~ spat_heterogeneity) +
  labs(x = "Emigration rate", 
       y = "Mean fitness effects of \n interspecific competition",
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
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "ad"), se = F) +
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
  summarize(demo_extinctions = mean(stoch_extinct_demo),
            demo_extinctions_se = sd(stoch_extinct_demo)/5) |> 
  ggplot(aes(x = emigration, y = demo_extinctions, 
             ymin = demo_extinctions-demo_extinctions_se,
             ymax = demo_extinctions+demo_extinctions_se,
             color = kernel_exp)) + 
  geom_point(alpha = 0.5) +
  geom_errorbar(alpha = .5) +
  geom_line() +
  scale_x_log10() +
  scale_color_viridis_d(option = "B", end = .9) +
  facet_grid(. ~ spat_heterogeneity) +
  labs(x = "Emigration rate", 
       y = "Mean stochastic \nlocal extinctions",
       color = "Dispersal kernel\n exponent",
       subtitle = "20 species in regional species pool")


# combine abiotic and biotic in one graph
fig_env_bio_filter <- 
  fig_env_costs + theme(legend.position = "null") + labs(subtitle = "75 species in regional species pool") +
  fig_competition_intra + theme(legend.position = "null") + 
  fig_competition_inter + 
  plot_layout(nrow = 3, guides = "collect") +
  plot_annotation(tag_levels = "A")
ggsave(filename = paste0("figures/",comp_scenario,"_env_bio_filtering_nodisturb_75species.png"),
       plot = fig_env_bio_filter, width = 7, height = 7, dpi = 500, bg = "white")
ggsave(filename = paste0("figures/",comp_scenario,"_env_bio_filtering_nodisturb_75species.pdf"),
       plot = fig_env_bio_filter, width = 7, height = 7, dpi = 500)

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


fig_dispersal_fitness_effects <- patches_over_time |> 
  filter(comp == comp_scenario, extirp_prob == 0) |> 
  #mutate(demo_stoch = (W_max + delta_stoch_demo) / abundance_mean) |> 
  mutate(dispersal_effect = (delta_dispersal) / abundance_mean) |> 
  group_by(emigration, kernel_exp, spat_heterogeneity, extirp_prob) |> 
  summarize(dispersal_effect_mean = mean(dispersal_effect),
            dispersal_effect_sd = sd(dispersal_effect)/10) |> 
  ggplot(aes(x = emigration, y = dispersal_effect_mean,
             ymin = dispersal_effect_mean-dispersal_effect_sd,
             ymax = dispersal_effect_mean+dispersal_effect_sd,
             color = kernel_exp)) + 
  geom_point(alpha = 0.5) +
  geom_errorbar(alpha = .5) +
  geom_line() +
  scale_x_log10() +
  scale_color_viridis_d(option = "B", end = .9) +
  facet_grid(. ~ spat_heterogeneity) +
  labs(x = "Emigration rate", 
       y = "Mean fitness effects \n of dispersal",
       color = "Dispersal kernel\n exponent", 
       subtitle = "75 species in regional species pool")

fig_demog <- fig_demo_extinctions + fig_dispersal_fitness_effects +
  plot_layout(nrow = 2, guides = "collect") +
  plot_annotation(tag_levels = "A")
ggsave(filename = paste0("figures/",comp_scenario,"_stoch_dispersal_75species.png"),
       plot = fig_demog, width = 8, height = 6, dpi = 500, bg = "white")
ggsave(filename = paste0("figures/",comp_scenario,"_stoch_dispersal_75species.pdf"),
       plot = fig_demog, width = 8, height = 6)

