library(tidyverse)
library(here)
library(patchwork)
library(data.table)

# Make figures for the simulation output from simulations/metacom_dispersal-kernel.R"
theme_set(theme_bw())

subfolder <- "sim_output/2026-04-29/"
comp_scenario <- "stable"

# load and combine reps
diversity <- data.table()
i = 1
for(file in list.files(path = subfolder, pattern = "summary\\.csv$")){
  this_run <- fread(paste(subfolder, file, sep = "/"))
  this_run$rep <- i
  diversity <- bind_rows(diversity, this_run)
  i = i + 1
} 

div_df <- diversity |> 
  mutate(beta_div = replace_na(beta_div, 0),
         kernel_exp = round(kernel_exp, 4)) 

kernel_exps <- sort(unique(div_df$kernel_exp))
disp_rates <- sort(unique(div_df$disp_rate))
disturb_rates <- sort(unique(div_df$disturb_rate))

# Diversity patterns
div_long <- div_df |> 
  filter(condition == comp_scenario) |>  
  pivot_longer(cols = c(alpha_div, gamma_div, beta_spatial, beta_temporal),
               names_to = "div_type", values_to = "diversity") |> 
  mutate(
    div_type = factor(div_type, 
                      levels = c("alpha_div", "gamma_div", "beta_spatial", "beta_temporal"),
                      labels = c("alpha diversity", "gamma diversity", str_wrap("spatial beta diversity", width = 15), str_wrap("temporal beta diversity", width = 15)))) |> 
  mutate(spat_heterogeneity = factor(spat_heterogeneity, 
                                levels = c(0, 0.1, 1000), 
                                labels = c("temporal\nvariation", "spatial + temporal\nvariation", "spatial\nvariation")))


# Figs 2-3: make diversity heatmaps
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
  labs(x = "Emigration rate", y = "Dispersal kernel exponent", fill = "Diversity", title = "Alpha diversity")

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
  labs(x = "Emigration rate", y = "Dispersal kernel exponent", fill = "Diversity", title = "Gamma diversity")


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
  labs(x = "Emigration rate", y = "Dispersal kernel exponent", fill = "Diversity", title = "Spatial beta diversity")

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
  labs(x = "Emigration rate", y = "Dispersal kernel exponent", fill = "Diversity", title = "Temporal beta diversity")


# Fig 2
heat_alpha_gamma_combined <- heat_alpha + heat_gamma +
  plot_layout(ncol = 1) +
  plot_annotation(tag_levels = "A")
ggsave(paste0("figures/Fig2.pdf"), plot = heat_alpha_gamma_combined, width = 8, height = 8)
ggsave(paste0("figures/Fig2.tif"), plot = heat_alpha_gamma_combined, width = 8, height = 8, dpi = 700, bg = "white")

# Fig 3
heat_beta_combined <- heat_beta_spatial + heat_beta_temporal +
  plot_layout(ncol = 1)  +
  plot_annotation(tag_levels = "A")
ggsave(paste0("figures/Fig3.pdf"), plot = heat_beta_combined, width = 8, height = 8)
ggsave(paste0("figures/Fig3.tif"), plot = heat_beta_combined, width = 8, height = 8, dpi = 700, bg = "white")


# Figs 4-5 Partition fitness effects
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
nrep <- i-1 # how many replicates are here (for se calculation)

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
                                     levels = c(0, 0.1, 1000), 
                                     labels = c("temporal\nvariation", "spatial + temporal\nvariation", "spatial\nvariation")))


# analyze environmental filtering
fig_env_costs <- patches_over_time |> 
  filter(comp == comp_scenario, extirp_prob == 0) |> 
  mutate(tot_comp = (delta_env) / abundance_mean) |> 
  group_by(emigration, kernel_exp, spat_heterogeneity) |> 
  summarize(comp_mean = mean(tot_comp),
            comp_se = sd(tot_comp)/nrep) |> 
  ggplot(aes(x = emigration, y = comp_mean, 
             ymin = comp_mean-comp_se,
             ymax = comp_mean+comp_se,
             color = kernel_exp)) + 
  geom_point(alpha = 0.5) +
  geom_errorbar(alpha = 0.5) +
  geom_line() +
  scale_x_log10() +
  scale_color_viridis_d(option = "B", end = .9) +
  facet_grid(. ~ spat_heterogeneity) +
  labs(x = "Emigration rate", 
       y = "Mean fitness effects \n of env. filtering",
       color = "Dispersal kernel\n exponent")

# analyze biotic filtering
fig_competition_intra <- patches_over_time |> 
  filter(comp == comp_scenario, extirp_prob == 0) |> 
  mutate(intra_comp = (delta_bio_intra)/abundance_mean,
         inter_comp = (delta_bio_inter)/abundance_mean) |> 
  group_by(emigration, kernel_exp, spat_heterogeneity) |> 
  summarize(intra_comp_mean = mean(intra_comp),
            inter_comp_mean = mean(inter_comp),
            intra_comp_se = sd(intra_comp)/nrep) |> 
  ggplot(aes(x = emigration, color = kernel_exp)) + 
  geom_point(aes(y = intra_comp_mean), alpha = 0.5) +
  geom_errorbar(aes(ymin = intra_comp_mean-intra_comp_se,
                    ymax = intra_comp_mean+intra_comp_se), alpha = 0.5) +
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
            inter_comp_mean = mean(inter_comp),
            inter_comp_se = sd(inter_comp)/nrep) |> 
  ggplot(aes(x = emigration, color = kernel_exp)) + 
  geom_point(aes(y = inter_comp_mean), alpha = 0.5) +
  geom_errorbar(aes(ymin = inter_comp_mean-inter_comp_se,
                    ymax = inter_comp_mean+inter_comp_se)) +
  geom_line(aes(y = inter_comp_mean), linetype = "solid") +
  
  scale_x_log10() +
  scale_color_viridis_d(option = "B", end = .9) +
  facet_grid(. ~ spat_heterogeneity) +
  labs(x = "Emigration rate", 
       y = "Mean fitness effects of \n interspecific competition",
       color = "Dispersal kernel\n exponent")

# Fig 4: combine abiotic and biotic in one graph
fig_env_bio_filter <- 
  fig_env_costs + theme(legend.position = "null") + 
  fig_competition_intra + theme(legend.position = "null") + 
  fig_competition_inter + 
  plot_layout(nrow = 3, guides = "collect") +
  plot_annotation(tag_levels = "A")
ggsave(filename = paste0("figures/Fig4.png"),
       plot = fig_env_bio_filter, width = 7, height = 7, dpi = 700, bg = "white")
ggsave(filename = paste0("figures/Fig4.pdf"),
       plot = fig_env_bio_filter, width = 7, height = 7)



# Fig 5: now quantify stochastic and dispersal effects
fig_demo_extinctions <- patches_over_time |> 
  filter(comp == comp_scenario, extirp_prob == 0) |> 
  group_by(emigration, kernel_exp, spat_heterogeneity) |> 
  summarize(demo_extinctions = mean(stoch_extinct_demo),
            demo_extinctions_se = sd(stoch_extinct_demo)/nrep) |> 
  ggplot(aes(x = emigration, y = demo_extinctions, 
             ymin = demo_extinctions-demo_extinctions_se,
             ymax = demo_extinctions+demo_extinctions_se,
             color = kernel_exp)) + 
  geom_point(alpha = 0.5) +
  geom_errorbar(alpha = 0.5) +
  geom_line() +
  scale_x_log10() +
  scale_color_viridis_d(option = "B", end = .9) +
  facet_grid(. ~ spat_heterogeneity) +
  labs(x = "Emigration rate", 
       y = "Mean stochastic \nlocal extinctions",
       color = "Dispersal kernel\n exponent")


fig_dispersal_fitness_effects <- patches_over_time |> 
  filter(comp == comp_scenario, extirp_prob == 0) |> 
  mutate(dispersal_effect = (delta_dispersal) / abundance_mean) |> 
  group_by(emigration, kernel_exp, spat_heterogeneity, extirp_prob) |> 
  summarize(dispersal_effect_mean = mean(dispersal_effect),
            dispersal_effect_se = sd(dispersal_effect)/nrep) |> 
  ggplot(aes(x = emigration, y = dispersal_effect_mean,
             ymin = dispersal_effect_mean-dispersal_effect_se,
             ymax = dispersal_effect_mean+dispersal_effect_se,
             color = kernel_exp)) + 
  geom_point(alpha = 0.5) +
  geom_errorbar(alpha = .5) +
  geom_line() +
  scale_x_log10() +
  scale_color_viridis_d(option = "B", end = .9) +
  facet_grid(. ~ spat_heterogeneity) +
  labs(x = "Emigration rate", 
       y = "Mean fitness effects \n of dispersal",
       color = "Dispersal kernel\n exponent")


fig_demog <- fig_demo_extinctions + fig_dispersal_fitness_effects +
  plot_layout(nrow = 2, guides = "collect") +
  plot_annotation(tag_levels = "A")
ggsave(filename = paste0("figures/Fig5.png"),
       plot = fig_demog, width = 8, height = 6, dpi = 700, bg = "white")
ggsave(filename = paste0("figures/Fig5.pdf"),
       plot = fig_demog, width = 8, height = 6)

