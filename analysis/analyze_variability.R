library(tidyverse)
library(here)


theme_set(theme_bw())

variability <- read_csv(here("data/variability_data.csv"))


var_long_stable <- variability %>% 
  na.omit() %>% 
  filter(disturb_rate < 0.1) %>% 
  filter(condition == "stable") %>% 
  pivot_longer(cols = c(CV_S_L, CV_C_L, CV_S_R, CV_C_R), 
               names_to = "variability", values_to = "CV") %>% 
  pivot_longer(cols = c(phi_S_L2R, phi_C_L2R, phi_S2C_L, phi_S2C_R), 
               names_to = "synchrony", values_to = "phi") %>% 
  mutate(variability = factor(variability, levels = c("CV_S_L", "CV_C_L", "CV_S_R", "CV_C_R")))

var_long_priority <- variability %>% 
  na.omit() %>% 
  filter(disturb_rate < 0.1) %>% 
  filter(condition == "priority") %>% 
  pivot_longer(cols = c(CV_S_L, CV_C_L, CV_S_R, CV_C_R), 
               names_to = "variability", values_to = "CV") %>% 
  pivot_longer(cols = c(phi_S_L2R, phi_C_L2R, phi_S2C_L, phi_S2C_R), 
               names_to = "synchrony", values_to = "phi") 


# Plot figures
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

var_long_stable %>% 
  filter(disturb_rate > 0) %>% 
  ggplot(aes(x = disp_rate, y = CV, color = as.factor(round(kernel_exp, 2)))) + 
  geom_point(alpha = 0.2) + 
  geom_smooth(se = FALSE) + 
  facet_grid(~variability, scales = "free_y") +
  scale_x_log10() +
  scale_color_viridis_d(option = "A") +
  labs(color = "Dispersal kernel \nexponent",
       x = "Emigration rate")

var_long_stable %>% 
  ggplot(aes(x = disp_rate, y = phi, color = as.factor(round(kernel_exp, 2)))) + 
  geom_point(alpha = 0.2) + 
  geom_smooth(se = FALSE) + 
  facet_grid(synchrony~round(disturb_rate,2), scales = "free_y") +
  scale_x_log10() +
  labs(color = "Dispersal kernel \nexponent",
       x = "Emigration rate")



var_long_priority %>% 
  ggplot(aes(x = disp_rate, y = CV, color = as.factor(round(kernel_exp, 2)))) + 
  geom_point(alpha = 0.2) + 
  geom_smooth(se = FALSE) + 
  facet_grid(variability~round(disturb_rate,2), scales = "free_y") +
  scale_x_log10() +
  labs(color = "Dispersal kernel \nexponent",
       x = "Emigration rate")

var_long_priority %>% 
  ggplot(aes(x = disp_rate, y = phi, color = as.factor(round(kernel_exp, 2)))) + 
  geom_point(alpha = 0.2) + 
  geom_smooth(se = FALSE) + 
  facet_grid(synchrony~round(disturb_rate,2), scales = "free_y") +
  scale_x_log10() +
  labs(color = "Dispersal kernel \nexponent",
       x = "Emigration rate")






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



variability %>% 
  filter(condition == "stable") %>% 
  filter(disturb_rate == 0) %>% 
  ggplot(aes(x = alpha_div, y = CV_C_L, color = as.factor(disp_rate))) + 
  geom_point(alpha = 0.3) + 
  geom_smooth(method = "lm", se=FALSE) +
  facet_wrap(~as.factor(round(kernel_exp,2)), scales = "free") + 
  scale_color_viridis_d() +
  theme_minimal() +
  labs(color = "Emigration rate", 
       x = "Mean alpha-diversity", 
       y = "Local community variability")

       