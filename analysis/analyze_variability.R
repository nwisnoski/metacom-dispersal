library(tidyverse)
library(here)

variability <- read_csv(here("data/variability_data.csv"))


variability %>% 
  na.omit() %>% 
  filter(disturb_rate < 0.1) %>% 
  filter(condition == "stable") %>% 
  pivot_longer(cols = c(CV_S_L, CV_C_L, CV_S_R, CV_C_R), 
               names_to = "variability", values_to = "CV") %>% 
  pivot_longer(cols = c(phi_S_L2R, phi_C_L2R), 
               names_to = "synchrony", values_to = "phi") %>% 
  ggplot(aes(x = disp_rate, y = CV, color = as.factor(kernel_exp))) + 
  geom_point(alpha = 0.2) + 
  geom_smooth(se = FALSE) + 
  facet_wrap(variability~disturb_rate, scales = "free_y") +
  scale_x_log10()

variability %>% 
  na.omit() %>% 
  filter(disturb_rate < 0.1) %>% 
  filter(condition == "stable") %>% 
  pivot_longer(cols = c(CV_S_L, CV_C_L, CV_S_R, CV_C_R), 
               names_to = "variability", values_to = "CV") %>% 
  pivot_longer(cols = c(phi_S_L2R, phi_C_L2R), 
               names_to = "synchrony", values_to = "phi") %>% 
  ggplot(aes(x = disp_rate, y = phi, color = as.factor(kernel_exp))) + 
  geom_point(alpha = 0.2) + 
  geom_smooth(se = FALSE) + 
  facet_wrap(synchrony~disturb_rate, scales = "free_y") +
  scale_x_log10()



variability %>% 
  na.omit() %>% 
  filter(disturb_rate < 0.1) %>% 
  filter(condition == "priority") %>% 
  pivot_longer(cols = c(CV_S_L, CV_C_L, CV_S_R, CV_C_R), 
               names_to = "variability", values_to = "CV") %>% 
  pivot_longer(cols = c(phi_S_L2R, phi_C_L2R), 
               names_to = "synchrony", values_to = "phi") %>% 
  ggplot(aes(x = disp_rate, y = CV, color = as.factor(kernel_exp))) + 
  geom_point(alpha = 0.2) + 
  geom_smooth(se = FALSE) + 
  facet_wrap(variability~disturb_rate, scales = "free_y") +
  scale_x_log10()

variability %>% 
  na.omit() %>% 
  filter(disturb_rate < 0.1) %>% 
  filter(condition == "priority") %>% 
  pivot_longer(cols = c(CV_S_L, CV_C_L, CV_S_R, CV_C_R), 
               names_to = "variability", values_to = "CV") %>% 
  pivot_longer(cols = c(phi_S_L2R, phi_C_L2R), 
               names_to = "synchrony", values_to = "phi") %>% 
  ggplot(aes(x = disp_rate, y = phi, color = as.factor(kernel_exp))) + 
  geom_point(alpha = 0.2) + 
  geom_smooth(se = FALSE) + 
  facet_wrap(synchrony~disturb_rate, scales = "free_y") +
  scale_x_log10()
