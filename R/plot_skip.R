#plot skip spawning inputs
library(patchwork)
library(tidyverse)
source(here::here("R", "functions.R"))
skip <- readRDS(here::here("data", "skip.RDS"))
data <- readRDS(here::here("data", "goa_nork_dat.RDS"))
rpt <- readRDS(here::here("output", "rpt.RDS"))

max_age = length(rpt$waa)
skip_levels = c(0, 0.1, 0.2, 0.3)
smin = 3
smax = 50

expand.grid(age = 1:50,
            skip = skip_levels) %>%
  mutate(dome = mapply(flexi_curve, age, skip, smin, smax, 'dome'),
         skewed_dome = mapply(flexi_curve, age, skip, smin, smax, 'skewed_dome', skew=1.5),
         increasing = mapply(flexi_curve, age, skip, smin, smax, 'increase'),
         decreasing = mapply(flexi_curve, age, skip, smin, smax, 'decrease'),
         constant = ifelse(age %in% smin:smax, skip, 0),
         # shift =  mapply(flexi_curve, age, skip, smin, smax, 'decrease'),
         inverse_dome = mapply(flexi_curve, age, skip, smin, smax, 'inverse_dome', skew = 0.2)) -> p_skip

p_skip %>% 
  filter(skip ==0.2) %>% 
  select(age, inverse_dome) %>% 
  mutate(inverse_dome = (1-inverse_dome) * data$maa) %>% 
  bind_cols(data.frame(skip = skip)) %>% 
  pivot_longer(-age) %>% 
  ggplot(aes(age, value, color = name)) + 
  geom_line()

p2 = p_skip %>%
  pivot_longer(-c(age, skip)) %>% 
  left_join(data.frame(age = 2:51, maa = data$maa)) %>% 
  mutate(Skip = factor(skip),
         value = (1 - value) * maa) %>%
  ggplot(aes(age, value, color = Skip)) + 
  geom_line() +
  facet_wrap(~name) +
  theme(legend.position = c(0.8, 0.2)) +
  scico::scale_color_scico_d(palette = 'roma') +
  xlab("Age") +
  ylab("Proportion mature") +
  afscassess::theme_report() +
  geom_line(data = data.frame(age = 1:50, value = skip, Skip = "observed"), color = "gray")

p1 = p_skip %>%
  pivot_longer(-c(age, skip)) %>% 
  mutate(Skip = factor(skip)) %>% 
  ggplot(aes(age, value, color = Skip)) + 
  geom_line() +
  facet_wrap(~name) +
  theme(legend.position = c(0.8, 0.2)) +
  scico::scale_color_scico_d(palette = 'roma') +
  xlab("Age") +
  ylab("Proportion skip spawning") +
  afscassess::theme_report()

p1 + p2 +
  plot_layout(guides = 'collect') 
 
p_skip %>%
  pivot_longer(-c(age, skip)) %>% 
  left_join(data.frame(age = 2:51, maa = data$maa)) %>% 
  mutate(Skip = factor(skip),
         value = (1 - value) * maa) %>%
  filter(name=="inverse_dome") %>% 
  ggplot(aes(age, value, color = Skip)) + 
  geom_line(alpha = 0.5) +
  facet_wrap(~name) +
  theme(legend.position = c(0.8, 0.2)) +
  scico::scale_color_scico_d(palette = 'roma') +
  xlab("Age") +
  ylab("Proportion mature") +
  afscassess::theme_report() +
  geom_line(data = data.frame(age = 1:50, value = skip, Skip = "observed"), color = "black")

ggsave(here::here("figs", "skip2.png"), width = 6.5, height = 6.5, dpi = 200)
  
####################################################
results_df <- extract_results(sim1, start_year = 2024)
results_df2 <- extract_results(avg_sim_inverse_dome_0.3 , start_year = 2024)
  
d1 <- results_df %>%
  select(iteration, year, scenario, B0, B40) %>%
  tidyr::pivot_wider(names_from = scenario, values_from = c(B0, B40)) %>%
  select(iteration, year, B0_i, B40_i)

d2 <- results_df2 %>%
  select(iteration, year, scenario, B0, B40) %>%
  tidyr::pivot_wider(names_from = scenario, values_from = c(B0, B40)) %>%
  select(iteration, year, B0_ii, B40_ii)

comparison = left_join(d1, d2) %>% 
  mutate(
    # How much larger is the Inverse Dome B0 compared to the Truth?
    B0_Inflation = B0_ii / B0_i,
    # How much higher is the Management Target?
    Target_Inflation = B40_ii / B40_i
  )

# 2. Check the mean inflation
mean(comparison$B0_Inflation, na.rm=TRUE) 
# If this is > 1.0 (e.g., 1.2 or 1.5), it confirms the Inverse Dome 
# is hallucinating extra productivity and setting the bar too high.


comparison %>% 
  # filter(iteration == 3) %>% 
  ggplot(aes(year, B40_ii/B0_i, group = iteration)) + 
  geom_line() +
  geom_line(aes(y = B0_ii/B40_ii), color = 4) 
  geom_line(aes(y = B0_ii), color = 4) +
  geom_line(aes(y = B0_i))
  

# 3. Visualize the "High Bar" effect
results_df %>% 
  filter(iteration == 3) %>% # Pick one iteration to see the trajectory
  ggplot(aes(x = year)) +
  # Plot the Estimated B0 for both scenarios
  geom_line(aes(y = B0, color = scenario, linetype = "Est. Virgin Biomass (B0)")) +
  # Plot the Target B40
  geom_line(aes(y = B40, color = scenario, linetype = "Target (B40)")) +
  labs(title = "Why Inverse Dome is Conservative",
       subtitle = "Inverse Dome (ii) calculates a much higher Baseline (B0) & Target (B40)",
       y = "Biomass") +
  afscassess::theme_report() +
  expand_limits(y=c(0, 90000))
