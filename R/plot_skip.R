#plot skip spawning inputs
library(patchwork)
library(tidyverse)
source(here::here("R", "functions.R"))
n_dat <- readRDS(here::here("data", "n_dat.RDS"))
d_dat <- readRDS(here::here("data", "d_dat.RDS"))
data <- readRDS(here::here("data", "goa_nork_dat.RDS"))
rpt <- readRDS(here::here("output", "rpt.RDS"))

skip_levels = c(0, 0.1, 0.2, 0.3)

# northerns
smin = 3
smax = 50

expand.grid(age = 1:50,
            skip = skip_levels) %>%
  mutate(dome = mapply(flexi_curve, age, skip, smin, smax, 'dome'),
         "skewed dome" = mapply(flexi_curve, age, skip, smin, smax, 'skewed_dome', skew=1.5),
         increasing = mapply(flexi_curve, age, skip, smin, smax, 'increase'),
         decreasing = mapply(flexi_curve, age, skip, smin, smax, 'decrease'),
         constant = ifelse(age %in% smin:smax, skip, 0),
         # shift =  mapply(flexi_curve, age, skip, smin, smax, 'decrease'),
         "concave" = mapply(flexi_curve, age, skip, smin, smax, 'inverse_dome', skew = 0.2))  -> p_skip

p4 = p_skip %>%
  pivot_longer(-c(age, skip)) %>% 
  left_join(data.frame(age = 2:51, maa = filter(n_dat, age %in% 2:51) %>% pull(biological))) %>% 
  mutate("Skip rate"= factor(skip),
         value = (1 - value) * maa) %>%
  ggplot(aes(age, value, color = `Skip rate`)) + 
  geom_line() +
  facet_wrap(~name) +
  theme(legend.position = c(0.8, 0.2)) +
  scico::scale_color_scico_d(palette = 'roma') +
  xlab("Age") +
  ylab("Proportion mature") +
  afscassess::theme_report() +
  # tickr::scale_x_tickr(data = p_skip, var = age, by = 20) +
  geom_line(data = data.frame(age = 2:51, value = filter(n_dat, age %in% 2:51) %>% pull(functional), Skip = "observed"), color = "gray") 

p3 = p_skip %>%
  pivot_longer(-c(age, skip)) %>% 
  mutate("Skip rate"= factor(skip)) %>% 
  ggplot(aes(age, value, color = `Skip rate`)) + 
  geom_line() +
  facet_wrap(~name) +
  theme(legend.position = c(0.8, 0.2)) +
  scico::scale_color_scico_d(palette = 'roma') +
  xlab("Age") +
  ylab("Proportion skip spawning") +
  # tickr::scale_x_tickr(data = p_skip, var = age, by = 10) +
  afscassess::theme_report()

# dusky
smin = 4
smax = 33

expand.grid(age = 1:33,
            skip = skip_levels) %>%
  mutate(dome = mapply(flexi_curve, age, skip, smin, smax, 'dome'),
         "skewed dome" = mapply(flexi_curve, age, skip, smin, smax, 'skewed_dome', skew=1.5),
         increasing = mapply(flexi_curve, age, skip, smin, smax, 'increase'),
         decreasing = mapply(flexi_curve, age, skip, smin, smax, 'decrease'),
         constant = ifelse(age %in% smin:smax, skip, 0),
         # shift =  mapply(flexi_curve, age, skip, smin, smax, 'decrease'),
         "concave" = mapply(flexi_curve, age, skip, smin, smax, 'inverse_dome', skew = 0.2)) -> p_skip_d

p1 = p_skip_d %>%
  pivot_longer(-c(age, skip)) %>% 
  left_join(data.frame(age = 4:33, maa = filter(d_dat, age %in% 4:33) %>% pull(biological))) %>% 
  mutate("Skip rate"= factor(skip),
         value = (1 - value) * maa) %>%
  ggplot(aes(age, value, color = `Skip rate`)) + 
  geom_line() +
  facet_wrap(~name) +
  theme(legend.position = c(0.8, 0.2)) +
  scico::scale_color_scico_d(palette = 'roma') +
  xlab("Age") +
  ylab("Proportion mature") +
  afscassess::theme_report() +
  # tickr::scale_x_tickr(data = p_skip_d, var = age, by = 20) +
  geom_line(data = d_dat %>% 
    				mutate(Skip = "observed") %>% 
    				filter(age %in% 4:33), aes(y = functional), color = "gray")

p2 = p_skip_d %>%
  filter(age %in% 4:33) %>% 
  pivot_longer(-c(age, skip)) %>% 
  mutate("Skip rate"= factor(skip)) %>% 
  ggplot(aes(age, value, color = `Skip rate`)) + 
  geom_line() +
  facet_wrap(~name) +
  theme(legend.position = c(0.8, 0.2)) +
  scico::scale_color_scico_d(palette = 'roma') +
  xlab("Age") +
  ylab("Proportion skip spawning") +
  # tickr::scale_x_tickr(data = p_skip_d, var = age, by = 10) +
  afscassess::theme_report()

(p2 / p3) 
(p1 / p4)
  
 
ggsave(here::here("figs", "skip.png"), width = 6.5, height = 6.5, dpi = 200)
  