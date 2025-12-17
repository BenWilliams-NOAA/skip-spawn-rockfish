# updating northern rockfish maturity 
# ben.williams@noaa.gov
# 2025-12

# notes ----
# northern rockfish have historically been assessed using biological maturity with a weighted logistic
# regression - in order to keep the origin near 0 as there is limited info for young fish
# two projects were used to contribute to maturity 
# Lunsford 1999 & Chilton 2003 <----- check dates...
# Conrath 2017 examined skip spawning in northern rockfish (by length)
# age data are also available for these samples
# two analyses are performed.
# 1) update biological maturity with these samples
# 2) estimate skip spawning at age
# load ----
library(readxl)
library(tidytable)
library(ggplot2)
library(mgcv)

# current maturity estimate
og = read_xlsx(here::here("data", "Conrath.xlsx"), sheet = "orig")
og %>%
  tidytable::uncount(n_lunsford, .id = 'id') %>% 
  tidytable::mutate(mat = ifelse(mat_lunsford>0 & mat_lunsford>=id, 1, 0), .by =age) %>% 
  tidytable::select(age, mat, wt = wt_lunsford) -> lunds

og %>%
  tidytable::uncount(n_chilton, .id = 'id') %>% 
  tidytable::mutate(mat = ifelse(mat_chilton>0 & mat_chilton>=id, 1, 0), .by =age) %>% 
  tidytable::select(age, mat, wt = wt_chilton) -> chilt

bind_rows(lunds, chilt) -> og

# skip spawning data
# step 1. examine northern rockfish skip spawning
read_xlsx(here::here("data", "Conrath.xlsx"), sheet = "Northern RF") %>% 
  select(length = FL, weight = WT, skip = SS, age) %>% 
  mutate(length = length /10,
         mat = 1,
         wt = ifelse(age<=6, 50, 1)) -> update

update %>% 
  select(age, mat, wt) %>% 
  bind_rows(og) -> dat


m0 = glm(mat ~ age, data = og, family = 'binomial', weights =wt) 
m1 = glm(mat ~ age, data = dat, family = 'binomial', weights =wt) 
summary(m0)
summary(m1)


data.frame(age = 2:51) %>% 
  mutate(fit = predict.glm(m0, ., type = "response"),
         fit1 = predict.glm(m1, ., type = "response")) %>% 
  pivot_longer(-age) %>% 
  ggplot(aes(age, value, color = name)) + 
  geom_line() 

# fit1 is the biological maturity that is used 
# in the assessment (updated in 2024) & for this study

# skip spawning 
read_xlsx(here::here("data", "Conrath.xlsx"), sheet = "Northern RF") %>% 
  select(length = FL, weight = WT, skip = SS, age) %>% 
  mutate(length = length /10,
         mat = ifelse(skip==1, 0, 1),
         wt = ifelse(age<=6, 50, 1)) %>% 
  arrange(age) -> nr_skip

ggplot(nr_skip, aes(age, skip)) + 
  geom_point() + 
  stat_smooth()

# explore skip spawning by age
m2 <- gam(skip ~ s(age), data = nr_skip, family="binomial")
summary(m2)
preds <- predict.gam(m2, data.frame(age = 2:51), type = "response", se.fit = TRUE)

# have no skip spawning information below age-12, so fixed values at age-12 level
# otherwise will have a large "bump" in the maturity curve
data.frame(age = 2:51) %>% 
  mutate(skip = c(as.numeric(preds$fit)),
         skip = ifelse(age<12, 0.32921833 , skip),
         biological = predict.glm(m1, ., type = "response"),
         functional = (1-skip) * predict.glm(m1, ., type = "response")) -> dat1


dat1 %>% pull(functional) -> skip
saveRDS(skip, here::here("data", "skip.RDS"))

dat1 %>% 
  pivot_longer(-age, names_to = "Model") %>% 
  ggplot(aes(age, value, color = Model)) + 
  annotate("rect", xmin = 0, xmax = 11.5, ymin = 0, ymax = 1, 
           fill = "gray", alpha = 0.3) +
  geom_line() +
  scico::scale_color_scico_d(palette = 'roma') +
  afscassess::theme_report(base_size = 16) +
  theme(legend.position = c(0.8, 0.4)) +
  xlab("Age") +
  ylab("Proportion mature")
  


