# updating northern and dusky rockfish maturity
# ben.williams@noaa.gov
# 2025-12

# notes ----
# these rockfish have historically been assessed using biological maturity with a weighted logistic
# regression - in order to keep the origin near 0 as there is limited info for young fish
# two projects were used to contribute to maturity 
# Lunsford 1999 & Chilton 2003 <----- check dates...
# Conrath 2017 examined skip spawning in rockfish (by length)
# age data are also available for these samples
# two analyses are performed.
# 1) update biological maturity with these samples
# 2) estimate skip spawning at age

# load ----
library(readxl)
library(tidytable)
library(ggplot2)
library(mgcv)
theme_set(afscassess::theme_report())

# data ----
# dusky
og = read_xlsx(here::here("data", "Conrath.xlsx"), sheet = "orig_dusk")
og %>%
  tidytable::uncount(n_lunsford, .id = 'id') %>% 
  tidytable::mutate(mat = ifelse(mat_lunsford>0 & mat_lunsford>=id, 1, 0), .by =age) %>% 
  tidytable::select(age, mat, wt) -> lunds

og %>%
  tidytable::uncount(n_chilton, .id = 'id') %>% 
  tidytable::mutate(mat = ifelse(mat_chilton>0 & mat_chilton>=id, 1, 0), .by =age) %>% 
  tidytable::select(age, mat, wt) -> chilt

bind_rows(lunds, chilt) -> dusk

# skip spawning data 
# note <= age-5 has greater weight
read_xlsx(here::here("data", "Conrath.xlsx"), sheet = "Dusky RF") %>% 
  select(length = FL, weight = WT, skip = SS, age) %>% 
  mutate(length = length /10,
         mat = 1,
         wt = ifelse(age<=5, 50, 1)) -> d_skip

d_skip %>% 
  select(age, mat, wt) %>% 
  bind_rows(dusk) -> d_dat

# northern
og = read_xlsx(here::here("data", "Conrath.xlsx"), sheet = "orig")
og %>%
  tidytable::uncount(n_lunsford, .id = 'id') %>% 
  tidytable::mutate(mat = ifelse(mat_lunsford>0 & mat_lunsford>=id, 1, 0), .by =age) %>% 
  tidytable::select(age, mat, wt = wt_lunsford) -> lunds

og %>%
  tidytable::uncount(n_chilton, .id = 'id') %>% 
  tidytable::mutate(mat = ifelse(mat_chilton>0 & mat_chilton>=id, 1, 0), .by =age) %>% 
  tidytable::select(age, mat, wt = wt_chilton) -> chilt

bind_rows(lunds, chilt) -> nork

# skip spawning data 
# note <= age-6 has greater weight
read_xlsx(here::here("data", "Conrath.xlsx"), sheet = "Northern RF") %>% 
  select(length = FL, weight = WT, skip = SS, age) %>% 
  mutate(length = length /10,
         mat = 1,
         wt = ifelse(age<=6, 50, 1)) -> n_skip

n_skip %>% 
  select(age, mat, wt) %>% 
  bind_rows(nork) -> n_dat


# basic comparisons
# dusky
m0 = glm(mat ~ age, data = dusk, family = 'binomial', weights =wt) 
m1 = glm(mat ~ age, data = d_dat, family = 'binomial', weights =wt) 
summary(m0)
summary(m1)


data.frame(age = 4:33) %>% 
  mutate(fit = predict.glm(m0, ., type = "response"),
         fit1 = predict.glm(m1, ., type = "response")) %>% 
  pivot_longer(-age) %>% 
  ggplot(aes(age, value, color = name)) + 
  geom_line() 

# old pars - dusky
coef(m0)
# (Intercept)         age 
#  -6.8159296   0.6593501 

# -coef(m0)[1] / coef(m0)[2]
# A50 = 10.33734 

# new pars - dusky
coef(m1)
# (Intercept)         age 
#  -7.4957748   0.7626374 

# -coef(m1)[1] / coef(m1)[2]
# A50 =  9.828754

# m0 is the biological maturity that is used in the assessment (should be updated to m1)
# m1 is used in this study

# northern
m2 = glm(mat ~ age, data = nork, family = 'binomial', weights =wt) 
m3 = glm(mat ~ age, data = n_dat, family = 'binomial', weights =wt) 
summary(m2)
summary(m3)


data.frame(age = 2:51) %>% 
  mutate(fit = predict.glm(m2, ., type = "response"),
         fit1 = predict.glm(m3, ., type = "response")) %>% 
  pivot_longer(-age) %>% 
  ggplot(aes(age, value, color = name)) + 
  geom_line() 

# old pars - dusky
coef(m2)
# (Intercept)         age 
#  -6.8849136   0.6488762 

# -coef(m2)[1] / coef(m2)[2]
# A50 = 10.61052

# new pars - dusky
coef(m3)
# (Intercept)         age 
#  -7.1627793   0.6915622  

# -coef(m3)[1] / coef(m3)[2]
# A50 =   10.35739 

# m3 is the biological maturity that is used in the assessment & this study

# skip spawning ---
# dusky
ggplot(d_skip, aes(age, skip)) + 
  geom_point() + 
  stat_smooth(method = "gam", method.args = list(family = "binomial"))

# explore skip spawning by age
m1s <- gam(skip ~ s(age), data = d_skip, family="binomial")
summary(m1s)
d_preds <- predict.gam(m1s, data.frame(age = 1:38), type = "response", se.fit = TRUE)

# have no skip spawning information below age-10 so fixed values at age-10 level
# otherwise will have a large "bump" in the maturity curve
min(d_skip$age)
data.frame(age = 1:38, 
  skip = d_preds$fit) %>% 
  mutate(skip = ifelse(age < min(d_skip$age), skip[age==min(d_skip$age)] , skip),
				 biological = predict.glm(m1, ., type = "response"),
         functional = (1-skip) * biological) -> d_dat

ggplot(d_dat, aes(age, skip)) + 
  geom_point() 

# skip spawning ---
# northern
ggplot(n_skip, aes(age, skip)) + 
  geom_point() + 
  stat_smooth(method = "gam", method.args = list(family = "binomial"))

# explore skip spawning by age
m3s <- gam(skip ~ s(age), data = n_skip, family="binomial")
summary(m3s)
n_preds <- predict.gam(m3s, data.frame(age = 1:51), type = "response", se.fit = TRUE)

# have no skip spawning information below age-12 so fixed values at age-12 level
# otherwise will have a large "bump" in the maturity curve
min(n_skip$age)
data.frame(age = 1:51, 
  skip = n_preds$fit) %>% 
  mutate(skip = ifelse(age < min(n_skip$age), skip[age== min(n_skip$age)] , skip),
				 biological = predict.glm(m3, ., type = "response"),
         functional = (1-skip) * biological) -> n_dat

ggplot(n_dat, aes(age, skip)) + 
  geom_point() 

max(n_skip$age)
max(d_skip$age)

p = d_dat %>% 
  mutate(id = "dusky rockfish") %>% 
  bind_rows(n_dat %>% 
    						mutate(id = "northern rockfish")) %>% 
  pivot_longer(-c(age, id), names_to = "Model") %>% 
  mutate(Model = case_when(Model=="skip" ~ "skip spawning",
                           TRUE ~ Model)) %>% 
  ggplot(aes(age, value, color = Model)) + 
  annotate("rect", xmin = 0, xmax = 11.5, ymin = 0, ymax = 1, 
           fill = "gray", alpha = 0.3) +
  geom_line() +
  scico::scale_color_scico_d(palette = 'grayC', end = 0.8) +
  facet_wrap(~id) +
  # afscassess::theme_report(base_size = 16) +
  theme(legend.position = c(0.8, 0.4)) +
  xlab("Age") +
  ylab("Proportion") +
  tickr::scale_x_tickr(data=dat1, var = age, var_min=0)

ggsave(here::here("figs", "maturity.png"), p, grDevices::png, width = 6.5, height = 5.5, dpi = 200)


d_dat %>% 
  filter(age %in% 4:33) %>% 
  pull(skip) -> dusk_skip

saveRDS(dusk_skip, here::here("data", "dusk_skip.RDS"))

n_dat %>% 
  filter(age %in% 2:51) %>% 
  pull(skip) -> nork_skip

saveRDS(nork_skip, here::here("data", "nork_skip.RDS"))
