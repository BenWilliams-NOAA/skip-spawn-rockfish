# notes ----
# MSE self test - check to replicate the last assessment report
# making sure the 'project_step' function is correct
# ben.williams@noaa.gov
# 2025-12

# load ----
source(here::here("R", "utils.R"))
source(here::here("R", "functions.R"))
source(here::here("R", "em.R"))
# data ----
data <- readRDS(here::here("data", "goa_nork_dat.RDS"))
pars <- readRDS(here::here("data", "goa_nork_pars.RDS"))

# base maturity
# duplicate the biological maturity for this test case, also add necessary parameters
data$wt_mature_f = data$wt_mature
pars$log_F50_f = pars$log_F40_f = pars$log_F35_f = 0

# run base model
obj <- try(RTMB::MakeADFun(cmb(f1, data),
                           parameters = pars,
                           map = list(sigmaR = factor(NA))))
fit = nlminb(start = obj$par,
             objective = obj$fn,
             gradient = obj$gr,
             control = list(iter.max = 100000,
                            eval.max = 20000))
rpt <- obj$report(obj$env$last.par.best)

# check mechanics of project_step
test_project <- self_test(rpt, project_step)

test_project %>% 
  pivot_longer(c(report_SSB, new_SSB)) %>% 
  mutate(name = case_when(name=="report_SSB" ~ "Original",
                          TRUE ~ "Simulated")) %>% 
  ggplot(aes(year, value, color = name, linetype = name)) + 
  geom_line(linewidth = 1) +
  labs(title = "Model verification",
       y = "Spawning biomass (t)", 
       x = "Year") +
  scale_color_scico_d(palette = 'roma')
  
