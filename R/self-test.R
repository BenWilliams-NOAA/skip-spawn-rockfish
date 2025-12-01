# notes ----
# MSE self test - check to replicate the last assessment report
# making sure the 'project_step' function is correct
# ben.williams@noaa.gov
# 2025-12

# load ----
source(here::here("R", "utils.R"))
source(here::here("R", "em.R"))
# data ----
data <- readRDS(here::here("data", "dat.RDS"))
pars <- readRDS(here::here("data", "pars.RDS"))

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
test_project <- verify_mechanics(rpt, project_step)

test_project %>% 
  # filter(year!=2024) %>% 
  ggplot(aes(x = year)) +
  geom_line(aes(y = report_SSB, color = "Original")) +
  geom_line(aes(y = replay_SSB, color = "MSE"), linetype = "dashed") +
  labs(title = "project_step verification",
       y = "Spawning Biomass", 
       x = "Year", color = "Source")
