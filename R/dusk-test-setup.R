# notes ----
# MSE test - check to see that is running and reporting the desired values
# this timwe for dusky rockfish
# ben.williams@noaa.gov
# 2025-12

# load ----
source(here::here("R", "utils.R"))
source(here::here("R", "em.R"))
source(here::here("R", "functions.R"))
source(here::here("R", "hcr.R"))
# data ----
dusk = readRDS(here::here("data", "dusk.rds"))
data = dusk$dat
fit = dusk$fit
pars = dusk$obj$env$parList()
map = dusk$obj$env$map
nms = unique(names(fit$par))
split_list = split(fit$par, names(fit$par))
pars = lapply(split_list, unname)
pars = pars[nms]
# put any mapped items back into the pars
if(!is.null(map)) {
  pars[[names(map)]] = str()[[names(map)]]
}
data <- readRDS(here::here("data", "dat.RDS"))
pars <- readRDS(here::here("data", "pars.RDS"))

# base maturity
# duplicate the biological maturity for this test case, also add necessary parameters
data$wt_mature_f = data$wt_mature
pars$log_F50_f = pars$log_F40_f = pars$log_F35_f = 0

# run base model
obj <- try(RTMB::MakeADFun(cmb(f1, data),
                           parameters = pars,
                           map = list(log_M = factor(NA))))
fit = nlminb(start = obj$par,
             objective = obj$fn,
             gradient = obj$gr,
             control = list(iter.max = 100000,
                            eval.max = 20000))
rpt <- obj$report(obj$env$last.par.best)

saveRDS(rpt, here::here("output", "d_rpt.RDS"))

