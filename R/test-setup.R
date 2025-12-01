# notes ----
# MSE test - check to see that is running and reporting the desired values
# ben.williams@noaa.gov
# 2025-12

# load ----
source(here::here("R", "utils.R"))
source(here::here("R", "em.R"))
source(here::here("R", "functions.R"))
source(here::here("R", "hcr.R"))
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

saveRDS(rpt, here::here("output", "rpt.RDS"))

# globals ----
n_iter <- 15
n_years <- 35 
bio_mat <- matrix(rep(data$maa, n_years), ncol = n_years)
skip <- c(0.0016, 0.0033, 0.0069, 0.0115, 0.023, 0.0435, 0.0841, 0.1561, 0.2604, 0.323,
          0.4572, 0.5772, 0.668, 0.7289, 0.7688, 0.8019, 0.8111, 0.8178, 0.816, 0.824,
          0.8351, 0.8297, 0.8335, 0.8412, 0.8418, 0.8443, 0.8493, 0.8505, 0.8481, 0.8532,
          0.8542, 0.8572, 0.8597, 0.8666, 0.858, 0.8697, 0.8726, 0.8687, 0.8761, 0.8742,
          0.8788, 0.8873, 0.8847, 0.89, 0.8873, 0.8939, 0.8976, 0.8941, 0.8958, 0.9066)
saveRDS(skip, here::here("data", "skip.RDS"))
func_mat = bio_mat * skip
data$wt_mature_f = data$waa * 0.5 * skip

# recruitment matrix
collapse = 0 # number of years of population collapse (if desired)
set.seed(309)
rec_matrix <- matrix(c(sample(1:3, n_iter * collapse, replace = TRUE),
                       sample(rpt$recruits, n_iter * (n_years - collapse), replace = TRUE)),
                     nrow = n_iter, ncol = n_years)
rec_matrix = rec_matrix * 1.4 # increase recruitment
sim1 = omem_parallel(rpt, data, bio_mat, func_mat, rec_matrix, obj_f = f1)
s1 <- extract_results(sim1)

# plots
plot_mgmt(s1)
plot_risk(s1)
table_risk(s1)
plot_F_shock(rpt, s1)
plot_ssb_shock(rpt, s1)
plot_catch(rpt, s1)
plot_recruits(rpt, s1)


