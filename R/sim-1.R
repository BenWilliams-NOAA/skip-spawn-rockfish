# notes ----
# average recruitment from recruit timeseries
# observed skip spawning

# load ----
source(here::here("R", "utils.R"))
source(here::here("R", "em.R"))
source(here::here("R", "functions.R"))
source(here::here("R", "hcr.R"))

# data ----
data <- readRDS(here::here("data", "dat.RDS"))
rpt <- readRDS(here::here("output", "rpt.RDS"))
skip <- readRDS(here::here("data", "skip.RDS"))
data$wt_mature_f = data$waa * 0.5 * skip

# globals ----
n_iter <- 15
n_years <- 35 
bio_mat <- matrix(rep(data$maa, n_years), ncol = n_years)
func_mat = matrix(rep(skip, n_years), ncol = n_years)
# recruitment matrix
collapse = 0 # number of years of population collapse (if desired)
set.seed(309)
rec_matrix <- matrix(c(rlnorm(n_iter * collapse, 
                              quantile(log(rpt$recruits), 0.05), 0.1),
                       rlnorm(n_iter * (n_years - collapse), 
                              mean(log(rpt$recruits)), 
                              sd(log(rpt$recruits)))),
                     nrow = n_iter, ncol = n_years)

sim1 = omem_parallel(rpt, data, bio_mat, func_mat, rec_matrix, obj_f = f1)
s1 <- extract_results(sim1)
saveRDS(sim1, here::here("output", "sim1.RDS"))
saveRDS(s1, here::here("output", "s1.RDS"))

# save results
plot_mgmt(s1)
plot_risk(s1)
table_risk(s1)
plot_F_shock(rpt, s1)
plot_ssb_shock(rpt, s1)
plot_catch(rpt, s1)
plot_recruits(rpt, s1)
