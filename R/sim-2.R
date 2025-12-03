# notes ----
# high recruitment from recruit timeseries
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
n_iter <- 50
n_years <- 50 
bio_mat <- matrix(rep(data$maa, n_years), ncol = n_years)
func_mat = matrix(rep(skip, n_years), ncol = n_years)
# recruitment matrix
collapse = 0 # number of years of population collapse (if desired)
set.seed(309)
rec_matrix <- matrix(c(rlnorm(n_iter * collapse, 
                              quantile(log(rpt$recruits), 0.05), 0.1),
                       rlnorm(n_iter * (n_years - collapse), 
                              quantile(log(rpt$recruits), 0.80), 0.2)),
                     nrow = n_iter, ncol = n_years)

sim2 = omem_parallel(rpt, data, bio_mat, func_mat, rec_matrix, obj_f = f1)
s2 <- extract_results(sim2)
saveRDS(sim2, here::here("output", "sim2.RDS"))
saveRDS(s2, here::here("output", "s2.RDS"))

# save results
plot_risk(s2, ssb_col = "spawn_bio_r", bio_ref_col = "B35", F_col = "F35", rmv_yrs = 0) # mgmt view
plot_risk(s2, rmv_yrs = 0) # reality view
table_risk(s2, ssb_col = "spawn_bio_r", bio_ref_col = "B35", F_col = "F35", rmv_yrs = 0) # mgmt view
table_risk(s2, rmv_yrs = 0)
plot_F_shock(rpt, s2)
plot_ssb_shock(rpt, s2)
plot_catch(rpt, s2)
plot_recruits(rpt, s2)
