# notes ----
# crash recruitment from recruit timeseries - followed by median recruits
# observed skip spawning

# load ----
source(here::here("R", "utils.R"))
source(here::here("R", "em.R"))
source(here::here("R", "functions.R"))
source(here::here("R", "hcr.R"))

# data ----
data <- readRDS(here::here("data", "goa_nork_dat.RDS"))
rpt <- readRDS(here::here("output", "rpt.RDS"))
skip <- readRDS(here::here("data", "skip.RDS"))
data$wt_mature_f = data$waa * 0.5 * skip

# globals ----
n_iter <- 50
n_years <- 100 
bio_mat <- matrix(rep(data$maa, n_years), ncol = n_years)
func_mat = matrix(rep(skip, n_years), ncol = n_years)
# recruitment matrix
collapse = 20 # number of years of population collapse (if desired)
set.seed(309)
rec_matrix <- matrix(c(rlnorm(n_iter * collapse, 
                              quantile(log(rpt$recruits), 0.05), 0.1),
                       rlnorm(n_iter * (n_years - collapse), 
                              quantile(log(rpt$recruits), 0.50), 0.2)),
                     nrow = n_iter, ncol = n_years)

sim5 = omem_parallel(rpt, data, bio_mat, func_mat, rec_matrix, obj_f = f1)
s5 <- extract_results(sim5)
saveRDS(sim5, here::here("output", "sim5.RDS"))
saveRDS(s5, here::here("output", "s5.RDS"))
rm(sim5)
gc()
gc()
s5 = readRDS(here::here("output", "s5.RDS"))

# save results
plot_risk(s5, ssb_col = "spawn_bio_r", bio_ref_col = "B35", F_col = "F35", rmv_yrs = 0) +
  ggtitle("Manager's View")# mgmt view
plot_risk(s5, rmv_yrs = 0) # reality view
table_risk(s5, ssb_col = "spawn_bio_r", bio_ref_col = "B35", F_col = "F35", rmv_yrs = 15) # mgmt view
table_risk(s5, rmv_yrs = 0)
plot_F_shock(rpt, s5)
plot_ssb_shock(rpt, s5)
plot_catch(rpt, s5)
plot_recruits(rpt, s5)
plot_management_divergence(s5, rmv_yrs = 0)
calc_risk_prob(s5, rmv_yrs = 0)
plot_risk_worms(s5, rmv_yrs = 0)
calc_bias(s5)
