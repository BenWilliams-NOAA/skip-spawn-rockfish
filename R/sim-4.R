# notes ----
# crash recruitment from recruit timeseries
# observed skip spawning

# load ----
source(here::here("R", "utils.R"))
source(here::here("R", "em.R"))
source(here::here("R", "functions.R"))
source(here::here("R", "hcr.R"))

# data ----
data <- readRDS(here::here("data", "dat.RDS"))
rpt <- readRDS(here::here("output", "rpt.RDS"))
# skip <- readRDS(here::here("data", "skip.RDS"))
# data$wt_mature_f = data$waa * 0.5 * skip

# globals ----
max_age = length(rpt$waa)
skip_levels = c(0.02, 0.05, 0.1, 0.2, 0.3)
smin = 1
smax = 50

expand.grid(age = 1:max_age,
            skip = skip_levels) %>%
  mutate(dome = mapply(flexi_curve, age, skip, smin, smax, 'dome'),
         skewed_dome = mapply(flexi_curve, age, skip, smin, smax, 'skewed_dome', skew=1.5),
         increasing = mapply(flexi_curve, age, skip, smin, smax, 'increase'),
         decreasing = mapply(flexi_curve, age, skip, smin, smax, 'decrease'),
         constant = ifelse(age %in% smin:smax, skip, 0)) -> p_skip


skip =  (1 -( p_skip %>% filter(skip==.02) %>% pull(skewed_dome) )) * data$maa
data$wt_mature_f = data$waa * 0.5 * skip

n_iter <- 50
n_years <- 50 
bio_mat <- matrix(rep(data$maa, n_years), ncol = n_years)
func_mat = matrix(rep(skip, n_years), ncol = n_years)
# recruitment matrix
collapse = 20 # number of years of population collapse (if desired)
set.seed(309)
rec_matrix <- matrix(c(rlnorm(n_iter * collapse, 
                              quantile(log(rpt$recruits), 0.05), 0.1),
                       rlnorm(n_iter * (n_years - collapse), 
                              quantile(log(rpt$recruits), 0.80), 0.2)),
                     nrow = n_iter, ncol = n_years)

shapes <- c("dome", "skewed_dome", "increasing", "decreasing", "constant")

for (shp in shapes) {
  for (lvl in skip_levels) {
    run_id <- paste0(shp, "_", lvl) # e.g., "dome_0.02"
    message(paste("Running simulation for:", run_id))
    
    # extract curve values
    current_curve_vals <- p_skip %>% 
      filter(skip == lvl) %>% 
      pull(all_of(shp))
    
    current_skip_vec = (1 - current_curve_vals) * data$maa
    
    # update the data object for this specific iteration
    data$wt_mature_f = data$waa * 0.5 * current_skip_vec
    func_mat = matrix(rep(current_skip_vec, n_years), ncol = n_years)
    
    # run model
    sim_res = omem_parallel(rpt, data, bio_mat, func_mat, rec_matrix, obj_f = f1)
    s_res = extract_results(sim_res)
    
    # results
    saveRDS(sim_res, here::here("output", paste0("crash_sim_", run_id, ".RDS")))
    saveRDS(s_res, here::here("output", paste0("crash_res_", run_id, ".RDS")))
    gc()
    gc()
  }
}

