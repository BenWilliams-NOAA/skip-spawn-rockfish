# dusky rockfish simulations
# load ----
source(here::here("R", "utils.R"))
source(here::here("R", "em.R"))
source(here::here("R", "functions.R"))
source(here::here("R", "hcr.R"))

# data ----
dusk = readRDS(here::here("data", "dusk.rds"))
data = dusk$dat
rpt = readRDS(here::here("output", "d_rpt.rds"))
skip <- readRDS(here::here("data", "d_dat.RDS")) %>% 
  filter(age %in% 4:33)


# globals ----
n_iter <- 50
n_years <- 100
bio_mat <- matrix(rep(skip$biological, n_years), ncol = n_years)
func_mat = matrix(rep(skip$functional, n_years), ncol = n_years)

data$wt_mature_f = data$waa * 0.5 * skip$functional
max_age = length(data$waa)
smin = 4
smax = 30
skip_levels = c(0.1, 0.2, 0.3)
# plot(4:33, (data$waa * 0.5))
# lines(4:33, data$wt_mature_f)

expand.grid(age = 1:max_age,
            skip = skip_levels) %>%
  mutate(dome = mapply(flexi_curve, age, skip, smin, smax, 'dome'),
         skewed_dome = mapply(flexi_curve, age, skip, smin, smax, 'skewed_dome', skew=1.5),
         increasing = mapply(flexi_curve, age, skip, smin, smax, 'increase'),
         decreasing = mapply(flexi_curve, age, skip, smin, smax, 'decrease'),
         concave = mapply(flexi_curve, age, skip, smin, smax, 'inverse_dome'),
         constant = ifelse(age %in% smin:smax, skip, 0)) -> p_skip
shapes <- c("dome", "skewed_dome", "concave", "increasing", "decreasing", "constant")

# median recruitment ----
collapse = 0 # number of years of population collapse 
set.seed(309)
rec_matrix <- matrix(c(rlnorm(n_iter * collapse, 
                              quantile(log(rpt$recruits), 0.05), 0.1),
                       rlnorm(n_iter * (n_years - collapse), 
                              quantile(log(rpt$recruits), 0.50), 0.2)),
                     nrow = n_iter, ncol = n_years)


# setup multicore 
cl = parallel::makeCluster(parallel::detectCores() - 2)
doParallel::registerDoParallel(cl)
sim1d = omem_parallel(rpt, data, bio_mat, func_mat, rec_matrix, obj_f = f1)
# parallel::stopCluster(cl)
s1d <- extract_results(sim1d)

# save results
saveRDS(sim1d, here::here("output", "dusky", "sim1d.RDS"))
saveRDS(s1d, here::here("output", "dusky", "s1d.RDS"))

# run sims
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
    saveRDS(sim_res, here::here("output", "dusky", paste0("d_avg_sim_", run_id, ".RDS")))
    saveRDS(s_res, here::here("output", "dusky",  paste0("d_avg_res_", run_id, ".RDS")))
  }
}

gc()
gc()

# high recruitment ----
# recruitment matrix
collapse = 0 # number of years of population collapse (if desired)
set.seed(309)
rec_matrix <- matrix(c(rlnorm(n_iter * collapse, 
                              quantile(log(rpt$recruits), 0.05), 0.1),
                       rlnorm(n_iter * (n_years - collapse), 
                              quantile(log(rpt$recruits), 0.80), 0.2)),
                     nrow = n_iter, ncol = n_years)

# cl = parallel::makeCluster(parallel::detectCores() - 2)
# doParallel::registerDoParallel(cl)
sim2d = omem_parallel(rpt, data, bio_mat, func_mat, rec_matrix, obj_f = f1)
# parallel::stopCluster(cl)
s2d <- extract_results(sim2d)

# save results
saveRDS(sim2d, here::here("output", "dusky", "sim2d.RDS"))
saveRDS(s2d, here::here("output", "dusky", "s2d.RDS"))

# run sims
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
    saveRDS(sim_res, here::here("output", paste0("d_high_sim_", run_id, ".RDS")))
    saveRDS(s_res, here::here("output", paste0("d_high_res_", run_id, ".RDS")))
  }
}

gc()
gc()
# mean/crash/high recruitment ----
# recruitment matrix
collapse = 20 # number of years of population collapse (if desired)
set.seed(309)
rec_matrix <- matrix(c(rlnorm(n_iter * 20, 
                              quantile(log(rpt$recruits), 0.50), 0.2),
                      rlnorm(n_iter * collapse, 
                              quantile(log(rpt$recruits), 0.05), 0.1),
                       rlnorm(n_iter * (n_years - collapse - 20), 
                              quantile(log(rpt$recruits), 0.80), 0.2)),
                     nrow = n_iter, ncol = n_years) 
# cl = parallel::makeCluster(parallel::detectCores() - 2)
# doParallel::registerDoParallel(cl)
sim3d = omem_parallel(rpt, data, bio_mat, func_mat, rec_matrix, obj_f = f1)
# parallel::stopCluster(cl)
s3d <- extract_results(sim3d)

# save results
saveRDS(sim3d, here::here("output", "dusky", "sim3d.RDS"))
saveRDS(s3d, here::here("output", "dusky", "s3d.RDS"))

# run sims
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
    saveRDS(sim_res, here::here("output", "dusky", paste0("d_mch_sim_", run_id, ".RDS")))
    saveRDS(s_res, here::here("output", "dusky", paste0("d_mch_res_", run_id, ".RDS")))
  }
}

gc()
gc()

# mean/crash/mean recruitment ----
# recruitment matrix
collapse = 20 # number of years of population collapse (if desired)
set.seed(309)
rec_matrix <- matrix(c(rlnorm(n_iter * 20, 
                              quantile(log(rpt$recruits), 0.50), 0.2),
                      rlnorm(n_iter * collapse, 
                              quantile(log(rpt$recruits), 0.05), 0.1),
                       rlnorm(n_iter * (n_years - collapse - 20), 
                              quantile(log(rpt$recruits), 0.50), 0.2)),
                     nrow = n_iter, ncol = n_years) 

sim4d = omem_parallel(rpt, data, bio_mat, func_mat, rec_matrix, obj_f = f1)
s4d <- extract_results(sim4d)

# save results
saveRDS(sim4d, here::here("output", "dusky", "sim4d.RDS"))
saveRDS(s4d, here::here("output", "dusky", "s4d.RDS"))

# run sims
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
    saveRDS(sim_res, here::here("output", paste0("d_mcm_sim_", run_id, ".RDS")))
    saveRDS(s_res, here::here("output", paste0("d_mcm_res_", run_id, ".RDS")))
  }
}
