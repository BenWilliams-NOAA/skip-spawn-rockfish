# northern rockfish simulations
# load ----
source(here::here("R", "utils.R"))
source(here::here("R", "em.R"))
source(here::here("R", "functions.R"))
source(here::here("R", "hcr.R"))

# data ----
data <- readRDS(here::here("data", "goa_nork_dat.RDS"))
rpt <- readRDS(here::here("data", "rpt.RDS"))
skip <- readRDS(here::here("data", "n_dat.RDS")) %>% 
  filter(age %in% 2:51)
data$wt_mature_f = data$waa * 0.5 * skip$functional
max_age = length(data$waa)
smin = 3
smax = 50
skip_levels = c(0.1, 0.2, 0.3)
# plot(2:51, (data$waa * 0.5))
# lines(2:51, data$wt_mature_f)

expand.grid(age = 1:max_age,
            skip = skip_levels) %>%
  mutate(dome = mapply(flexi_curve, age, skip, smin, smax, 'dome'),
         skewed_dome = mapply(flexi_curve, age, skip, smin, smax, 'skewed_dome', skew=1.5),
         increasing = mapply(flexi_curve, age, skip, smin, smax, 'increase'),
         decreasing = mapply(flexi_curve, age, skip, smin, smax, 'decrease'),
         concave = mapply(flexi_curve, age, skip, smin, smax, 'inverse_dome'),
         constant = ifelse(age %in% smin:smax, skip, 0)) -> p_skip

# median recruitment ----
# globals 
n_iter <- 50
n_years <- 100 
bio_mat <- matrix(rep(skip$biological, n_years), ncol = n_years)
func_mat = matrix(rep(skip$functional, n_years), ncol = n_years)

# recruitment matrix
collapse = 0 # number of years of population collapse (if desired)
set.seed(309)
rec_matrix <- matrix(c(rlnorm(n_iter * collapse, 
                              quantile(log(rpt$recruits), 0.05), 0.1),
                       rlnorm(n_iter * (n_years - collapse), 
                              quantile(log(rpt$recruits), 0.50), 0.2)),
                     nrow = n_iter, ncol = n_years)


cl = parallel::makeCluster(parallel::detectCores() - 2)
doParallel::registerDoParallel(cl)
sim1 = omem_parallel(rpt, data, bio_mat, func_mat, rec_matrix, obj_f = f1)
parallel::stopCluster(cl)
s1 <- extract_results(sim1)

# save results
saveRDS(sim1, here::here("output", "sim1.RDS"))
saveRDS(s1, here::here("output", "s1.RDS"))



run_sims(data, 
       shapes = c("dome"), 
              rec_type = "avg",
          skip_levels = 0.1,
         p_skip = p_skip) 



# high recruitment ----
# recruitment matrix
collapse = 0 # number of years of population collapse (if desired)
set.seed(309)
rec_matrix <- matrix(c(rlnorm(n_iter * collapse, 
                              quantile(log(rpt$recruits), 0.05), 0.1),
                       rlnorm(n_iter * (n_years - collapse), 
                              quantile(log(rpt$recruits), 0.80), 0.2)),
                     nrow = n_iter, ncol = n_years)

cl = parallel::makeCluster(parallel::detectCores() - 2)
doParallel::registerDoParallel(cl)
sim2 = omem_parallel(rpt, data, bio_mat, func_mat, rec_matrix, obj_f = f1)
parallel::stopCluster(cl)
s2 <- extract_results(sim2)
glimpse(s2)

# save results
saveRDS(sim2, here::here("output", "sim2.RDS"))
saveRDS(s2, here::here("output", "s2.RDS"))
s2 = readRDS(here::here("output", "s2.RDS"))

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
cl = parallel::makeCluster(parallel::detectCores() - 2)
doParallel::registerDoParallel(cl)
sim3 = omem_parallel(rpt, data, bio_mat, func_mat, rec_matrix, obj_f = f1)
parallel::stopCluster(cl)
s3 <- extract_results(sim3)
glimpse(s3)

# save results
saveRDS(sim3, here::here("output", "sim3.RDS"))
saveRDS(s3, here::here("output", "s3.RDS"))

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

sim4 = omem_parallel(rpt, data, bio_mat, func_mat, rec_matrix, obj_f = f1)
s4 <- extract_results(sim4)
glimpse(s4)

# save results
saveRDS(sim4, here::here("output", "sim4.RDS"))
saveRDS(s4, here::here("output", "s4.RDS"))


  