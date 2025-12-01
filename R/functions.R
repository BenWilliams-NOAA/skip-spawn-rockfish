run_model <- function(obj_f, data, pars, map_list) {
  # create AD function
  obj <- RTMB::MakeADFun(cmb(obj_f, data), parameters = pars, map = map_list, silent = TRUE)
  # optimize
  fit <- try(nlminb(obj$par, obj$fn, obj$gr, 
                    control = list(iter.max = 100000, eval.max = 20000)), 
             silent = TRUE)
  # return success object
  list(report = obj$report(obj$env$last.par.best), 
       pars = obj$env$parList())
}

project_step <- function(N, F, M, sel, R0) {
  # project population one time step
  n_age = length(N)
  Z <- M + F * sel
  S <- exp(-Z)
  
  # calculate next N
  N_next = numeric(n_age)
  N_next[1] = R0
  N_next[2:(n_age-1)] = N[1:(n_age-2)] * S[1:(n_age-2)]
  # Plus group
  N_next[n_age] = N[n_age-1] * S[n_age-1] + N[n_age] * S[n_age]
  
  # catch in numbers
  C = (F * sel / Z) * N * (1-S)
  
  list(N = N_next, C = C, Z = Z)
}

mismatch_iter <- function(i, rpt, data, bio_mat, func_mat, rec_matrix, obj_f) {
  # setup
  n_years = ncol(rec_matrix)
  n_age = length(rpt$waa)
  proj_years = (max(data$years) + 1):(max(data$years) + n_years)
  M = rpt$M
  q = rpt$q
  ae = diag(nrow = ncol(data$age_error), ncol = n_age)
  cv = 0.1 
  sigma_obs = sqrt(log(cv^2 + 1))
  
  # initial parameter list 
  init_pars <- list(
    log_M = log(M), log_a50C = log(rpt$a50C), deltaC = rpt$deltaC,
    log_a50S = log(rpt$a50S), deltaS = rpt$deltaS, log_q = log(q),
    log_mean_R = rpt$log_mean_R, init_log_Rt = rpt$init_log_Rt, log_Rt = rpt$log_Rt,
    log_mean_F = rpt$log_mean_F, log_Ft = rpt$log_Ft,
    log_F35 = log(rpt$F35), log_F40 = log(rpt$F40), log_F50 = log(rpt$F50),
    log_F35_f = log(rpt$F35), log_F40_f = log(rpt$F40), log_F50_f = log(rpt$F50),
    sigmaR = rpt$sigmaR
  )
  
  # initialize scenarios (i = func maturity, ii = mis-specified)
  # list to hold state for scenarios 
  scen <- list(
    i  = list(data = data, rpt = rpt, pars = init_pars, N_proj = matrix(0, n_age, n_years), C_proj = matrix(0, n_age, n_years)),
    ii = list(data = data, rpt = rpt, pars = init_pars, N_proj = matrix(0, n_age, n_years), C_proj = matrix(0, n_age, n_years))
  )
  # fixed parameters 
  map_list <- list(sigmaR = factor(NA), 
                   log_M = factor(NA), 
                   log_q = factor(NA),
                   log_a50C = factor(NA), deltaC = factor(NA),
                   log_a50S = factor(NA), deltaS = factor(NA),
                   log_F40 = factor(NA),
                   log_F40_f = factor(NA))
  
  # results containers
  res_template <- list(data = vector("list", n_years),
                       report = vector("list", n_years), # stores full report
                       projection = list(tot_bio = numeric(n_years), spawn_bio = numeric(n_years),
                                         spawn_bio_f = numeric(n_years), catch = numeric(n_years), 
                                         recruits = numeric(n_years)))
  results <- list(i = res_template, ii = res_template)
  
  # time loop
  for (y in 1:n_years) {
    # management advice
    F_targets <- list()
    
    for(s in names(scen)) {
      
      rpt = scen[[s]]$rpt 
      
      if(s=="i") {
        ssb_est = as.numeric(tail(rpt$spawn_bio_f, 1))
        F_targets[[s]] = tier_3(ssb_est, as.numeric(rpt$B40_f), as.numeric(rpt$F40_f)) * .4
      } else {
        ssb_est = as.numeric(tail(rpt$spawn_bio, 1))
        F_targets[[s]] = tier_3(ssb_est, as.numeric(rpt$B40), as.numeric(rpt$F40)) * .4
      }
      
      # project population
      # determine biology
      # N from previous year
      if (y == 1) {
        N_prev = rpt$Nat[, ncol(rpt$Nat)]
      } else {
        N_prev = scen[[s]]$N_proj[, y-1]
      }
      spawn_fract = rpt$spawn_fract 
      spawn_adj = exp(-M)^(spawn_fract)
      
      # save projections
      results[[s]]$projection$tot_bio[y] <- sum(N_prev * rpt$waa)
      results[[s]]$projection$spawn_bio[y] <- sum(N_prev * spawn_adj * rpt$waa * bio_mat[, y] * 0.5)
      results[[s]]$projection$spawn_bio_f[y] <- sum(N_prev * spawn_adj * rpt$waa * func_mat[, y] * 0.5)
      
      # project 1-step
      step = project_step(N_prev, F_targets[[s]], M, rpt$slx[, 1], rec_matrix[i, y])
      # update state
      scen[[s]]$N_proj[, y] = step$N
      scen[[s]]$C_proj[, y] = step$C
      
      # store results
      results[[s]]$projection$catch[y] <- sum(step$C * rpt$waa)
      results[[s]]$projection$recruits[y] <- rec_matrix[i, y]
      
      # generate data 
      N_true = step$N
      true_bio = sum(N_true * rpt$waa)
      is_survey = (proj_years[y] %% 2 == 1)
      srv_bio_true = sum(N_true * rpt$waa * rpt$slx[, 2])
      
      if (is_survey) {
        sim_srv_obs = rlnorm(1, log(srv_bio_true * q) - sigma_obs^2/2, sigma_obs)
        p_srv = N_true * rpt$slx[, 2] * q
        sim_srv_acomp = ae %*% (rmultinom(1, 100, p_srv / sum(p_srv)) / 100)
      }
      
      p_fish = N_true * rpt$slx[, 1]
      sim_fish_acomp = ae %*% (rmultinom(1, 100, p_fish / sum(p_fish)) / 100)
      
      # prepare data 
      d = scen[[s]]$data
      d$years = c(d$years, proj_years[y])
      # update catch & comps
      d$catch_obs = c(d$catch_obs, results[[s]]$projection$catch[y])
      d$catch_ind = c(d$catch_ind, 1)
      d$catch_wt = c(d$catch_wt, tail(d$catch_wt, 1))
      d$fish_age_obs = cbind(d$fish_age_obs, sim_fish_acomp)
      d$fish_age_iss = c(d$fish_age_iss, 100)
      d$fish_size_ind = c(d$fish_size_ind, 0)
      
      if (is_survey) {
        d$srv_obs = c(d$srv_obs, sim_srv_obs)
        d$srv_sd = c(d$srv_sd, sim_srv_obs * cv)
        d$srv_ind = c(d$srv_ind, 1)
        d$srv_age_obs = cbind(d$srv_age_obs, sim_srv_acomp)
        d$srv_age_ind = c(d$srv_age_ind, 1)
        d$srv_age_iss = c(d$srv_age_iss, 100)
        d$fish_age_ind = c(d$fish_age_ind, 0) 
      } else {
        d$srv_ind = c(d$srv_ind, 0)
        d$srv_age_ind = c(d$srv_age_ind, 0)
        d$fish_age_ind = c(d$fish_age_ind, 1)
      }
      
      # define maturity for the assessment
      if (s == "i") {
        # Scenario i: Assessment uses Functional (Matches Truth)
        d$wt_mature = 0.5 * rpt$waa * func_mat[, y]
      } else {
        d$wt_mature = 0.5 * rpt$waa * bio_mat[, y]
      }
      d$wt_mature_f = 0.5 * rpt$waa * func_mat[, y]
      
      scen[[s]]$data <- d
      
      # update parameters
      p <- scen[[s]]$pars
      p$log_Rt <- c(p$log_Rt, 0)
      p$log_Ft <- c(p$log_Ft, tail(p$log_Ft, 1))
      scen[[s]]$pars <- p
      
      # run model
      outcome <- run_model(obj_f, scen[[s]]$data, pars = scen[[s]]$pars, map_list = map_list)
      
      if (!is.null(outcome$report)) {
        scen[[s]]$rpt <- outcome$report
        scen[[s]]$pars <- outcome$pars
        results[[s]]$report[[y]] <- outcome$report
      } else {
        results[[s]]$report[[y]] <- NA
      }
      results[[s]]$data[[y]] <- scen[[s]]$data
    }
  }
  
  # Attach final matrices
  results$i$projection$N <- scen$i$N_proj
  results$ii$projection$N <- scen$ii$N_proj
  
  return(results)
}

omem_parallel <- function(rpt, data, bio_mat, func_mat, rec_matrix, obj_f) {
  require(doParallel)
  require(foreach)
  require(doRNG)
  require(RTMB) # Make sure RTMB is loaded on workers
  
  n_iter = nrow(rec_matrix)
  cl = parallel::makeCluster(parallel::detectCores()-2)
  doParallel::registerDoParallel(cl)
  
  # Export necessary functions
  func_exports <- c("mismatch_iter", "cmb", "run_model", "project_step", "tier_3", "f1")
  
  results <- foreach(i = 1:n_iter,
                     .packages = c("RTMB"),
                     .export = func_exports) %dorng% {
                       
                       tryCatch({
                         mismatch_iter(
                           i = i,
                           rpt = rpt,
                           data = data,
                           bio_mat = bio_mat,
                           func_mat = func_mat,
                           rec_matrix = rec_matrix,
                           obj_f = obj_f  # Pass 'f' (the TMB object) here
                         )
                       }, error = function(e) {
                         return(list(iter = i, error = e$message))
                       })
                     }
  
  parallel::stopCluster(cl)
  return(results)
}

extract_results <- function(res_list, start_year = NULL) {
  require(tidytable)
  require(purrr) 
  
  out <- map_dfr(seq_along(res_list), function(iter_idx) {
    
    iter_res <- res_list[[iter_idx]]
    
    # 1. Error Handling
    if (!is.null(iter_res$error) || !is.null(iter_res$error_msg)) {
      warning(paste("Iteration", iter_idx, "failed and was skipped."))
      return(NULL) 
    }
    
    # 2. Extract Data
    map_dfr(names(iter_res), function(scen_name) {
      
      proj <- iter_res[[scen_name]]$projection
      rpt_list <- iter_res[[scen_name]]$report
      n_years <- length(proj$tot_bio)
      
      # --- Helper 1: For scalar values (Reference Points) ---
      extract_scalar <- function(var_name) {
        map_dbl(rpt_list, function(x) {
          if (is.list(x) && !is.null(x[[var_name]])) {
            return(as.numeric(x[[var_name]]))
          } else {
            return(NA_real_)
          }
        })
      }
      
      # --- NEW Helper 2: For Vector Terminal Values (Ft) ---
      # This extracts the LAST value of the vector for the current year
      extract_terminal <- function(var_name) {
        map_dbl(rpt_list, function(x) {
          if (is.list(x) && !is.null(x[[var_name]])) {
            # Take the last value (tail) of the vector
            return(as.numeric(tail(x[[var_name]], 1)))
          } else {
            return(NA_real_)
          }
        })
      }
      
      # Extract Variables
      ts_B35 <- extract_scalar("B35_f")
      ts_F35 <- extract_scalar("F35_f")
      ts_B40 <- extract_scalar("B40_f")
      ts_F40 <- extract_scalar("F40_f")
      
      B35 <- extract_scalar("B35")
      F35 <- extract_scalar("F35")
      B40 <- extract_scalar("B40")
      F40 <- extract_scalar("F40")
      
      # Extract Terminal F
      ts_Ft  <- extract_terminal("Ft")
      sp  <- extract_terminal("spawn_bio")
      spf  <- extract_terminal("spawn_bio_f")
      
      # Handle Years
      if(is.null(start_year)) {
        yrs <- 1:n_years
      } else {
        yrs <- seq(from = start_year, by = 1, length.out = n_years)
      }
      
      # Create Tidy Table
      tidytable(
        iteration = iter_idx,
        scenario = scen_name,
        year = yrs,
        tot_bio = proj$tot_bio,
        spawn_bio = proj$spawn_bio, 
        spawn_bio_r = sp, 
        spawn_bio_f = proj$spawn_bio_f, 
        spawn_bio_fr = spf, 
        catch = proj$catch,
        recruits = proj$recruits,
        
        # Reference Points
        B35 = B35,
        F35 = F35,
        B40 = B40,
        F40 = F40,
        
        B35_f = ts_B35,
        F35_f = ts_F35,
        B40_f = ts_B40,
        F40_f = ts_F40,
        
        # Fishing Mortality
        Ft = ts_Ft
      )
    })
  })
  
  return(out)
}

flexi_curve <- function(age, skip, smin, smax, type = 'dome', skew = 0, width = 4) {
  # Return 0 for ages outside the [smin, smax] range
  if (age < smin || age > smax) {
    return(0)
  }
  
  # Calculate common values
  midpoint <- (smin + smax) / 2
  spread <- (smax - smin) / width
  
  if (type == 'dome') {
    # Standard dome-shaped curve using a Gaussian-like form
    return(skip * exp(-((age - midpoint) / spread)^2))
    
  } else if (type == 'skewed_dome') {
    # Skewed dome: modifies a symmetric dome by applying a skew factor
    # skew affects how quickly the curve drops off on either side
    dome_value <- skip * exp(-((age - midpoint) / spread)^2)
    
    # Apply skew using a logistic-like tilt around midpoint, bounded for stability
    skew_modifier <- 1 / (1 + exp(-skew * (age - midpoint) / spread))
    
    # Rescale skewed dome back to original scale (optional)
    skewed_value <- dome_value * skew_modifier * 2  # *2 to re-normalize the 0-1 logistic range to 0-2
    return(max(0, min(skewed_value, skip)))
    
  } else if (type == 'increase') {
    # Linearly increasing curve from smin to smax
    return(skip * (age - smin) / (smax - smin))
    
  } else if (type == 'decrease') {
    # Linearly decreasing curve from smin to smax
    return(skip * (1 - (age - smin) / (smax - smin)))
    
  } else {
    stop("Unknown curve type specified. Must be one of: 'dome', 'skewed_dome', 'increase', 'decrease'")
  }
}

# self testing ----
verify_mechanics <- function(rpt, project_step_func) {
  
  n_ages = nrow(rpt$Nat)  
  n_years = ncol(rpt$Nat)
  
  # historical drivers
  hist_N_init = rpt$Nat[, 1] # initial population vector
  hist_rec = rpt$recruits # historical recruitment
  hist_F = rpt$Ft # historical F
  hist_M = rpt$M # natural mortality
  hist_sel = rpt$slx[, 1] # fishery selectivity (use col 1)
  
  # biology 
  waa <- rpt$waa 
  wt_mature <- rpt$wt_mature
  
  # containers
  replay_N <- matrix(NA, nrow = n_ages, ncol = n_years)
  replay_SSB <- numeric(n_years)
  
  # initialize 
  replay_N[, 1] <- hist_N_init
  
  # loop up to n-1, project_step calculates t+1
  for (y in 1:(n_years - 1)) {
    
    # calculate SSB for year y (before projection)
    # check if inputs are vectors (static) or matrices (time-varying)
    waa_y <- if(is.matrix(waa)) waa[, y] else waa
    mat_y <- if(is.matrix(wt_mature)) wt_mature[, y] else wt_mature
    
    spawn_adj <- exp(-hist_M)^(rpt$spawn_fract)
    
    # calculate SSB
    replay_SSB[y] <- sum(replay_N[, y] * mat_y * spawn_adj) 
    
    # run the projection step
    # Note: R0 usually fills the FIRST age of the NEXT year
    # grab recruits[y+1] because that is the recruitment resulting from this step
    step <- project_step_func(N = replay_N[, y],
                              F = hist_F[y],
                              M = hist_M,
                              sel = hist_sel,
                              R0 = hist_rec[y+1])
    
    replay_N[, y+1] <- step$N
  }
  
  # final year SSB
  y_final <- n_years
  waa_last <- if(is.matrix(waa)) waa[, y_final] else waa
  mat_last <- if(is.matrix(wt_mature)) wt_mature[, y_final] else wt_mature
  
  replay_SSB[y_final] <- sum(replay_N[, y_final] *  mat_last * spawn_adj)
  
  # compare
  df <- data.frame(year = rpt$years,
                   report_SSB = rpt$spawn_bio, 
                   replay_SSB = replay_SSB) %>%
    mutate(diff = replay_SSB - report_SSB,
           percent_error = round((diff / report_SSB) * 100, 5)) 
  
  return(df)
}
