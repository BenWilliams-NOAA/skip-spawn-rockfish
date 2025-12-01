cmb <- function(f, d) {
  function(p) f(p, d)
}

f1 <- function(pars, data) {
  require(RTMB)
  RTMB::getAll(pars, data, warn=FALSE)
  
  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")
  
  # setup -------------
  # transform
  M = exp(log_M)
  a50C = exp(log_a50C)
  a50S = exp(log_a50S)
  q = exp(log_q)
  F50 = exp(log_F50)
  F40 = exp(log_F40)
  F35 = exp(log_F35)
  
  F50_f = exp(log_F50_f)
  F40_f = exp(log_F40_f)
  F35_f = exp(log_F35_f)
  
  # spawning adjustments
  spawn_fract = (spawn_mo - 1) / 12
  spawn_adj = exp(-M)^(spawn_fract)
  # index values
  A = nrow(age_error) # number of ages in model
  A1 = length(ages) # number of ages in comps
  T = sum(catch_ind)  # time (fishery years)
  Ts = sum(srv_ind) # time (survey)
  Tfa = sum(fish_age_ind) # time (fishery age comp)
  Tsa = sum(srv_age_ind) # time (survey age comp)
  Tfs = sum(fish_size_ind) # time (fishery size comp)
  L = length(length_bins) # number of length bins
  g = 0.00001 # small number
  
  # containers ---------
  slx = matrix(NA, A, 2) # selectivity - fishery and survey
  Bat = Cat = Nat = Fat = Zat = Sat = matrix(0, A, T)
  initNat = rep(0, A)
  catch_pred = rep(0, T)
  srv_pred = rep(0, Ts)
  srv_var = rep(0,Ts)
  fish_age_pred = matrix(0, A1, Tfa)
  srv_age_pred = matrix(0, A1, Tsa)
  fish_size_pred = matrix(0, L, Tfs)
  spawn_bio_f = spawn_bio = tot_bio = rep(0, T)
  N_spr = sb_spr =  N_spr_f = sb_spr_f = matrix(1, A, 4)
  
  # priors -----------------
  # admb priors do not include a constant that is in dnorm, so removing them
  # nll_M = dnorm(log(M), log(mean_M), cv_M, TRUE)
  # nll_M = -(ll + log(sqrt(2 * pi) * cv_M))
  # nll_q = dnorm(log(q), log(mean_q), cv_q, TRUE)
  # nll_q = -(ll + log(sqrt(2 * pi) * cv_q))
  # nll_q = dnorm(q, mean_q, cv_q, TRUE)
  # nll_sigmaR = dnorm(sigmaR, mean_sigmaR, cv_sigmaR, TRUE)
  
  nll_M = (log(M) - log(mean_M))^2 / (2 * cv_M^2)
  nll_q = (log(q) - log(mean_q))^2 / (2 * cv_q^2)
  nll_sigmaR = (log(sigmaR / mean_sigmaR))^2 / (2 * cv_sigmaR^2)
  
  # analysis ----------
  ## selectivity ----
  adj = ifelse(ages[1]<=1, 1, ages[1]-1)
  sel <- function(a, a50, delta) {
    1. / (1. + exp(-2.944438979 * ((a+adj) - a50) / delta))
  }
  
  for(a in 1:A) {
    slx[a,1] = sel(a, a50C, deltaC)
    slx[a,2] = sel(a, a50S, deltaS)
  }
  
  ## mortality ----
  Ft = exp(log_mean_F + log_Ft)
  for(t in 1:T){
    for(a in 1:A) {
      Fat[a,t] = Ft[t] * slx[a,1]
      Zat[a,t] = Fat[a,t] + M
      Sat[a,t] = exp(-Zat[a,t])
    }
  }
  ## population ----
  # not using this
  ## Bzero
  # initNat[1] = exp(log_mean_R)
  # for (a in 2:A) {
  #   initNat[a] <- initNat[a-1] * exp(-M)
  # }
  # initNat[A] <- initNat[A] / (1 - exp(-M))
  # Bzero = sum(initNat * wt_mature * spawn_adj)
  
  ## Nat ----
  # populate first row
  # need to use correct log_Rt to match ADMB model (the init_log_Rt are in reverse order)
  for(t in 1:T) {
    Nat[1,t] = exp(log_mean_R + log_Rt[t])
  }
  # populate first column
  for(a in 2:(A-1)) {
    Nat[a,1] = exp(log_mean_R - (a-1) * M + init_log_Rt[a-1])
  }
  Nat[A,1] = exp(log_mean_R - (A-1) * M) / (1 - exp(-M))
  
  for(t in 2:T) {
    for(a in 2:A) {
      Nat[a,t] = Nat[a-1,t-1] * Sat[a-1,t-1]
    }
    Nat[A,t] = Nat[A,t] + Nat[A,t-1] * Sat[A,t-1]
  }
  
  recruits = Nat[1,]
  ### spawn_bio ----
  for(t in 1:T) {
    spawn_bio[t] = sum(Nat[,t] * spawn_adj * wt_mature)
    spawn_bio_f[t] = sum(Nat[,t] * spawn_adj * wt_mature_f)
    tot_bio[t] = sum(Nat[,t] * waa)
  }
  
  
  # flag - should be:
  # spawn_bio[T] = sum(Nat[,T] * spawn_adj * wt_mature)
  
  
  ## catch ----
  for(t in 1:T){
    for(a in 1:A){
      Cat[a,t] = Fat[a,t] / Zat[a,t] * Nat[a,t] * (1.0 - Sat[a,t])
    }
    catch_pred[t] = sum(Cat[,t] * waa)
  }
  ssqcatch = sum(catch_wt * (log(catch_obs + g) - log(catch_pred + g))^2)
  
  ## survey biomass ----
  isrv = 1
  srv_like = 0.0
  # survey biomass & likelihood
  for(t in 1:T) {
    if(srv_ind[t]==1) {
      srv_pred[isrv] = sum(Nat[,t] * slx[,2] * waa) * q
      srv_var[isrv] = sqrt(log(1 + srv_sd[isrv]^2 / srv_obs[isrv]^2))
      srv_like = srv_like + log(srv_var[isrv]) + 0.5 * (log(srv_obs[isrv] / srv_pred[isrv]) / srv_var[isrv])^2;
      
      # srv_like = srv_like + sum((log(srv_obs[isrv]) - log(srv_pred[isrv]))^2 /
      # (2 * (srv_sd[isrv] / srv_obs[isrv])^2))
      # srv_like = srv_like + sum((srv_obs[isrv]-srv_pred[isrv])^2/ (2.*(srv_sd[isrv]^2)))
      isrv = isrv + 1
    }
  }
  like_srv = srv_like * srv_wt
  
  ## fishery age comp ----
  fish_age_lk = 0.0
  offset = 0.0
  icomp = 1
  for(t in 1:T) {
    if(fish_age_ind[t] == 1) {
      fish_age_pred[,icomp] = as.numeric(colSums((Cat[,t] / sum(Cat[,t])) * age_error))
      offset = offset - fish_age_iss[icomp] * sum((fish_age_obs[,icomp] + g) *
                                                    log(fish_age_obs[,icomp] + g))
      fish_age_lk = fish_age_lk - sum(fish_age_iss[icomp] * (fish_age_obs[,icomp] + g) *
                                        log(fish_age_pred[,icomp] + g))
      icomp = icomp + 1
    }
  }
  fish_age_lk = fish_age_lk - offset
  like_fish_age = fish_age_lk * fish_age_wt
  
  ## survey age comp ----
  srv_age_lk = 0.0
  offset_sa = 0.0
  icomp = 1
  for(t in 1:T) {
    if(srv_age_ind[t] == 1) {
      srv_age_pred[,icomp] = as.numeric(colSums((Nat[,t] * slx[,2]) /
                                                  sum(Nat[,t] * slx[,2]) * age_error))
      offset_sa = offset_sa - srv_age_iss[icomp] * sum((srv_age_obs[,icomp] + g) *
                                                         log(srv_age_obs[,icomp] + g))
      srv_age_lk = srv_age_lk - srv_age_iss[icomp] * sum((srv_age_obs[,icomp] + g) *
                                                           log(srv_age_pred[,icomp] + g))
      icomp = icomp + 1
    }
  }
  srv_age_lk = srv_age_lk - offset_sa
  like_srv_age = srv_age_lk * srv_age_wt
  
  ## fishery size comp ----
  icomp = 1
  fish_size_lk = 0.0
  offset_fs = 0.0
  for(t in 1:T) {
    if(fish_size_ind[t] == 1) {
      fish_size_pred[,icomp] = as.numeric(colSums((Cat[,t] / sum(Cat[,t])) * size_age))
      offset_fs = offset_fs - fish_size_iss[icomp] * sum((fish_size_obs[,icomp] + g) *
                                                           log(fish_size_obs[,icomp] + g))
      fish_size_lk = fish_size_lk - fish_size_iss[icomp] * sum((fish_size_obs[,icomp] + g) *
                                                                 log(fish_size_pred[,icomp] + g))
      icomp = icomp + 1
    }
  }
  fish_size_lk = fish_size_lk - offset_fs
  like_fish_size = fish_size_lk * fish_size_wt
  
  # SPR ------------------------
  data.frame(log_Rt = log_Rt,
             pred_rec = Nat[1,],
             year = years) -> df
  # filter years 1979:
  df = df[years>=(1977+ages[1]) & years<=(max(years)-ages[1]),]
  n_rec = nrow(df)
  yrs_rec = df$year
  pred_rec = mean(df$pred_rec)
  stdev_rec = sqrt(sum((df$log_Rt - mean(df$log_Rt))^2) / (length(df$log_Rt) - 1))
  
  for(a in 2:A) {
    N_spr[a,1] = N_spr[a-1,1] * exp(-M)
    N_spr[a,2] = N_spr[a-1,2] * exp(-(M + F50 * slx[a-1,1]))
    N_spr[a,3] = N_spr[a-1,3] * exp(-(M + F40 * slx[a-1,1]))
    N_spr[a,4] = N_spr[a-1,4] * exp(-(M + F35 * slx[a-1,1]))
  }
  # plus group
  N_spr[A,1] = N_spr[A-1,1] * exp(-M) / (1 - exp(-M))
  N_spr[A,2] = N_spr[A-1,2] * exp(-(M + F50 * slx[A-1,1])) /
    (1 - exp(-(M + F50 * slx[A,1])))
  N_spr[A,3] = N_spr[A-1,3] * exp(-(M + F40 * slx[A-1,1])) /
    (1 - exp(-(M + F40 * slx[A,1])))
  N_spr[A,4] = N_spr[A-1,4] * exp(-(M + F35 * slx[A-1,1])) /
    (1 - exp(-(M + F35 * slx[A,1])))
  
  # spawning spr
  for(a in 1:A) {
    sb_spr[a,1] = N_spr[a,1] * wt_mature[a] * exp(-spawn_fract * M)
    sb_spr[a,2] = N_spr[a,2] * wt_mature[a] * exp(-spawn_fract * (M + F50 * slx[a,1]))
    sb_spr[a,3] = N_spr[a,3] * wt_mature[a] * exp(-spawn_fract * (M + F40 * slx[a,1]))
    sb_spr[a,4] = N_spr[a,4] * wt_mature[a] * exp(-spawn_fract * (M + F35 * slx[a,1]))
  }
  
  # spr reference points
  SB0 = sum(sb_spr[,1])
  SBF50 = sum(sb_spr[,2])
  SBF40 = sum(sb_spr[,3])
  SBF35 = sum(sb_spr[,4])
  
  # spr penalities
  sprpen = 100. * (SBF50 / SB0 - 0.5)^2
  sprpen = sprpen + 100. * (SBF40 / SB0 - 0.4)^2
  sprpen = sprpen + 100. * (SBF35 / SB0 - 0.35)^2
  
  # scale spr reference points
  B0 = SB0 * pred_rec
  B40 = SBF40 * pred_rec
  B35 = SBF35 * pred_rec
  
  # SPR ------------------------
  
  for(a in 2:A) {
    N_spr_f[a,1] = N_spr_f[a-1,1] * exp(-M)
    N_spr_f[a,2] = N_spr_f[a-1,2] * exp(-(M + F50_f * slx[a-1,1]))
    N_spr_f[a,3] = N_spr_f[a-1,3] * exp(-(M + F40_f * slx[a-1,1]))
    N_spr_f[a,4] = N_spr_f[a-1,4] * exp(-(M + F35_f * slx[a-1,1]))
  }
  # plus group
  N_spr_f[A,1] = N_spr_f[A-1,1] * exp(-M) / (1 - exp(-M))
  N_spr_f[A,2] = N_spr_f[A-1,2] * exp(-(M + F50_f * slx[A-1,1])) /
    (1 - exp(-(M + F50_f * slx[A,1])))
  N_spr_f[A,3] = N_spr_f[A-1,3] * exp(-(M + F40_f * slx[A-1,1])) /
    (1 - exp(-(M + F40_f * slx[A,1])))
  N_spr_f[A,4] = N_spr_f[A-1,4] * exp(-(M + F35_f * slx[A-1,1])) /
    (1 - exp(-(M + F35_f * slx[A,1])))
  
  # spawning spr
  for(a in 1:A) {
    sb_spr_f[a,1] = N_spr_f[a,1] * wt_mature_f[a] * exp(-spawn_fract * M)
    sb_spr_f[a,2] = N_spr_f[a,2] * wt_mature_f[a] * exp(-spawn_fract * (M + F50_f * slx[a,1]))
    sb_spr_f[a,3] = N_spr_f[a,3] * wt_mature_f[a] * exp(-spawn_fract * (M + F40_f * slx[a,1]))
    sb_spr_f[a,4] = N_spr_f[a,4] * wt_mature_f[a] * exp(-spawn_fract * (M + F35_f * slx[a,1]))
  }
  
  # spr reference points
  SB0_f = sum(sb_spr_f[,1])
  SBF50_f = sum(sb_spr_f[,2])
  SBF40_f = sum(sb_spr_f[,3])
  SBF35_f = sum(sb_spr_f[,4])
  
  # spr penalities
  sprpen_f = 100. * (SBF50_f / SB0_f - 0.5)^2
  sprpen_f = sprpen_f + 100. * (SBF40_f / SB0_f - 0.4)^2
  sprpen_f = sprpen_f + 100. * (SBF35_f / SB0_f - 0.35)^2
  
  # scale spr reference points
  B0_f = SB0_f * pred_rec
  B40_f = SBF40_f * pred_rec
  B35_f = SBF35_f * pred_rec
  
  # likelihood/penalties --------------------
  # like_rec = sum(log_Rt^2) / (2 * sigmaR^2) + length(log_Rt) * log(sigmaR) * wt_rec_var
  # like_rec = sum((c(log_Rt, init_log_Rt) + sigmaR * sigmaR / 2)^2) / (2 * sigmaR^2)
  like_rec = sum((c(log_Rt, init_log_Rt) + sigmaR * sigmaR / 2)^2) / (2 * sigmaR^2) * wt_rec_var
  f_regularity = wt_fmort_reg * sum(log_Ft^2)
  
  nll = 0.0
  nll = nll + ssqcatch
  nll = nll + like_srv
  nll = nll + like_fish_age
  nll = nll + like_srv_age
  nll = nll + like_fish_size
  nll = nll + like_rec
  nll = nll + f_regularity
  nll = nll + nll_M
  nll = nll + nll_q
  nll = nll + nll_sigmaR
  nll = nll + sprpen
  nll = nll + sprpen_f
  
  # reports -------------------
  RTMB::REPORT(ages)
  RTMB::REPORT(years)
  RTMB::REPORT(M)
  RTMB::ADREPORT(M)
  RTMB::REPORT(a50C)
  RTMB::REPORT(deltaC)
  RTMB::REPORT(a50S)
  RTMB::REPORT(deltaS)
  RTMB::REPORT(q)
  RTMB::ADREPORT(q)
  RTMB::REPORT(sigmaR)
  RTMB::REPORT(init_log_Rt)
  RTMB::REPORT(log_mean_R)
  RTMB::REPORT(log_Rt)
  RTMB::ADREPORT(log_Rt)
  RTMB::REPORT(log_mean_F)
  RTMB::REPORT(log_Ft)
  RTMB::REPORT(waa)
  RTMB::REPORT(maa)
  RTMB::REPORT(wt_mature)
  RTMB::REPORT(yield_ratio)
  RTMB::REPORT(Fat)
  RTMB::REPORT(Zat)
  RTMB::REPORT(Sat)
  RTMB::REPORT(Cat)
  RTMB::REPORT(Nat)
  RTMB::REPORT(slx)
  RTMB::REPORT(Ft)
  RTMB::REPORT(catch_pred)
  RTMB::REPORT(srv_pred)
  
  RTMB::REPORT(fish_age_pred)
  RTMB::REPORT(srv_age_pred)
  RTMB::REPORT(fish_size_pred)
  
  RTMB::REPORT(tot_bio)
  RTMB::REPORT(spawn_bio)
  RTMB::REPORT(spawn_bio_f)
  RTMB::REPORT(recruits)
  RTMB::ADREPORT(srv_pred)
  RTMB::ADREPORT(tot_bio)
  RTMB::ADREPORT(spawn_bio)
  RTMB::ADREPORT(spawn_bio_f)
  RTMB::ADREPORT(recruits)
  RTMB::REPORT(spawn_fract)
  RTMB::REPORT(B0)
  RTMB::REPORT(B40)
  RTMB::REPORT(B35)
  RTMB::REPORT(F35)
  RTMB::REPORT(F40)
  RTMB::REPORT(F50)
  RTMB::REPORT(B0_f)
  RTMB::REPORT(B40_f)
  RTMB::REPORT(B35_f)
  RTMB::REPORT(F35_f)
  RTMB::REPORT(F40_f)
  RTMB::REPORT(F50_f)
  RTMB::REPORT(pred_rec)
  RTMB::REPORT(n_rec)
  RTMB::REPORT(yrs_rec)
  RTMB::REPORT(stdev_rec)
  
  RTMB::REPORT(ssqcatch)
  RTMB::REPORT(like_srv)
  RTMB::REPORT(like_fish_age)
  RTMB::REPORT(like_srv_age)
  RTMB::REPORT(like_fish_size)
  RTMB::REPORT(like_rec)
  RTMB::REPORT(f_regularity)
  RTMB::REPORT(sprpen)
  RTMB::REPORT(sprpen_f)
  RTMB::REPORT(nll_q)
  RTMB::REPORT(nll_M)
  RTMB::REPORT(nll_sigmaR)
  RTMB::REPORT(nll)
  # nll = 0.0
  return(nll)
}