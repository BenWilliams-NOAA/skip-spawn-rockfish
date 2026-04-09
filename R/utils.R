library(tidytable)
library(ggplot2)
library(foreach)
library(doParallel)
library(RTMB)
library(scico)
# remotes::install_git("BenWilliams-NOAA/afscassess")
theme_set(afscassess::theme_report())
library(ggpubr)

# plot functions ----

plot_F_shock <- function(rpt, sims) {
  bind_rows(data.frame(year = rpt$years, 
                       F = rpt$Ft, 
                       type = "historical",
                       iteration = 0),
            data.frame(year = max(rpt$years) + sims$year,
                       F = sims$Ft, 
                       type = "projection", 
                       scenario = sims$scenario,
                       iteration = sims$iteration) %>% 
              filter(scenario=="i")) %>% 
    ggplot(aes(year, F, color = type, group = iteration)) + 
    geom_line() + 
    geom_vline(xintercept = max(rpt$years), linetype =3) +
    labs(title = "The 'F' Shock", subtitle = "Did F jump at the start of the simulation?")
}

plot_ssb_shock <- function(rpt, sims, ssb_col = "spawn_bio_f") {
  bind_rows(data.frame(year = rpt$years, 
                       spawn_bio = rpt[[ssb_col]], 
                       scenario = "historical",
                       iteration = 0),
            data.frame(year = max(rpt$years) + sims$year,
                       spawn_bio = sims[[ssb_col]], 
                       # type = "projection", 
                       scenario = sims$scenario,
                       iteration = sims$iteration)) %>% 
    ggplot(aes(year, spawn_bio, group = interaction(scenario, iteration), color = scenario)) + 
    geom_line() + 
    geom_vline(xintercept = max(rpt$years), linetype =3) +
    labs(title = "The 'SSB' Shock", subtitle = "Did SSB drop at the start of the simulation?") +
    expand_limits(y = 0)
}

plot_catch <- function(rpt, sims) {
  bind_rows(data.frame(year = rpt$years, 
                       catch = rpt$catch_pred, 
                       scenario = "historical",
                       iteration = 0),
            data.frame(year = max(rpt$years) + sims$year,
                       catch = sims$catch, 
                       scenario = sims$scenario,
                       iteration = sims$iteration)) %>% 
    mutate(scenario = ifelse(scenario == "i", "Functional", "Biological")) %>% 
    ggplot(aes(year, catch, group = interaction(scenario, iteration), color = scenario)) + 
    geom_line() + 
    geom_vline(xintercept = max(rpt$years), linetype =3) +
    scale_color_scico_d(name = "", palette='roma')
}

plot_recruits <- function(rpt, sims) {
  bind_rows(data.frame(year = rpt$years, 
                       recruits = rpt$recruits, 
                       scenario = "historical",
                       iteration = 0),
            data.frame(year = max(rpt$years) + sims$year,
                       recruits = sims$recruits, 
                       scenario = "simulated",
                       iteration = sims$iteration)) %>% 
    ggplot(aes(year, recruits, group = interaction(scenario, iteration), color = scenario)) + 
    geom_line() + 
    stat_summary(data = . %>% filter(scenario=='simulated'), 
                 fun = median, 
                 geom = "line",
                 aes(group = scenario), 
                 color = 1) +
    geom_vline(xintercept = max(rpt$years), linetype =3) +
    scico::scale_color_scico_d(palette = "grayC", end = 0.8)
}

plot_risk <- function(sims, rmv_yrs = 0, ssb_col = "spawn_bio_f", bio_ref_col = "B35_f", F_col = "F35_f") {
  # manager's view
  sims %>%
    mutate(scenario = ifelse(scenario == "i", 
                             "Correct", 
                             "Misspecified")) %>% 
    filter(year >= (min(year) + rmv_yrs)) %>% # Analyze long-term equilibrium
    arrange(scenario, iteration, year) %>%
    mutate(true_status = .data[[ssb_col]] / .data[[bio_ref_col]], 
           true_intensity = Ft / .data[[F_col]]) %>%
    ggplot(aes(true_status, true_intensity)) + 
    annotate("rect", xmin = -Inf, xmax = 1, ymin = 1, ymax = Inf, 
             fill = "#a50f15", alpha = 0.5) +
    annotate("rect",xmin = 1, xmax = Inf, ymin = -Inf, ymax = 1, 
             fill = "#b8de64", alpha = 0.5) +
    annotate("rect",xmin = 1, xmax = Inf, ymin = 1, ymax = Inf, 
             fill = "#f5f77c", alpha = 0.5) +
    annotate("rect",xmin = -Inf, xmax = 1, ymin = -Inf, ymax = 1, 
             fill = "#E3AF27", alpha = 0.5) +
    # geom_density_2d_filled(contour_var = "ndensity", alpha = 0.2) +
    geom_density_2d(color="gray40") +# target crosshairs
    geom_vline(xintercept = 1, linetype = "dashed", color = "white") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "white") +
    geom_point(aes(color=iteration)) +
    facet_wrap(~scenario) +
    labs(title = "Risk Comparison",
         # subtitle = "Consequences of ignoring skip spawning",
         x = "True SSB / True B35",
         y = "Realized F / True F35") +
    theme(legend.position = "none") + 
    expand_limits(x = c(0, 2), y = c(0, 2))
}

plot_risk_true <- function(sims, rmv_yrs = 0, ssb_col = "spawn_bio_f", bio_ref_col = "B35_f", F_col = "F35_f") {
  
  sims %>%
    mutate(scenario = ifelse(scenario == "i", 
                             "Mgmt: Functional (Correct)", 
                             "Mgmt: Biological (Misspecified)")) %>% 
    filter(year >= (min(year) + rmv_yrs)) %>% 
    arrange(scenario, iteration, year) %>%
    
    # --- KEY CHANGE HERE ---
    # Numerator changed from 'Ft' (Estimated) to 'F_target_val' (True)
    mutate(true_status = .data[[ssb_col]] / .data[[bio_ref_col]], 
           true_intensity = F_target_val / .data[[F_col]]) %>% 
    # -----------------------
  
  ggplot(aes(true_status, true_intensity)) + 
    # (Kobe Plot Background)
    annotate("rect", xmin = -Inf, xmax = 1, ymin = 1, ymax = Inf, fill = "red", alpha = 0.1) +
    annotate("rect",xmin = 1, xmax = Inf, ymin = -Inf, ymax = 1, fill = "#CFE3A3", alpha= 0.2) +
    annotate("rect",xmin = 1, xmax = Inf, ymin = 1, ymax = Inf, fill = "yellow", alpha = 0.1) +
    annotate("rect",xmin = -Inf, xmax = 1, ymin = -Inf, ymax = 1, fill = "orange", alpha = 0.1) +
    
    geom_density_2d(color="gray40") +
    geom_vline(xintercept = 1, linetype = "dashed", color = "white") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "white") +
    geom_point(aes(color=iteration)) +
    facet_wrap(~scenario) +
    
    # Updated Labels
    labs(title = "True Biological Risk (The 'God's Eye' View)",
         subtitle = "True Stock Status (X) vs. True Implemented F (Y)",
         x = "True SSB / Limit (B35)",
         y = "True F / Limit (F35)") +
    theme(legend.position = "none") + 
    expand_limits(x = c(0, 2), y = c(0, 2))
}
# risk table
# calculate true ratios (same as plot)
# scenario i reference points for the denominator
table_risk <- function(sims, rmv_yrs = 0, ssb_col = "spawn_bio_f", bio_ref_col = "B35_f", F_col = "F35_f") {
  sims %>%
    filter(year >= (min(year) + rmv_yrs)) %>% 
    mutate(x_ratio = .data[[ssb_col]] / .data[[bio_ref_col]],
           y_ratio = Ft / .data[[F_col]],
           zone = case_when(
             x_ratio < 1 & y_ratio > 1  ~ "Red (Overfished & Overfishing)",
             x_ratio > 1 & y_ratio < 1  ~ "Green (Safe)",
             x_ratio < 1 & y_ratio <= 1 ~ "Orange (Overfished)",
             x_ratio >= 1 & y_ratio >= 1 ~ "Yellow (Overfishing)")) %>%
    summarise(count = n(), .by=c(scenario, zone)) %>%
    mutate(total = sum(count),
           percentage = round((count / total) * 100, 1),
           .by = scenario) %>%
    select(scenario, zone, percentage) %>%
    pivot_wider(names_from = scenario, values_from = percentage, values_fill = 0)
}

plot_true_reality <- function(sims, rmv_yrs = 0) {
  sims %>%
    filter(year >= (min(year) + rmv_yrs)) %>%
    # Numerator: F_target_val (True Implemented F)
    # Denominator: F35_f (We use the Est. Functional limit as the best proxy for True Limit 
    # unless you have the True OM limit saved elsewhere)
    mutate(true_intensity = F_target_val / F35_f, 
           true_status = spawn_bio_f / B35_f) %>% 
    ggplot(aes(true_status, true_intensity)) +
    # ... (rest of plotting code same as before) ...
    geom_point(aes(color=iteration)) +
    labs(title = "True Implementation Reality",
         subtitle = "True Implemented F / Limit",
         y = "True F / Functional Limit") +
    facet_wrap(~scenario) +
    expand_limits(y=1)
}
# plot_true_reality(s1, rmv_yrs = 0)

calc_bias <- function(sims, rmv_yrs = 0) {

  sims %>%
    filter(year >= (min(year) + rmv_yrs)) %>%
    mutate(true_intensity = Ft / F35_f) %>%
    summarise(
      mean_intensity = mean(true_intensity, na.rm = TRUE),
      median_intensity = median(true_intensity, na.rm = TRUE),
      .by = scenario
    ) %>%
    mutate(
      bias_pct = (mean_intensity - mean_intensity[scenario == "i"]) / mean_intensity[scenario == "i"] * 100
    ) %>% 
    filter(scenario == "ii")

}

calc_bias2 <- function(sims, rmv_yrs = 0) {
  
  sims %>%
    filter(year >= (min(year) + rmv_yrs)) %>%
    mutate(true_intensity = Ft / F35_f) %>%
    summarise(
      mean_intensity = mean(true_intensity, na.rm = TRUE),
      median_intensity = median(true_intensity, na.rm = TRUE),
      .by = c(scenario, year, iteration)
    ) %>%
    mutate(
      bias_pct = (mean_intensity - mean_intensity[scenario == "i"]) / 
        mean_intensity[scenario == "i"] * 100,
    ) %>% 
    filter(scenario == "ii")
  
}


# calc_bias2(s1) %>% 
#   mutate(recruits = "average") %>% 
#   bind_rows(calc_bias2(s2) %>% 
#   mutate(recruits = "high")) %>% 
#   bind_rows(calc_bias2(s3) %>% 
#   mutate(recruits = "crash")) %>% 
#   ggplot(aes(year, bias_pct, group = interaction(iteration, recruits), color = recruits)) + 
#            # geom_line(alpha = 0.00) +
#   geom_hline(yintercept = 0, lty = 3)  +
#   stat_summary(aes(group = recruits), fun.y = median, geom = "line", size = 1) +
#   scico::scale_color_scico_d(palette = 'roma', end = 0.7) +
#   coord_cartesian(ylim = c(-10, 30))

# calc_bias2(avg_res_constant_0.3) %>% 
#   mutate(recruits = "average") %>% 
#   bind_rows(calc_bias2(high_res_constant_0.3) %>% 
#               mutate(recruits = "high")) %>% 
#   bind_rows(calc_bias2(crash_res_constant_0.3) %>% 
#               mutate(recruits = "crash")) %>% 
#   ggplot(aes(year, bias_pct, group = interaction(iteration, recruits), color = recruits)) + 
#   geom_line() +
#   geom_hline(yintercept = 0, lty = 3) 


# (calc_bias2(s1) %>% 
#   ggplot(aes(year, bias_pct, group = iteration)) + 
#   geom_line(alpha = 0.2) +
#   geom_hline(yintercept = 0, lty = 3) +
#   stat_summary(aes(group = 1), fun.y = median, geom = "line", size = 1)) / 
# 
# (calc_bias2(s2) %>% 
#   ggplot(aes(year, bias_pct, group = iteration)) + 
#   geom_line(alpha = 0.2) +
#   geom_hline(yintercept = 0, lty = 3) +
#   stat_summary(aes(group = 1), fun.y = median, geom = "line", size = 1)) /
#   
#   (calc_bias2(s3) %>% 
#      ggplot(aes(year, bias_pct, group = iteration)) + 
#      geom_line(alpha = 0.2) +
#      geom_hline(yintercept = 0, lty = 3) +
#      stat_summary(aes(group = 1), fun.y = median, geom = "line", size = 1)) +
#   plot_layout(axis = 'collect') 

plot_management_divergence <- function(sims, rmv_yrs = 0) {
  
  sims %>%
    filter(year >= (min(year) + rmv_yrs)) %>%
    # Standardize everything to the TRUE Functional Reality
    mutate(true_status = spawn_bio_f / B35_f,
           true_intensity = Ft / F35_f,
           scenario_label = ifelse(scenario == "i", "Correct (Functional)", "Misspecified (Biological)")) %>%
    ggplot(aes(true_status, true_intensity, color = scenario_label, fill = scenario_label)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") + # B35 limit
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") + # F35 limit
    geom_point(alpha = 0.3) +
    geom_smooth(method = "gam", formula = y ~ s(x), se = TRUE, alpha = 0.2) +
    scale_color_manual(values = c("Correct (Functional)" = "steelblue", 
                                  "Misspecified (Biological)" = "firebrick")) +
    scale_fill_manual(values = c("Correct (Functional)" = "steelblue", 
                                 "Misspecified (Biological)" = "firebrick")) +
    labs(title = "Realized Control Rules",
         x = "True Stock Status",
         y = "True Fishing Intensity",
         color = "Management Scenario",
         fill = "Management Scenario") +
    theme_bw() +
    theme(legend.position = "top") 
    # coord_cartesian(ylim = c(0, 1.2), xlim = c(0.5, 2.0)) # Focus on the relevant range
}

calc_risk_prob <- function(sims, rmv_yrs = 0) {
  require(dplyr)
  require(tidytable)
  
  risk_table <- sims %>%
    filter(year >= (min(year) + rmv_yrs)) %>%
    mutate(
      # Define breaches based on True Functional Biomass
      # Breach Limit (B35): Stock is below the functional limit
      breach_B40 = ifelse(spawn_bio_f < B40_f, 1, 0),
      breach_B35 = ifelse(spawn_bio_f < B35_f, 1, 0),
      
      # Breach Hard Limit (B20): Stock is severely depleted (approaching collapse)
      # Note: We estimate B20 as roughly B35 / 1.75 since B35 is proxy for BMSY
      # or explicitly if you have B20 in your report. Here we use B35 * (20/35)
      breach_B20 = ifelse(spawn_bio_f < (B35_f * (20/35)), 1, 0)
    ) %>%
    summarise(
      # Calculate probability (Mean of binary 0/1)
      prob_below_B40 = mean(breach_B40),
      prob_below_B35 = mean(breach_B35),
      prob_below_B20 = mean(breach_B20),
      .by = scenario
    ) %>%
    mutate(
      scenario = ifelse(scenario == "i", "Correct (Functional)", "Misspecified (Biological)"),
      # Convert to percentage for readability
      prob_below_B40 = round(prob_below_B40 * 100, 1),
      prob_below_B35 = round(prob_below_B35 * 100, 1),
      prob_below_B20 = round(prob_below_B20 * 100, 1)
    ) %>%
    select(scenario, prob_below_B40, prob_below_B35, prob_below_B20)
  
  return(risk_table)
}

plot_risk_worms <- function(sims, rmv_yrs = 0) {
  require(dplyr)
  require(ggplot2)
  require(tidytable)
  
  # 1. Prepare Data & Define Status
  df <- sims %>%
    filter(year >= (min(year) + rmv_yrs)) %>%
    mutate(
      # Define thresholds based on Functional B35 (The Truth)
      target_B40 = B40_f,
      limit_B35 = B35_f,
      target_F35 = F40_f,
      limit_F35 = F35_f,
      # limit_B20 = B35_f * (20/35),
      
      # Assign Categorical Status
      status = case_when(
        spawn_bio_f < limit_B35 & Ft > limit_F35 ~ "High Risk (<B35 & >F35)",
        spawn_bio_f < limit_B35 ~ "Overfished (< B35)",
        Ft > target_F35 ~  "Overfishing (> F35)",
        # target_B40 > spawn_bio_f & spawn_bio_f > limit_B35 ~ "Approaching Overfished",
        TRUE ~ "Safe (> B35)"
      ),
      scenario_label = ifelse(scenario == "i", "Correct", "Misspecified")
    )
  
  rankings <- df %>%
    summarise(severity_score = sum(ifelse(status == "Critical (< B20)", 10, 
                                          ifelse(status == "Overfished (< B35)", 1, 0))), 
              .by = c(scenario, iteration)) %>%
      mutate(iter_rank = row_number(), .by = scenario) %>% 
    arrange(iter_rank)
  
  df_plot <- df %>%
    left_join(rankings, by = c("scenario", "iteration"))
  
  ggplot(df_plot, aes(year, iter_rank, fill = status)) +
    geom_raster() + # Raster is faster than tile for large grids
    scale_fill_manual(values = c("Safe (> B35)" = "#b8de64",       # Pale Green
                                 "Overfished (< B35)" = "#E3AF27", # Orange
                                 # "Approaching Overfished" = "#f3cb7c" ,
                                 "Overfishing (> F35)" = "#f5f77c",
                                 "High Risk (<B35 & >F35)" = "#a50f15")) + # Deep Red
    facet_wrap(~scenario_label) +
    labs(title = "Temporal Patterns",
         # subtitle = "Each row is one 50-year simulation run",
         x = "Year",
         y = "Iteration",
         fill = "True Stock Status") +
    afscassess::theme_report() +
    theme(
      legend.position = "bottom", legend.box="vertical", legend.margin=margin()
    ) 
}

summarise_performance <- function(sims, rmv_yrs = 0) {
  require(tidytable)
  
  sims %>%
    filter(year >= (min(year) + rmv_yrs)) %>%
    # Calculate AAV (Average Annual Variability in Catch) per iteration first
    arrange(scenario, iteration, year) %>% 
    mutate(
      catch_diff = abs(catch - lag(catch)),
      # Avoid division by zero if catch is 0
      catch_lag = ifelse(lag(catch) == 0, 0.001, lag(catch)), 
      ann_var = catch_diff / catch_lag
    ) %>%
    summarise(
      # --- Yield Metrics ---
      avg_catch = mean(catch),
      median_catch = median(catch),
      
      # --- Stability Metrics (AAV) ---
      # Mean of the annual % change in catch
      catch_cv = mean(ann_var, na.rm = TRUE), 
      
      # --- Status Metrics (True Functional Reality) ---
       # avg_depletion = mean(spawn_bio_f / B0_f), # Relative to Virgin

      # --- Risk Probability ---
      # Prob SSB < True B35
      risk_b35 = mean(spawn_bio_f < B35_f),
      # Prob SSB < True B20 (Severe Depletion)
      risk_b20 = mean(spawn_bio_f < (B35_f * (20/35))),
      
      .by = c(scenario)
    ) %>%
    mutate(
      scenario_label = ifelse(scenario == "i", "Correct (Functional)", "Misspecified (Biological)"),
      # Format percentages for readability
      catch_cv_pct = paste0(round(catch_cv * 100, 1), "%"),
      risk_b35_pct = paste0(round(risk_b35 * 100, 1), "%"),
      risk_b20_pct = paste0(round(risk_b20 * 100, 1), "%")
    ) %>%
    select(scenario_label, avg_catch, catch_cv_pct,  risk_b35_pct, risk_b20_pct)
}

calc_estimation_error <- function(sims, rmv_yrs = 0) {
  sims %>%
    filter(year >= (min(year) + rmv_yrs)) %>%
    mutate(
      # Relative Error for SSB
      err_ssb = (spawn_bio_r - spawn_bio_f) / spawn_bio_f,
      # Relative Error for F
      err_f   = (Ft - F_target_val) / F_target_val 
      # Note: For F, we compare Estimated F (Ft) to the Realized F (F_target_val)
    ) %>%
    summarise(
      # MRE = Directional Bias (Positive = Overestimation)
      MRE_SSB = mean(err_ssb, na.rm = TRUE),
      MRE_F   = mean(err_f, na.rm = TRUE),
      
      # MARE = Magnitude of Error (Absolute value)
      MARE_SSB = mean(abs(err_ssb), na.rm = TRUE),
      
      # RMSE = Penalizes large errors more heavily
      RMSE_SSB = sqrt(mean(err_ssb^2, na.rm = TRUE)),
      
      .by = scenario
    ) %>%
    mutate(
      scenario = ifelse(scenario == "i", "Correct", "Misspecified"),
      across(where(is.numeric), ~ round(.x, 3))
    )
}

compare_time_blocks <- function(sims, rmv_yrs = 20) {
  sims %>%
    filter(year >= (min(year) + rmv_yrs)) %>%
    mutate(
      period = case_when(
        year <= 30 ~ "Short Term (1-10 yrs)",
        year > 30  ~ "Long Term (30+ yrs)",
        TRUE ~ "Transition"
      )
    ) %>% 
    filter(period != "Transition") %>%
    summarise(
      avg_catch = mean(catch),
      avg_ssb = mean(spawn_bio_f),
      risk_b35 = mean(spawn_bio_f < B35_f),
      risk_f40 = mean(Ft > F40_f),
      risk_f35 = mean(Ft > F35_f),
      .by = c(scenario, period)
    ) %>%
    arrange(period, scenario) %>%
    mutate(
      scenario = ifelse(scenario == "i", "Correct", "Misspecified")
    )
}

calc_rebuilding <- function(sims, rmv_yrs = 0, time_frame = 10) {
  sims %>%
    filter(year >= (min(year) + rmv_yrs)) %>%
    mutate(
      is_depleted = spawn_bio_f < B35_f,
      is_rebuilt  = spawn_bio_f >= B40_f,
      yr_crashed_temp = as.numeric(suppressWarnings(min(year[is_depleted], na.rm = TRUE))),
      .by = c(scenario, iteration)
    ) %>%
    summarise(
      yr_crashed = first(yr_crashed_temp),
      yr_recovered = as.numeric(suppressWarnings(min(year[year > first(yr_crashed_temp) & is_rebuilt], na.rm = TRUE))),
      .by = c(scenario, iteration)
    ) %>%
    mutate(
      # Use NA_real_ to maintain the numeric/double column type
      yr_crashed = ifelse(is.infinite(yr_crashed), NA_real_, yr_crashed),
      yr_recovered = ifelse(is.infinite(yr_recovered), NA_real_, yr_recovered),
      
      time_to_rebuild = yr_recovered - yr_crashed,
      rebuilt_in_time = ifelse(!is.na(time_to_rebuild) & time_to_rebuild <= time_frame, 1, 0)
    ) %>%
    filter(!is.na(yr_crashed)) %>%
    summarise(
      n_crashes = n(),
      mean_yrs_to_rebuild = round(mean(time_to_rebuild, na.rm = TRUE), 1),
      prob_rebuilt_in_time = round(mean(rebuilt_in_time, na.rm = TRUE) * 100, 1),
      .by = scenario
    ) %>%
    mutate(
      scenario = ifelse(scenario == "i", "Correct (Functional)", "Misspecified (Biological)"),
      prob_rebuilt_in_time = paste0(prob_rebuilt_in_time, "%")
    )
}
# calc_rebuilding(d3, time_frame = 21)
# calc_rebuilding(s1, rmv_yrs = 15)
# calc_rebuilding(s2, rmv_yrs = 15)
# calc_rebuilding(s3, rmv_yrs = 15)
# 
# calc_rebuilding(ac3)
# calc_rebuilding(ad3, rmv_yrs = 15)

plot_tradeoff_frontier <- function(sims, rmv_yrs = 0) {
  
  # Calculate Risk and Yield per iteration
  tradeoff_data <- sims %>%
    filter(year >= (min(year) + rmv_yrs)) %>%
    summarise(
      avg_catch = mean(catch, na.rm = TRUE),
      # Risk: Probability of severe depletion (using B20 as a hard limit proxy)
      risk_b20 = mean(spawn_bio_f < (B35_f * (20/35)), na.rm = TRUE), 
      .by = c(scenario, iteration)
    ) %>%
    mutate(scenario_label = ifelse(scenario == "i", "Correct (Functional)", "Misspecified (Biological)"))
  
  ggplot(tradeoff_data, aes(risk_b20, avg_catch, color = scenario_label)) +
    # Plot individual iterations as transparent points
    geom_point(alpha = 0.5, size = 2) +
    # Add 95% confidence ellipses to show the spread/variance
    stat_ellipse(level = 0.95, linetype = "dashed", alpha = 0.7) + 
    # Add a large diamond centroid for the overall scenario mean
    stat_summary(fun = mean, geom = "point", size = 5, shape = 18, color = "black") +
    stat_summary(fun = mean, geom = "point", size = 4, shape = 18, aes(color = scenario_label)) +
    
    scale_color_manual(values = c("Correct (Functional)" = "steelblue", 
                                  "Misspecified (Biological)" = "firebrick")) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      title = "Management Trade-off Frontier",
      subtitle = "Average Yield vs. Biological Risk (50 Iterations)",
      x = "Probability of Severe Depletion (True SSB < B20)",
      y = "Average Annual Catch",
      color = "Management Scenario"
    ) +
    # theme_bw() +
    theme(legend.position = "top")
}
# plot_tradeoff_frontier(crash_res_inverse_dome_0.3)
# plot_tradeoff_frontier(d3)
# compare_time_blocks(s1)
# compare_time_blocks(s2)
# compare_time_blocks(s3)
