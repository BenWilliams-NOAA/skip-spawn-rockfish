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

plot_spr_curves <- function(rpt, data, skip_vec) {
  
  # 1. Extract Biological Inputs
  waa <- rpt$waa
  ages <- 1:length(waa)
  n_ages <- length(ages)
  M <- rpt$M
  sel <- rpt$slx[, 1]   # Fishery selectivity
  spawn_fract <- rpt$spawn_fract # e.g., 0.0 or 0.5 depending on spawning month
  
  # 2. Define Maturity Vectors
  # Scenario ii: Biological Maturity (The "All Mature" assumption)
  mat_bio <- data$maa 
  # Scenario i: Functional Maturity (The "Skipping" reality)
  mat_func <- mat_bio * (1 - skip_vec)
  
  # 3. Define Range of F to test
  f_range <- seq(0, 0.4, by = 0.001) # Increased range slightly to ensure capture of F35/F40
  
  calc_spr <- function(f_val, mat_vector) {
    
    # initialize
    N = numeric(n_ages)
    N[1] = 1 # Per recruit basis
    
    # --- Step A: Calculate Numbers at Age (Start of Year) ---
    # Matches: N_spr[a,1] = N_spr[a-1,1] * exp(-M) ... (but with Z)
    
    # 1. Loop through ages up to A-1
    for(a in 2:(n_ages)) {
      # Z for the previous age (a-1)
      Z_prev <- M + f_val * sel[a-1]
      N[a] <- N[a-1] * exp(-Z_prev)
    }
    
    # 2. Plus Group Correction
    # Matches: N_spr[A] = N_spr[A-1] * exp(-Z_prev) / (1 - exp(-Z_last))
    # Note: The loop above already calculated the numerator for N[n_ages] (N[A])
    # We just need to apply the geometric series division for the plus group
    Z_last <- M + f_val * sel[n_ages]
    N[n_ages] <- N[n_ages] / (1 - exp(-Z_last))
    
    # --- Step B: Calculate Spawning Biomass (Adjusted for timing) ---
    # Matches: sb_spr[a] = N_spr[a] * wt * exp(-fract * (M + F * sel))
    
    SB_a <- numeric(n_ages)
    
    for(a in 1:n_ages) {
      # Current total mortality for this age
      Z_curr <- M + f_val * sel[a]
      
      # Discount N to the time of spawning
      # Crucial change: This includes F in the exponent, matching EM code
      SB_a[a] <- N[a] * waa[a] * mat_vector[a] * exp(-spawn_fract * Z_curr)
    }
    
    # Return Sum (Total SPR)
    return(sum(SB_a))
  }
  
  # 5. Run Calculations
  results <- map_dfr(f_range, function(f) {
    spr_abs_i  <- calc_spr(f, mat_func)
    spr_abs_ii <- calc_spr(f, mat_bio)
    
    tibble::tibble(F_rate = f,
                   SPR_Absolute_i = spr_abs_i,
                   SPR_Absolute_ii = spr_abs_ii)
  })
  
  # 6. Calculate Relative SPR (SPR / SPR_0)
  # Get SPR at F=0 for normalization
  spr0_i  <- results$SPR_Absolute_i[1]
  spr0_ii <- results$SPR_Absolute_ii[1]
  
  plot_data <- results %>%
    mutate(
      # Normalize to get the % Ratio
      Ratio_i  = SPR_Absolute_i / spr0_i,
      Ratio_ii = SPR_Absolute_ii / spr0_ii
    ) %>%
    pivot_longer(cols = starts_with("Ratio"), names_to = "Scenario", values_to = "SPR_Ratio") %>%
    mutate(Scenario = ifelse(Scenario == "Ratio_i", "i (Functional)", "ii (Biological)"))
  
  # 7. Find the intersection with 0.4 (The F40 Reference Point)
  refs <- plot_data %>% 
    group_by(Scenario) %>% 
    filter(abs(SPR_Ratio - 0.5) == min(abs(SPR_Ratio - 0.5)))
  
  # 8. Plot
  ggplot(plot_data, aes(x = F_rate, y = SPR_Ratio, color = Scenario)) +
    geom_line(linewidth = 1.2) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "black") +
    geom_vline(data = refs, aes(xintercept = F_rate, color = Scenario), linetype="dotted") +
    geom_text(data = refs, aes(label = paste0("F40 = ", round(F_rate, 3)), y = 0.42), 
              show.legend = FALSE, hjust = -0.1) +
    scale_color_manual(values = c("i (Functional)" = "blue", "ii (Biological)" = "red")) +
    labs(
      title = "SPR Curves: Functional vs Biological Maturity",
      subtitle = "Matches EM Code: Includes F-mortality in spawning timing adjustment",
      x = "Fishing Mortality (F)",
      y = "SPR Ratio (Current / Unfished)"
    ) +
    theme_minimal()+
    expand_limits(y = 0) + 
    scale_y_continuous(breaks = seq(0, 1, 0.1))
}



p_skip %>% 
  left_join(data.frame(age=2:51, bio = rpt$maa)) %>% 
  mutate(dome = (1-dome) * bio,
         skewed_dome = (1-skewed_dome) * bio,
         increasing = (1-increasing) * bio,
         decreasing = (1-decreasing) * bio,
         inverse_dome = (1-inverse_dome) * bio,
         constant = (1-constant) * bio) %>% 
  pivot_longer(-c(age, skip, bio)) %>%
  ggplot(aes(age, value, color = factor(skip))) +
  geom_line() +
  facet_wrap(~name) +
  geom_line(aes(y=bio), color = 1)

skip_vector = p_skip %>% filter(skip == 0.3) %>% pull(decreasing)
plot_spr_curves(rpt = rpt, data = data, skip_vec = skip_vector) + 
  xlim(0, 0.15)







########################
scenarios <- p_skip %>%
  pivot_longer(cols = -c(age, skip), 
               names_to = "curve_type", 
               values_to = "skip_prob") %>%
  bind_rows(expand.grid(age = 1:50, 
                       skip = c(0, 0.3),
                       curve_type = "observed",
                       skip_prob = skip)) %>% 
  left_join(tidytable(age = 1:max_age, 
                      maa = rpt$maa, 
                      waa = rpt$waa, 
                      M = rpt$M,
                      sel = rpt$slx[,1])) %>%
  mutate(mat_eff = maa * (1 - skip_prob)) 

# 2. Run SPR Analysis over a range of Fishing Mortality (F)
spr_results <- scenarios %>%
  expand_grid(F_rate = seq(0, 0.25, by = 0.01)) %>%
  mutate(Z = M + (F_rate * sel)) %>% 
  mutate(
   N = exp(-cumsum(lag(Z, default = 0))),  # Cumulative survival
  ssb = N * waa * mat_eff * 0.5, # Spawning Biomass using the EFFECTIVE maturity (post-skip)
        .by = c(skip, curve_type, F_rate) ) %>% 
  summarise(
    ssb_total = sum(ssb),
    .by = c(skip, curve_type, F_rate)
  ) %>% 
  # Calculate Ratio relative to F=0 (Unfished) for THAT specific scenario
  mutate(
    ssb_virgin = first(ssb_total[F_rate == 0]),
    SPR = ssb_total / ssb_virgin,
    .by = c(skip, curve_type)
  ) 

spr_results %>%
  # filter(curve_type %in% c("constant")) %>%
  filter(skip %in% c(0,0.3,0.4)) %>% 
  mutate(skip = ifelse(skip==0, 0, 0.3),
         curve_type = factor(curve_type, levels = c("observed","decreasing", "inverse_dome", "constant", "dome", "increasing", "skewed_dome"))) %>% 
  ggplot(aes(F_rate, SPR, color = factor(curve_type), group = interaction(curve_type, factor(skip)))) + 
  geom_line() +
  geom_hline(yintercept = 0.4, linetype = "dashed", color = "gray50") +
  # scale_y_continuous(labels = scales::percent, limits = c(0.25,1)) +
  labs(
    title = "SPR curves under different skip spawning scenarios",
    x = "Fishing Mortality (F)",
    y = "Spawning Potential Ratio (SPR)",
    color = "Skip Level"
  ) +
  # facet_wrap(~curve_type) +
  afscassess::theme_report() +
  scico::scale_color_scico_d(palette = 'roma') +
  # scale_linetype_manual(values = c(1, 2, 2)) +
  xlim(0, 0.15) +
  ylim(0.25, 0.55)

p1 = spr_results %>%
  filter(curve_type %in% c("observed", "decreasing", "constant")) %>%
  filter(skip %in% c(0,0.3,0.4)) %>% 
  mutate(skip = ifelse(skip==0, 0, 0.3),
         curve_type = factor(curve_type, levels = c("observed","decreasing", "inverse_dome", "constant", "dome", "increasing", "skewed_dome"))) %>% 
  ggplot(aes(F_rate, SPR, color = factor(curve_type), group = interaction(curve_type, factor(skip)))) + 
  geom_line() +
  geom_hline(yintercept = 0.4, linetype = "dashed", color = "gray50") +
  # scale_y_continuous(labels = scales::percent, limits = c(0.25,1)) +
  labs(
    title = "SPR curves under different skip spawning scenarios",
    x = "Fishing Mortality (F)",
    y = "Spawning Potential Ratio (SPR)",
    color = "Skip Level"
  ) +
  # facet_wrap(~curve_type) +
  afscassess::theme_report() +
  scico::scale_color_scico_d(palette = 'roma') 


p2 = spr_results %>%
  # filter(curve_type %in% c("constant")) %>%
  filter(skip %in% c(0,0.3,0.4)) %>% 
  mutate(skip = ifelse(skip==0, 0, 0.3),
         curve_type = factor(curve_type, levels = c("observed","decreasing", "inverse_dome", "constant", "dome", "increasing", "skewed_dome"))) %>% 
  ggplot(aes(F_rate, SPR, color = factor(curve_type), group = interaction(curve_type, factor(skip)))) + 
  geom_line() +
  geom_hline(yintercept = 0.4, linetype = "dashed", color = "gray50") +
  # scale_y_continuous(labels = scales::percent, limits = c(0.25,1)) +
  labs(
    title = "SPR curves under different skip spawning scenarios",
    x = "Fishing Mortality (F)",
    y = "Spawning Potential Ratio (SPR)",
    color = "Skip Level"
  ) +
  # facet_wrap(~curve_type) +
  afscassess::theme_report() +
  scico::scale_color_scico_d(palette = 'roma')  +
  coord_cartesian(xlim = c(0.05, 0.15), ylim = c(0.3, 0.5)) +
  labs(x = NULL, y = NULL, title = NULL) + # Remove labels to save space
  theme(
    legend.position = "none",   # Hide legend in the zoom box
    plot.background = element_rect(fill = "white", color = "black"), # White background with black border
    axis.title = element_blank()
  )

p1 + 
  inset_element(p2, left = 0.5, bottom = 0.5, right = 0.99, top = 0.99)

library(ggforce)
p1 + 
  facet_zoom(xlim = c(0.05, 0.15), ylim = c(0.23, 0.5), horizontal = FALSE)
ggsave(here::here("figs", "spr.png"), width = 7.5, height = 5.5, dpi = 200)
