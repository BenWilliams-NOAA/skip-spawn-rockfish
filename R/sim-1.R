# notes ----
# average recruitment from recruit timeseries
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

plot(data$maa)
lines(skip)

# globals ----
n_iter <- 50
n_years <- 100 
bio_mat <- matrix(rep(data$maa, n_years), ncol = n_years)
func_mat = matrix(rep(skip, n_years), ncol = n_years)
# recruitment matrix
collapse = 0 # number of years of population collapse (if desired)
set.seed(309)
rec_matrix <- matrix(c(rlnorm(n_iter * collapse, 
                              quantile(log(rpt$recruits), 0.05), 0.1),
                       rlnorm(n_iter * (n_years - collapse), 
                              quantile(log(rpt$recruits), 0.50), 0.2)),
                     nrow = n_iter, ncol = n_years)

sim1 = omem_parallel(rpt, data, bio_mat, func_mat, rec_matrix, obj_f = f1)
s1 <- extract_results(sim1)
glimpse(s1)
# save results
saveRDS(sim1, here::here("output", "sim1.RDS"))
saveRDS(s1, here::here("output", "s1.RDS"))
s1 = readRDS(here::here("output", "s1.RDS"))


plot_risk(s1, ssb_col = "spawn_bio_r", bio_ref_col = "B35", F_col = "F35", rmv_yrs = 0) # mgmt view
plot_risk(s1, rmv_yrs = 0) # reality view
plot_risk_true(s1, rmv_yrs = 0) # reality view
plot_risk(s1, rmv_yrs = 0) # reality view
table_risk(s1, ssb_col = "spawn_bio_r", bio_ref_col = "B35", F_col = "F35", rmv_yrs = 0) %>% 
  select(-i) # mgmt view
table_risk(s1, rmv_yrs = 20)
plot_F_shock(rpt, s1)
plot_ssb_shock(rpt, s1) + scico::scale_color_scico_d(palette="roma")
plot_catch(rpt, s1) + scico::scale_color_scico_d(palette="roma")
plot_recruits(rpt, s1) + scico::scale_color_scico_d(palette="roma")
calc_bias(s1)

plot_true_reality(s1)
# Run the plot
plot_management_divergence(s1)
# Run the calculation
calc_risk_prob(s1, rmv_yrs = 0)

calc_risk_prob_corrected <- function(sims, rmv_yrs = 0) {
  require(dplyr)
  require(tidytable)
  
  # 1. Extract the "Gold Standard" Limits (from Scenario i only)
  true_limits <- sims %>%
    filter(scenario == "i") %>% 
    select(iteration, year, true_B40 = B40_f, true_B35 = B35_f)
  
  risk_table <- sims %>%
    filter(year > (min(year) + rmv_yrs)) %>%
    
    # 2. Join the True Limits to ALL rows
    left_join(true_limits, by = c("iteration", "year")) %>%
    
    mutate(
      # 3. Compare against the TRUE FIXED LIMITS (not the row-specific B40_f)
      breach_B40 = ifelse(spawn_bio_f < true_B40, 1, 0),
      breach_B35 = ifelse(spawn_bio_f < true_B35, 1, 0),
      
      # For B20, use the TRUE B35 as the base
      breach_B20 = ifelse(spawn_bio_f < (true_B35 * (20/35)), 1, 0)
    ) %>%
    summarise(
      prob_below_B40 = mean(breach_B40),
      prob_below_B35 = mean(breach_B35),
      prob_below_B20 = mean(breach_B20),
      .by = scenario
    ) %>%
    mutate(
      prob_below_B40 = round(prob_below_B40 * 100, 1),
      prob_below_B35 = round(prob_below_B35 * 100, 1),
      prob_below_B20 = round(prob_below_B20 * 100, 1)
    ) 
  
  return(risk_table)
}

calc_risk_prob_corrected(s1, rmv_yrs = 15)

# Create a small dataframe with your numbers
risk_data <- data.frame(
  Metric = rep(c("Below B40", "Below B35"), each = 2),
  Scenario = rep(c("i (Matched)", "ii (Mismatched)"), 2),
  Probability = c(91.3, 93.9, 26.3, 43.5)
)

ggplot(risk_data, aes(x = Metric, y = Probability, fill = Scenario)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("#2c7bb6", "#d7191c")) +
  geom_text(aes(label = paste0(Probability, "%")), 
            position = position_dodge(width = 0.7), 
            vjust = -0.5, fontface = "bold") +
  # Add an arrow or annotation for the B35 difference
  annotate("segment", x = 1.8, xend = 2.2, y = 46, yend = 46, 
           arrow = arrow(ends = "both", length = unit(0.2, "cm"))) +
  annotate("text", x = 2, y = 50, label = "+65% Relative Risk", color = "black", fontface = "italic") +
  labs(
    title = "Impact of Mismatch on Biological Risk",
    subtitle = "The mismatch nearly doubles the frequency of dropping below B35",
    y = "Probability of Occurrence (%)"
  ) +
  theme_minimal()
# Run the plot
plot_risk_worms(s1)


ggplot(s1, aes(x = year, y = spawn_bio_f, color = scenario)) +
  geom_line() + 
  scale_color_manual(values = c("i" = "#2c7bb6", "ii" = "#d7191c")) +
  labs(
    title = "Stock Health: Functional Spawning Biomass",
    subtitle = "Comparing Matched (i) vs. Mismatched (ii) Management",
    y = "True Spawning Biomass (Functional)",
    x = "Year"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

# 1. Pivot the data to get 'i' and 'ii' side-by-side
bias_analysis <- s1 %>%
  select(iteration, year, scenario, spawn_bio_f) %>%
  pivot_wider(
    names_from = scenario, 
    values_from = spawn_bio_f, 
    names_prefix = "ssb_"
  ) %>%
  # Calculate the Ratio: (Mismatched / Matched)
  mutate(relative_biomass = ssb_ii / ssb_i)

# 2. Plot the Ratio over time
ggplot(bias_analysis, aes(x = year, y = relative_biomass)) +
  # Add reference line at 1.0 (No Bias / Perfect Match)
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.8) +
  
  # Plot the median trajectory with a confidence ribbon
  stat_summary(fun = median, geom = "line", color = "#d7191c", linewidth = 1.2) +
  stat_summary(fun.data = median_hilow, geom = "ribbon", fill = "#d7191c", alpha = 0.2) +
  
  # Formatting
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "The Cost of Mismatch: Biomass Depletion",
    subtitle = "True Biomass of Scenario 'ii' as a percentage of Scenario 'i'",
    y = "Relative Biomass (Scenario ii / Scenario i)",
    x = "Projection Year",
    caption = "Values < 100% indicate the mismatched management caused a stock decline."
  ) +
  theme_minimal(base_size = 14)


# 1. Reshape your 's1' object into long format for plotting
plot_data <- s1 %>% 
  # Select only what we need
  select(year, scenario, spawn_bio_f, catch, F_target_val, iteration) %>%
  # Pivot to long format so we can facet
  pivot_longer(
    cols = c(spawn_bio_f, catch, F_target_val), 
    names_to = "Metric", 
    values_to = "Value"
  ) %>%
  # Rename metrics for pretty facet labels
  mutate(Metric = factor(Metric, 
                         levels = c("spawn_bio_f", "catch", "F_target_val"),
                         labels = c("True Spawning Biomass", "Catch", "Target F")
  ))

# 2. Create the faceted plot
ggplot(plot_data, aes(x = year, y = Value, color = scenario, fill = scenario)) +
  # stat_summary(fun.data = mean_se, geom = "ribbon", alpha = .2, color=NA) +
  stat_summary(fun = median, geom = "line") +
  facet_wrap(~Metric, scales = "free_y", ncol = 1) +
  scale_color_manual(values = c("i" = "#2c7bb6", "ii" = "#d7191c")) +
  labs(
    title = "MSE Simulation Results",
    y = "Value",
    x = "Year"
  ) +
  theme_bw() +
  theme(legend.position = "top")

# Summarize data first
s1_summary <- s1 %>%
  summarise(
    median_ssb = median(spawn_bio_f, na.rm = TRUE),
    lower_ssb = quantile(spawn_bio_f, 0.025, na.rm = TRUE),
    upper_ssb = quantile(spawn_bio_f, 0.975, na.rm = TRUE),
    .by = c(year, scenario)
  )

# Plot Median with Ribbon
ggplot(s1_summary, aes(x = year, y = median_ssb, color = scenario, fill = scenario)) +
  geom_ribbon(aes(ymin = lower_ssb, ymax = upper_ssb), alpha = 0.2, color = NA) +
  geom_line() +
  scale_color_manual(values = c("i" = "#2c7bb6", "ii" = "#d7191c")) +
  scale_fill_manual(values = c("i" = "#2c7bb6", "ii" = "#d7191c")) +
  labs(
    title = "Projected Stock Status (Median + 95% CI)",
    y = "True Functional SSB",
    x = "Year"
  ) +
  theme_minimal()

# Create a "Wide" dataset for direct math
comparison_df <- s1 %>%
  select(iteration, year, scenario, spawn_bio_f, catch, F_target_val) %>%
  pivot_wider(
    names_from = scenario, 
    values_from = c(spawn_bio_f, catch, F_target_val),
    names_glue = "{.value}_{scenario}" # Creates spawn_bio_f_i, spawn_bio_f_ii, etc.
  ) %>%
  mutate(
    # 1. Biological Impact: How much smaller is the stock in ii?
    ssb_ratio = spawn_bio_f_ii / spawn_bio_f_i,
    ssb_diff_pct = (spawn_bio_f_ii - spawn_bio_f_i) / spawn_bio_f_i * 100,
    # bb40_ii = as.numeric(spawn_bio_f_ii < B40),
    # bb40_i = spawn_bio_f_i < B40,
    # 2. Economic Impact: Did we catch more or less?
  )

# View the first few rows
glimpse(comparison_df)

kpi_summary <- comparison_df %>%
  mutate(period = case_when(
    year <= 10 ~ "Short Term (1-10)",
    year >= (max(year) - 10) ~ "Long Term (End)",
    TRUE ~ "Mid Term"
  )) %>% 
  group_by(period) %>%
  summarise(
    # Median Biological Loss: (e.g., -5% means ii is 5% lower than i)
    median_ssb_impact = median(ssb_diff_pct, na.rm = TRUE),
    # below_b40 = sum(bb40_ii, na.rm = T) / n(),
    # Probability of being worse: How often is ii < i?
    risk_of_loss = mean(spawn_bio_f_ii < spawn_bio_f_i) * 100,
    
    # Cumulative Catch Difference (Total tons gained/lost)
    avg_catch_diff_tons = mean(catch_ii - catch_i),
    
    # Intensity Bias: How much harder are we fishing?
    median_F_bias = median((F_target_val_ii - F_target_val_i)/F_target_val_i * 100, na.rm=TRUE))

print(kpi_summary)

truth_ref <- master_df %>%
  filter(scenario == "i") %>% 
  select(iteration, year, shape, skip_level, F35_true = F40_f)


# 1. Extract the "Gold Standard" Limits (from Scenario i)
# These are the actual biological safety lines we don't want to cross.
true_limits <- s1 %>%
  filter(scenario == "i") %>% # Ensure this matches your label (e.g., "Matched (i)")
  # Note: If your s1 object has 'shape' and 'skip_level' columns, include them in select()!
  select(iteration, year, true_B40 = B40_f, true_B35 = B35_f) 

# 2. Create the Wide Dataset with Risk Flags
comparison_df <- s1 %>%
  select(iteration, year, scenario, spawn_bio_f, catch, F_target_val) %>%
  pivot_wider(
    names_from = scenario, 
    values_from = c(spawn_bio_f, catch, F_target_val),
    names_glue = "{.value}_{scenario}"
  ) %>%
  
  # Join the "True Limits" so every row has the correct ruler
  # Note: Add 'shape' and 'skip_level' to 'by' if they exist in your data
  left_join(true_limits, by = c("iteration", "year")) %>%
  
  mutate(
    # --- Existing Metrics ---
    ssb_ratio = spawn_bio_f_ii / spawn_bio_f_i,
    ssb_diff_pct = (spawn_bio_f_ii - spawn_bio_f_i) / spawn_bio_f_i * 100,
    
    # --- NEW: Risk Flags (Boolean TRUE/FALSE) ---
    # We compare BOTH scenarios against the 'true_B40/B35'
    
    # Did Scenario i fail?
    below_B40_i = spawn_bio_f_i < true_B40,
    below_B35_i = spawn_bio_f_i < true_B35,
    
    # Did Scenario ii fail? (Tested against the True Limit, not its biased one)
    below_B40_ii = spawn_bio_f_ii < true_B40,
    below_B35_ii = spawn_bio_f_ii < true_B35
  )

# 3. (Optional) Quick Summary of Risk
# Calculate the probability of dropping below limits
risk_summary <- comparison_df %>%
  summarise(
    prob_B40_i  = mean(below_B40_i, na.rm = TRUE),
    prob_B40_ii = mean(below_B40_ii, na.rm = TRUE),
    prob_B35_i  = mean(below_B35_i, na.rm = TRUE),
    prob_B35_ii = mean(below_B35_ii, na.rm = TRUE)
  )

print(risk_summary)

# 2. Join it back using ALL matching keys
intensity_data_corrected <- master_df %>%
  left_join(truth_ref, by = c("iteration", "year", "shape", "skip_level")) %>%
  mutate(
    # Now we are certain we are dividing by the correct limit for this specific run
    realized_intensity = F_target_val / F35_true
  )

# 3. Plot - Now you should see the separation
ggplot(intensity_data_corrected, aes(x = year, y = realized_intensity, color = scenario, group = interaction(scenario, iteration))) +
  # stat_summary(fun = median, geom = "line", size = 1.2) +
  geom_line() +
  scale_color_manual(values = c("i" = "#2c7bb6", "ii" = "#d7191c")) +
  labs(
    title = "Realized Fishing Intensity (Fully Corrected)",
    subtitle = "Ratio of Target F to the Scenario-Specific True F35 Limit",
    y = "True Intensity (F_target / F35_true)"
  ) +
  theme_minimal()

# 3. Plot the corrected metric
ggplot(intensity_data_corrected, aes(x = year, y = realized_intensity, color = scenario)) +
  stat_summary(fun = median, geom = "line", size = 1.2) +
  labs(
    title = "Realized Fishing Intensity (Corrected)",
    subtitle = "Ratio of Target F to the TRUE Functional F35 Limit",
    y = "True Intensity (F_target / F35_true)"
  ) +
  theme_minimal()
