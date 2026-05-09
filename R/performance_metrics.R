# performance metrics

# load ----
library(tidytable)
library(ggplot2)
theme_set(afscassess::theme_report())

# functions ----
rel_dif <- function(est, obs) {
  dplyr::if_else(obs == 0, NA_real_, (est / obs) - 1)
}
# all performance metrics
calc_all_pm <- function(sims, rmv_yrs = 0) {
  
  sims %>%
    # remove any burn-in years 
    filter(year >= (min(year) + rmv_yrs)) %>%
    arrange(scenario, iteration, year) %>%
    # base metrics
    summarise(
      avg_catch      = mean(catch, na.rm = TRUE),
      avg_F          = mean(Ft, na.rm = TRUE),
      avg_SSB        = mean(spawn_bio_f, na.rm = TRUE),
      prob_F_excess  = mean(Ft > F35_f, na.rm = TRUE),
      prob_SSB_low   = mean(spawn_bio_f < B35_f, na.rm = TRUE),
      term_OFG       = last(Ft) / last(F35_f),
      term_OFD       = last(spawn_bio_f) / last(B35_f),
      aacv           = sum(abs(diff(catch)), na.rm = TRUE) / sum(head(catch, -1), na.rm = TRUE),
      .by = c(scenario, iteration)) %>%
    pivot_wider(
      names_from = scenario, 
      values_from = c(
        avg_catch, avg_F, avg_SSB, prob_F_excess, 
        prob_SSB_low, term_OFG, term_OFD, aacv
      )) %>%
    # relative differences (scenario ii relative to scenario i)
    mutate(
      rd_avg_catch     = rel_dif(avg_catch_ii, avg_catch_i),
      rd_avg_F         = rel_dif(avg_F_ii, avg_F_i),
      rd_avg_SSB       = rel_dif(avg_SSB_ii, avg_SSB_i),
      rd_prob_F_excess = rel_dif(prob_F_excess_ii, prob_F_excess_i),
      rd_prob_SSB_low  = rel_dif(prob_SSB_low_ii, prob_SSB_low_i),
      rd_term_OFG      = rel_dif(term_OFG_ii, term_OFG_i),
      rd_term_OFD      = rel_dif(term_OFD_ii, term_OFD_i),
      rd_aacv          = rel_dif(aacv_ii, aacv_i)
    )
}

calc_aacv <- function(sims, rmv_yrs = 0) {
  sims %>% 
  filter(year >= (min(year) + rmv_yrs)) %>%
    arrange(scenario, iteration, year) %>%
    mutate(
      catch_prev = lag(catch), # NOTE: Change 'catch' if your column is named differently (e.g., catch_f)
      catch_diff = abs(catch - catch_prev),
      .by = c(scenario, iteration)
    ) %>% 
    summarize(
      aacv = sum(catch_diff, na.rm = TRUE) / sum(catch_prev, na.rm = TRUE) * 100,
      avg_catch = mean(catch, na.rm = TRUE),
      .by = c(scenario, iteration)
    ) #%>% 
    # pivot_wider(
    #   names_from = scenario, 
    #   values_from = c(aacv, avg_catch)
    # )
}

# calc_aacv(s1)
# data ----

s1 = readRDS(here::here("output", "s1.rds"))
s2 = readRDS(here::here("output", "s2.rds"))
s3 = readRDS(here::here("output", "s3.rds"))
s4 = readRDS(here::here("output", "s5.rds"))

s5 = readRDS(here::here("output", "crash_res_constant_0.3.RDS"))
s6 = readRDS(here::here("output", "crash_res_increasing_0.3.RDS"))
s7 = readRDS(here::here("output", "crash_res_decreasing_0.3.RDS"))
s8 = readRDS(here::here("output", "crash_res_dome_0.3.RDS"))
s9 = readRDS(here::here("output", "crash_res_skewed_dome_0.3.RDS"))
s10 = readRDS(here::here("output", "crash_res_inverse_dome_0.3.RDS"))

s5_2 = readRDS(here::here("output", "crash_res_constant_0.2.RDS"))
s6_2 = readRDS(here::here("output", "crash_res_increasing_0.2.RDS"))
s7_2 = readRDS(here::here("output", "crash_res_decreasing_0.2.RDS"))
s8_2 = readRDS(here::here("output", "crash_res_dome_0.2.RDS"))
s9_2 = readRDS(here::here("output", "crash_res_skewed_dome_0.2.RDS"))
s10_2 = readRDS(here::here("output", "crash_res_inverse_dome_0.2.RDS"))

s5_1 = readRDS(here::here("output", "crash_res_constant_0.1.RDS"))
s6_1 = readRDS(here::here("output", "crash_res_increasing_0.1.RDS"))
s7_1 = readRDS(here::here("output", "crash_res_decreasing_0.1.RDS"))
s8_1 = readRDS(here::here("output", "crash_res_dome_0.1.RDS"))
s9_1 = readRDS(here::here("output", "crash_res_skewed_dome_0.1.RDS"))
s10_1 = readRDS(here::here("output", "crash_res_inverse_dome_0.1.RDS"))

s11 = readRDS(here::here("output", "high_res_constant_0.3.RDS"))
s12 = readRDS(here::here("output", "high_res_increasing_0.3.RDS"))
s13 = readRDS(here::here("output", "high_res_decreasing_0.3.RDS"))
s14 = readRDS(here::here("output", "high_res_dome_0.3.RDS"))
s15 = readRDS(here::here("output", "high_res_skewed_dome_0.3.RDS"))
s16 = readRDS(here::here("output", "high_res_inverse_dome_0.3.RDS"))

s11_2 = readRDS(here::here("output", "high_res_constant_0.2.RDS"))
s12_2 = readRDS(here::here("output", "high_res_increasing_0.2.RDS"))
s13_2 = readRDS(here::here("output", "high_res_decreasing_0.2.RDS"))
s14_2 = readRDS(here::here("output", "high_res_dome_0.2.RDS"))
s15_2 = readRDS(here::here("output", "high_res_skewed_dome_0.2.RDS"))
s16_2 = readRDS(here::here("output", "high_res_inverse_dome_0.2.RDS"))

s11_1 = readRDS(here::here("output", "high_res_constant_0.1.RDS"))
s12_1 = readRDS(here::here("output", "high_res_increasing_0.1.RDS"))
s13_1 = readRDS(here::here("output", "high_res_decreasing_0.1.RDS"))
s14_1 = readRDS(here::here("output", "high_res_dome_0.1.RDS"))
s15_1 = readRDS(here::here("output", "high_res_skewed_dome_0.1.RDS"))
s16_1 = readRDS(here::here("output", "high_res_inverse_dome_0.1.RDS"))


s17 = readRDS(here::here("output", "avg_res_constant_0.3.RDS"))
s18 = readRDS(here::here("output", "avg_res_increasing_0.3.RDS"))
s19 = readRDS(here::here("output", "avg_res_decreasing_0.3.RDS"))
s20 = readRDS(here::here("output", "avg_res_dome_0.3.RDS"))
s21 = readRDS(here::here("output", "avg_res_skewed_dome_0.3.RDS"))
s22 = readRDS(here::here("output", "avg_res_inverse_dome_0.3.RDS"))

s17_2 = readRDS(here::here("output", "avg_res_constant_0.2.RDS"))
s18_2 = readRDS(here::here("output", "avg_res_increasing_0.2.RDS"))
s19_2 = readRDS(here::here("output", "avg_res_decreasing_0.2.RDS"))
s20_2 = readRDS(here::here("output", "avg_res_dome_0.2.RDS"))
s21_2 = readRDS(here::here("output", "avg_res_skewed_dome_0.2.RDS"))
s22_2 = readRDS(here::here("output", "avg_res_inverse_dome_0.2.RDS"))

s17_1 = readRDS(here::here("output", "avg_res_constant_0.1.RDS"))
s18_1 = readRDS(here::here("output", "avg_res_increasing_0.1.RDS"))
s19_1 = readRDS(here::here("output", "avg_res_decreasing_0.1.RDS"))
s20_1 = readRDS(here::here("output", "avg_res_dome_0.1.RDS"))
s21_1 = readRDS(here::here("output", "avg_res_skewed_dome_0.1.RDS"))
s22_1 = readRDS(here::here("output", "avg_res_inverse_dome_0.1.RDS"))

# crash recruitment
df_crash <- list(
  "0.3" = list("constant" = s5, "increasing" = s6, "decreasing" = s7,
    "dome" = s8, "skewed dome" = s9, "inverse dome" = s10),
  "0.2" = list("constant" = s5_2, "increasing" = s6_2, "decreasing" = s7_2,
    "dome" = s8_2, "skewed dome" = s9_2, "inverse dome" = s10_2),
  "0.1" = list("constant" = s5_1, "increasing" = s6_1, "decreasing" = s7_1,
    "dome" = s8_1, "skewed dome" = s9_1, "inverse dome" = s10_1 )
) %>% 
  # outer map: iterates over "0.3", "0.2", "0.1" and creates the 'skip' column
  map_dfr(function(inner_list) {
    # inner map: applies calc_bias2 to each dataset and creates the 'id' column
    map_dfr(inner_list, calc_all_pm, .id = "skip_shape")
  }, .id = "skip") %>% 
  mutate(recruitment = "crash")


# high recruitment
df_high <- list(
  "0.3" = list("constant" = s11, "increasing" = s12, "decreasing" = s13, "dome" = s14, "skewed dome" = s15, "inverse dome" = s16),
  "0.2" = list("constant" = s11_2, "increasing" = s12_2, "decreasing" = s13_2, "dome" = s14_2, "skewed dome" = s15_2, "inverse dome" = s16_2),
  "0.1" = list("constant" = s11_1, "increasing" = s12_1, "decreasing" = s13_1, "dome" = s14_1, "skewed dome" = s15_1, "inverse dome" = s16_1)
) %>% 
  map_dfr(function(inner_list) {
    map_dfr(inner_list, calc_all_pm, .id = "skip_shape")
  }, .id = "skip") %>% 
  mutate(recruitment = "high")

# 2. Average Recruitment
df_avg <- list(
  "0.3" = list("constant" = s17, "increasing" = s18, "decreasing" = s19, "dome" = s20, "skewed dome" = s21, "inverse dome" = s22),
   "0.2" = list("constant" = s17_2, "increasing" = s18_2, "decreasing" = s19_2, "dome" = s20_2, "skewed dome" = s21_2, "inverse dome" = s22_2),
    "0.1" = list("constant" = s17_1, "increasing" = s18_1, "decreasing" = s19_1, "dome" = s20_1, "skewed dome" = s21_1, "inverse dome" = s22_1)
) %>% 
  map_dfr(~ map_dfr(.x, calc_all_pm, .id = "skip_shape"), .id = "skip") %>% 
  mutate(recruitment = "average")

# 3. Observed Data (No skip levels, so we just map once)
df_observed <- list(
  "average" = s1, "high" = s2, "crash" = s3) %>%
  map_dfr(calc_all_pm, .id = "recruitment") %>%
  mutate(
    skip_shape = "observed",
    skip = "obs" # filler so the column exists when binding
  )


all_pm <- bind_rows(df_crash, df_high, df_avg, df_observed) %>%
  mutate(
    recruitment = factor(recruitment, levels = c("average", "high", "crash")),
    skip_shape = factor(skip_shape, levels = c("observed", "constant", "increasing", "decreasing", "dome", "skewed dome", "inverse dome")),
    # Optional: factorize skip to control legend order
    skip = factor(skip, levels = c("0.1", "0.2", "0.3", "obs")) 
  )

# tables ----

# rd_avg_catch     = rel_dif(avg_catch_ii, avg_catch_i),
#       rd_avg_F         = rel_dif(avg_F_ii, avg_F_i),
#       rd_avg_SSB       = rel_dif(avg_SSB_ii, avg_SSB_i),
#       rd_prob_F_excess = rel_dif(prob_F_excess_ii, prob_F_excess_i),
#       rd_prob_SSB_low  = rel_dif(prob_SSB_low_ii, prob_SSB_low_i),
#       rd_term_OFG      = rel_dif(term_OFG_ii, term_OFG_i),
#       rd_term_OFD      = rel_dif(term_OFD_ii, term_OFD_i),
#       rd_aacv          = rel_dif(aacv_ii, aacv_i)

all_pm %>% 
  ggplot(aes(skip_shape, rd_term_OFD, fill = skip)) + 
  geom_hline(yintercept = 0, lty = 3, color = "gray30") +
  geom_boxplot(position = position_dodge2(width = 0.8, preserve = "single"), 
               alpha = 0.8, outlier.size = 1, outlier.alpha = 0.5) +
  facet_grid(recruitment ~ .) + 
  scico::scale_fill_scico_d(palette = "berlin", begin = 0.2, end = 0.8) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(face = "bold", size = 11)
  ) +
  labs(
    title = "Impact of Skip Spawning on Terminal Overfished Status",
    subtitle = "Relative difference between misspecified and baseline EM after 100 years",
    x = "Skip Spawning Shape",
    y = "Relative Difference in Terminal SSB / B35",
    fill = "Skip Level"
  )

all_pm %>%
  summarise(
    med_OFD = median(rd_term_OFD, na.rm = TRUE),
    med_OFG = median(rd_term_OFG, na.rm = TRUE),
    .by = c(recruitment, skip_shape, skip)
  ) -> temp
  
  all_pm %>% 
  ggplot(aes(x = rd_term_OFD, y = rd_term_OFG, color = skip)) +
  geom_vline(xintercept = 0, lty = 3, color = "gray50", linewidth = 0.5) +
  geom_hline(yintercept = 0, lty = 3, color = "gray50", linewidth = 0.5) +
  geom_point(alpha = 0.3) +
  geom_point(data = temp,
    aes(x=med_OFD, y= med_OFG, fill = skip),
    shape = 21,       # shape 21 allows both fill and color (outline)
    color = "black",  # black outline
    size = 2.5, 
    stroke = 0.8 ) +
   facet_grid(recruitment ~ skip_shape) +
  scico::scale_color_scico_d(palette = "berlin", begin = 0.2, end = 0.8) +
  scico::scale_fill_scico_d(palette = "berlin", begin = 0.2, end = 0.8) +
  theme( legend.position = "bottom" ) +
  labs(
    title = "Terminal Stock Status (Kobe Plot) by Skip Spawning Dynamics",
    subtitle = "Solid points represent scenario medians; transparent points represent simulation iterations",
    x = expression("Terminal SSB / B"[35]*" (term_OFD)"),
    y = expression("Terminal F / F"[35]*" (term_OFG)"),
    color = "Skip Level",
    fill  = "Skip Level"
  ) +
    scale_x_continuous(breaks = c(-.10, 0, .1))
    coord_cartesian(
    xlim = c(0.5, max(all_pm$term_OFD_ii, na.rm = TRUE) + 0.1), 
    ylim = c(0.5, 1.3) # Adjust this upper limit based on typical Kobe plot standards
  )

all_pm %>%
  summarise(
    # Median Catch with 95% Quantile Interval (rounded to whole numbers)
    Catch = sprintf("%.0f (%.0f - %.0f)", 
                    median(avg_catch_ii, na.rm = TRUE), 
                    quantile(avg_catch_ii, 0.025, na.rm = TRUE), 
                    quantile(avg_catch_ii, 0.975, na.rm = TRUE)),
    
    # Mean probability of dropping below B35 (formatted as a percentage)
    Prob_SSB_Low = sprintf("%.1f%%", mean(prob_SSB_low_ii, na.rm = TRUE) * 100),
    
    # Median Terminal OFD with 95% Quantile Interval (rounded to 2 decimals)
    Term_OFD = sprintf("%.2f (%.2f - %.2f)", 
                       median(term_OFD_ii, na.rm = TRUE), 
                       quantile(term_OFD_ii, 0.025, na.rm = TRUE), 
                       quantile(term_OFD_ii, 0.975, na.rm = TRUE)),
    
    .by = c(recruitment, skip_shape, skip)
  ) %>%
  # Sort for logical reading
  arrange(recruitment, skip_shape, skip)

all_pm %>%
  # filter(skip_shape != "observed") %>% 
  summarise(
    # Median Relative Difference in Catch (%)
    RD_Catch_Med = round(median(rd_avg_catch, na.rm = TRUE) * 100, 1),
    
    # Median Relative Difference in Terminal OFD (%)
    RD_Term_OFD_Med = round(median(rd_term_OFD, na.rm = TRUE) * 100, 1),
    
    # Risk: What proportion of iterations resulted in WORSE biomass than the baseline?
    # (i.e., rd_term_OFD < 0)
    Prob_Worse_SSB = round(mean(rd_term_OFD < 0, na.rm = TRUE) * 100, 1),
    
    .by = c(recruitment, skip_shape, skip)
  ) %>%
  arrange(recruitment, skip_shape, skip) %>% 
  # Clean up the output names for readability
  rename_with(~gsub("_ii", "", .x)) %>%
  arrange(recruitment, skip_shape, skip) %>% 
  pivot_wider(
    names_from = skip,
    values_from = c(RD_Catch_Med, RD_Term_OFD_Med, Prob_Worse_SSB)
  ) %>%
  # Reorder columns so that all 0.1 metrics are grouped together, then 0.2, etc.
  select(
    recruitment, 
    skip_shape, 
    ends_with("0.1"), 
    ends_with("0.2"), 
    ends_with("0.3")
  ) %>%
  arrange(recruitment, skip_shape)%>%
  kbl(
    format = "html",
    booktabs = TRUE,
    linesep = "", 
    col.names = c("Recruitment", "Skip Shape", 
                  "Catch", "SSB", "Risk", 
                  "Catch", "SSB", "Risk", 
                  "Catch", "SSB", "Risk"),
    caption = "Median relative difference in catch and terminal SSB, and the probability of worse terminal biomass (Risk) compared to the correctly specified baseline."
  ) %>%
  add_header_above(c(" " = 2, 
                     "Skip Level = 0.1" = 3, 
                     "Skip Level = 0.2" = 3, 
                     "Skip Level = 0.3" = 3)) %>%
  collapse_rows(columns = 1, valign = "top", latex_hline = "major") %>%
  kable_styling(latex_options = c("scale_down", "hold_position"))

all_pm %>%
  summarise(
    across(
      .cols = c(avg_catch_ii, avg_F_ii, avg_SSB_ii, aacv_ii, term_OFD_ii, term_OFG_ii),
      .fns = list(
        Median = ~median(.x, na.rm = TRUE),
        Mean   = ~mean(.x, na.rm = TRUE),
        SD     = ~sd(.x, na.rm = TRUE),
        Q025   = ~quantile(.x, 0.025, na.rm = TRUE),
        Q975   = ~quantile(.x, 0.975, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .by = c(recruitment, skip_shape, skip)
  )

all_pm %>% 
  ggplot(aes(skip_shape, rd_term_OFD/rd_term_OFG, fill = skip)) + 
  geom_hline(yintercept = 0, lty = 3, color = "gray30") +
  geom_boxplot(position = position_dodge2(width = 0.8, preserve = "single"), 
               alpha = 0.8, outlier.size = 1, outlier.alpha = 0.5) +
  facet_grid(recruitment ~ .) + 
  scico::scale_fill_scico_d(palette = "berlin", begin = 0.2, end = 0.8) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(face = "bold", size = 11)
  ) +
  labs(
    title = "Impact of Skip Spawning on Terminal Overfished Status",
    subtitle = "Relative difference between misspecified and baseline EM after 100 years",
    x = "Skip Spawning Shape",
    y = "Relative Difference in Terminal SSB / B35",
    fill = "Skip Level"
  ) +
  coord_cartesian(ylim = c(-25, 25))


all_pm %>% 
  # replace_na(list(rd_prob_F_excess =0)) %>% 
  ggplot(aes(skip_shape, rd_prob_F_excess, fill = skip)) + 
  geom_hline(yintercept = 0, lty = 3, color = "gray30") +
  geom_boxplot(position = position_dodge2(width = 0.8, preserve = "single"), 
               alpha = 0.8, outlier.size = 1, outlier.alpha = 0.5) +
  facet_grid(recruitment ~ ., scales = "free_y") + 
  scico::scale_fill_scico_d(palette = "berlin", begin = 0.2, end = 0.8) +

  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(face = "bold", size = 11)
  ) +
  labs(
    title = "Impact of Skip Spawning on Terminal Overfished Status",
    subtitle = "Relative difference between misspecified and baseline EM after 100 years",
    x = "Skip Spawning Shape",
    y = "Relative Difference in Terminal SSB / B35",
    fill = "Skip Level"
  )


all_pm %>% 
  # replace_na(list(rd_prob_F_excess =0)) %>% 
  ggplot(aes(skip_shape, rd_aacv, fill = skip)) + 
  geom_hline(yintercept = 0, lty = 3, color = "gray30") +
  geom_boxplot(position = position_dodge2(width = 0.8, preserve = "single"), 
               alpha = 0.8, outlier.size = 1, outlier.alpha = 0.5) +
  facet_grid(recruitment ~ .) + 
  scico::scale_fill_scico_d(palette = "berlin", begin = 0.2, end = 0.8) +

  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(face = "bold", size = 11)
  ) +
  labs(
    title = "Impact of Skip Spawning on Terminal Overfished Status",
    subtitle = "Relative difference between misspecified and baseline EM after 100 years",
    x = "Skip Spawning Shape",
    y = "Relative Difference in Terminal SSB / B35",
    fill = "Skip Level"
  )


all_pm %>% 
  ggplot(aes(skip_shape, rd_prob_SSB_low , fill = skip)) + 
  geom_hline(yintercept = 0, lty = 3, color = "gray30") +
  geom_boxplot(position = position_dodge2(width = 0.8, preserve = "single"), 
               alpha = 0.8, outlier.size = 1, outlier.alpha = 0.5) +
  facet_grid(recruitment ~ ., scales = "free_y") + 
  scico::scale_fill_scico_d(palette = "berlin", begin = 0.2, end = 0.8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    panel.grid.major.x = element_blank(),
  ) +
  labs(title = "Probability that SSB < B35",
    x = "Skip Spawning Shape",
    y = "SSB < B35",
    fill = "Skip Level"
  )


all_pm %>% 
  select(skip, skip_shape, iteration, avg_catch_i, avg_catch_ii, recruitment) %>% 
   mutate(aci = diff(avg_catch_i), .by = c(skip, skip_shape, recruitment, iteration))



# catch series ----
# crash recruitment
df_crash <- list(
  "0.3" = list("constant" = s5, "increasing" = s6, "decreasing" = s7,
    "dome" = s8, "skewed dome" = s9, "inverse dome" = s10),
  "0.2" = list("constant" = s5_2, "increasing" = s6_2, "decreasing" = s7_2,
    "dome" = s8_2, "skewed dome" = s9_2, "inverse dome" = s10_2),
  "0.1" = list("constant" = s5_1, "increasing" = s6_1, "decreasing" = s7_1,
    "dome" = s8_1, "skewed dome" = s9_1, "inverse dome" = s10_1 )
) %>% 
  # outer map: iterates over "0.3", "0.2", "0.1" and creates the 'skip' column
  map_dfr(function(inner_list) {
    # inner map: applies calc_bias2 to each dataset and creates the 'id' column
    map_dfr(inner_list, calc_aacv, .id = "skip_shape")
  }, .id = "skip") %>% 
  mutate(recruitment = "crash")


# high recruitment
df_high <- list(
  "0.3" = list("constant" = s11, "increasing" = s12, "decreasing" = s13, "dome" = s14, "skewed dome" = s15, "inverse dome" = s16),
  "0.2" = list("constant" = s11_2, "increasing" = s12_2, "decreasing" = s13_2, "dome" = s14_2, "skewed dome" = s15_2, "inverse dome" = s16_2),
  "0.1" = list("constant" = s11_1, "increasing" = s12_1, "decreasing" = s13_1, "dome" = s14_1, "skewed dome" = s15_1, "inverse dome" = s16_1)
) %>% 
  map_dfr(function(inner_list) {
    map_dfr(inner_list, calc_aacv, .id = "skip_shape")
  }, .id = "skip") %>% 
  mutate(recruitment = "high")

# 2. Average Recruitment
df_avg <- list(
  "0.3" = list("constant" = s17, "increasing" = s18, "decreasing" = s19, "dome" = s20, "skewed dome" = s21, "inverse dome" = s22),
   "0.2" = list("constant" = s17_2, "increasing" = s18_2, "decreasing" = s19_2, "dome" = s20_2, "skewed dome" = s21_2, "inverse dome" = s22_2),
    "0.1" = list("constant" = s17_1, "increasing" = s18_1, "decreasing" = s19_1, "dome" = s20_1, "skewed dome" = s21_1, "inverse dome" = s22_1)
) %>% 
  map_dfr(~ map_dfr(.x, calc_aacv, .id = "skip_shape"), .id = "skip") %>% 
  mutate(recruitment = "average")

# 3. Observed Data (No skip levels, so we just map once)
df_observed <- list(
  "average" = s1, "high" = s2, "crash" = s3
) %>%
  map_dfr(calc_aacv, .id = "recruitment") %>%
  mutate(
    skip_shape = "observed",
    skip = "obs" # filler so the column exists when binding
  )


all_catch <- bind_rows(df_crash, df_high, df_avg, df_observed) %>%
  mutate(
    recruitment = factor(recruitment, levels = c("average", "high", "crash")),
    skip_shape = factor(skip_shape, levels = c("observed", "constant", "increasing", "decreasing", "dome", "skewed dome", "inverse dome")),
    # Optional: factorize skip to control legend order
    skip = factor(skip, levels = c("0.1", "0.2", "0.3", "obs")) 
  )


all_catch %>% 
	ggplot(aes(skip_shape, aacv, fill = skip, color = scenario)) +
    geom_violin(
    trim = FALSE,           
    position = position_dodge(width = 0.8), 
    alpha = 0.7) +
    facet_wrap(~ recruitment, scales = "free_y", ncol = 1) +
  scale_y_continuous(labels = scales::comma) +
  scico::scale_fill_scico_d(palette = "berlin") +
  scico::scale_color_scico_d(palette='roma')

all_catch %>% drop_na() %>% 
	ggplot(aes(skip_shape, avg_catch, fill = scenario, color = scenario)) +
    geom_violin(
    trim = FALSE,           
    position = position_dodge(width = 0.8), 
    alpha = 0.7) +
    facet_grid(recruitment~skip, scales = "free_y") +
  scale_y_continuous(labels = scales::comma) +
  scico::scale_fill_scico_d(palette = "berlin") +
  scico::scale_color_scico_d(palette='roma') +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1))
  
  
all_in %>% 
  ggplot(aes(year, ssb_error_ii, group = interaction(iteration, skip), color = skip)) + 
  geom_line(alpha = 0.1) +
  facet_grid(recruitment ~ skip_shape)
