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
# relative error
# compare biological estimate (spawn_bio_r) to truth (spawn_bio_f)
# what did the assessment get wwrong - bias?
plot_bias <- function(sims, rmv_yrs=15) {
  sims %>%
    filter(year > (min(year) + rmv_yrs)) %>% # Skip the F-shock burn-in
    mutate(rel_error = (spawn_bio_r - spawn_bio_f) / spawn_bio_f) %>% 
    ggplot(aes(year, rel_error)) + 
    stat_summary(fun.data = median_mad,
                 geom = "ribbon",
                 alpha = 0.2,
                 fill = "cornflowerblue") +
    geom_line(aes(group = iteration), alpha = 0.2, color = "gray50") + # Individual worms
    stat_summary(fun.y = median, geom = "line", color = "cornflowerblue") + # Median trend
    geom_hline(yintercept = 0, linetype = 3) +
    facet_wrap(~scenario) +
    labs(title = "Assessment Bias",
         subtitle = "Does ignoring skip spawning cause over/under estimation?",
         y = "Relative Error (Est - True) / True", x = "Year")
}

# plot_mgmt <- function(sims, rmv_yrs = 15) {
#   # manager's view
#   sims %>%
#     mutate(scenario = ifelse(scenario == "i", 
#                              "Mgmt: Functional (Correct)", 
#                              "Mgmt: Biological (Misspecified)")) %>% 
#     filter(year > (min(year) + rmv_yrs)) %>% # Analyze long-term equilibrium
#     arrange(scenario, iteration, year) %>%
#     mutate(true_status = spawn_bio_r / B35, 
#            true_intensity = Ft / F35) %>% 
#     ggplot(aes(true_status, true_intensity)) + 
#     annotate("rect", xmin = -Inf, xmax = 1, ymin = 1, ymax = Inf, 
#              fill = "red", alpha = 0.1) +
#     annotate("rect",xmin = 1, xmax = Inf, ymin = -Inf, ymax = 1, 
#              fill = "green", alpha= 0.1) +
#     annotate("rect",xmin = 1, xmax = Inf, ymin = 1, ymax = Inf, 
#              fill = "yellow", alpha = 0.1) +
#     annotate("rect",xmin = -Inf, xmax = 1, ymin = -Inf, ymax = 1, 
#              fill = "yellow", alpha = 0.1) +
#     # geom_density_2d_filled(contour_var = "ndensity", alpha = 0.2) +
#     geom_density_2d(color="gray40") +# target crosshairs
#     geom_vline(xintercept = 1, linetype = "dashed", color = "white") +
#     geom_hline(yintercept = 1, linetype = "dashed", color = "white") +
#     geom_point(aes(color=iteration)) +
#     facet_wrap(~scenario) +
#     labs(title = "Manager's View",
#          subtitle = "Managing based on biological maturity",
#          x = "True SSB / True B35",
#          y = "Realized F / True F35") +
#     theme(legend.position = "none") +
#     expand_limits(x = 0.5, y = c(0.5, 1.5))
# }



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
    geom_vline(xintercept = max(hist_F$year), linetype =3) +
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
    geom_vline(xintercept = max(hist_F$year), linetype =3) +
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
    ggplot(aes(year, catch, group = interaction(scenario, iteration), color = scenario)) + 
    geom_line() + 
    geom_vline(xintercept = max(hist_F$year), linetype =3) 
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
                 color = 1,
                 linewidth = 1) +
    geom_vline(xintercept = max(hist_F$year), linetype =3) 
}

plot_risk <- function(sims, rmv_yrs = 15, ssb_col = "spawn_bio_f", bio_ref_col = "B35_f", F_col = "F35_f") {
  # manager's view
  sims %>%
    mutate(scenario = ifelse(scenario == "i", 
                             "Mgmt: Functional (Correct)", 
                             "Mgmt: Biological (Misspecified)")) %>% 
    filter(year > (min(year) + rmv_yrs)) %>% # Analyze long-term equilibrium
    arrange(scenario, iteration, year) %>%
    mutate(true_status = .data[[ssb_col]] / .data[[bio_ref_col]], 
           true_intensity = Ft / .data[[F_col]]) %>% 
    ggplot(aes(true_status, true_intensity)) + 
    annotate("rect", xmin = -Inf, xmax = 1, ymin = 1, ymax = Inf, 
             fill = "red", alpha = 0.1) +
    annotate("rect",xmin = 1, xmax = Inf, ymin = -Inf, ymax = 1, 
             fill = "green", alpha= 0.1) +
    annotate("rect",xmin = 1, xmax = Inf, ymin = 1, ymax = Inf, 
             fill = "yellow", alpha = 0.1) +
    annotate("rect",xmin = -Inf, xmax = 1, ymin = -Inf, ymax = 1, 
             fill = "orange", alpha = 0.1) +
    # geom_density_2d_filled(contour_var = "ndensity", alpha = 0.2) +
    geom_density_2d(color="gray40") +# target crosshairs
    geom_vline(xintercept = 1, linetype = "dashed", color = "white") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "white") +
    geom_point(aes(color=iteration)) +
    facet_wrap(~scenario) +
    labs(title = "Biological Risk Comparison",
         subtitle = "Consequences of ignoring skip spawning (Left) vs. Correctly accounting for it (Right)",
         x = "True SSB / True B35",
         y = "Realized F / True F35") +
    theme(legend.position = "none") + 
    expand_limits(x = c(0, 2), y = c(0, 2))
}

# risk table
# calculate true ratios (same as plot)
# scenario i reference points for the denominator
table_risk <- function(sims, rmv_yrs = 15, ssb_col = "spawn_bio_f", bio_ref_col = "B35_f", F_col = "F35_f") {
  sims %>%
    filter(year > (min(year) + rmv_yrs)) %>% 
    # left_join(
    #   sims %>% filter(scenario == "i") %>% 
    #     select(iteration, year, 
    #            TRUE_B35 = all_of(bio_ref_col), 
    #            TRUE_F35 = all_of(F_col)),
    #   by = c("iteration", "year")) %>%
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
