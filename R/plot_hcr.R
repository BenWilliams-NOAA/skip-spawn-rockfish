library(ggplot2)
library(patchwork)
# 1. Define your function
tier_3 <- function(b, b40, f40) {
  b <- as.numeric(b)
  b40 <- as.numeric(b40)
  f40 <- as.numeric(f40)
  
  if(length(b) == 0 || is.na(b)) return(0)
  
  bb40 <- b / b40
  
  if(bb40 > 1) {
    return(f40)
  } else if(bb40 > 0.05 && bb40 <= 1) {
    return(f40 * (bb40 - 0.05) / (1 - 0.05))
  } else {
    return(0)
  }
}

# 2. Vectorize the function so it works with ggplot
# This allows the function to accept a list of 'b' values rather than just one
tier_3_vec <- Vectorize(tier_3, vectorize.args = "b")

# 3. Set parameters for the plot
B40_ref <- 1   # Arbitrary reference biomass
F40_ref <- 0.2   # Arbitrary reference fishing mortality
max_b   <- 1.50   # Max biomass to plot (go slightly above B40 to show the flat top)

# 4. Generate data
# Create a sequence of biomass values from 0 to max_b
b_values <- seq(0, max_b, length.out = 200)

# Calculate F for each biomass value
f_values <- tier_3_vec(b = b_values, b40 = B40_ref, f40 = F40_ref)

# Create a data frame for plotting
hcr_data <- data.frame(
  Biomass = b_values,
  Fishing_Mortality = f_values,
  Status = b_values / B40_ref, # Calculate B/B40 ratio for the secondary axis
  Catch = b_values * f_values
)

# 5. Plot using ggplot2
ggplot(hcr_data, aes(x = Biomass, y = Fishing_Mortality)) +
  geom_line( size = 1.2) +
  
  # Add visual markers for the "Knee" and cut-off
  geom_vline(xintercept = B40_ref, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = B40_ref * 0.05, linetype = "dotted", color = "red") +
  
  # Labels and Theme
  labs(
    title = "Tier 3 HCR",
    # subtitle = paste0("Parameters: B40 = ", B40_ref, ", F40 = ", F40_ref),
    x = "Biomass (B40)",
    y = "Fishing Mortality (F)"
  ) +
  
  # Annotations to explain the shape
  annotate("text", x = B40_ref , y = F40_ref, label = "Max F (F40)", vjust = -0.5, hjust = -0.025) +
  annotate("text", x = B40_ref/2, y = F40_ref/4, label = "Linear Ramp") +
  expand_limits(y = 0.25) +
  
  afscassess::theme_report(base_size = 16) -> p1



# 5. Plot Catch
ggplot(hcr_data, aes(x = Biomass, y = Catch)) +
  geom_line(size = 1.2) +
  
  # Add markers
  geom_vline(xintercept = B40_ref, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = B40_ref * 0.05, linetype = "dotted", color = "red") +
  
  # Labels
  labs(
    title = "Projected Catch (Yield)",
    # subtitle = paste0("Derived from Tier 3 HCR (Catch = Biomass * F)"),
    x = "Biomass (B/B40)",
    y = "Catch (Yield)"
  ) +
  
  # Annotations explaining the geometry
  annotate("text", x = B40_ref, y = (B40_ref)*F40_ref, 
           label = "Constant F", vjust = 1, hjust = -0.025) +
  annotate("text", x = B40_ref/1.5, y = 0.1, 
           label = "Ramping F", hjust = 1) +
  
  afscassess::theme_report(base_size = 16) -> p2

p1 + p2


# Calculate the Ratios
hcr_data$BB40 <- hcr_data$Biomass / B40_ref
hcr_data$FF40 <- hcr_data$Fishing_Mortality / F40_ref

# 5. Plot the Normalized Curve
ggplot(hcr_data, aes(x = BB40, y = FF40)) +
  geom_line( size = 1.2) +
  
  # Reference lines at 1 (The Target)
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  # geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  
  # The Cutoff line at 0.05 (alpha)
  geom_vline(xintercept = 0.05, linetype = "dotted", color = "red") +
  
  # Labels with LaTeX formatting
  labs(
    title = "Tier 3 HCR",
    x = expression(Biomass~(B/B[40])),
    y = expression(Fishing~Mortality~(F/F[40]))
  ) +
  afscassess::theme_report(base_size = 16) -> p1


ggplot(hcr_data, aes(x = Biomass, y = Catch)) +
  geom_line(size = 1.2) +
  
  # Add markers
  geom_vline(xintercept = B40_ref, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = B40_ref * 0.05, linetype = "dotted", color = "red") +
  
  # Labels
  labs(
    title = "Projected Catch (Yield)",
    # subtitle = paste0("Derived from Tier 3 HCR (Catch = Biomass * F)"),
    x = expression(Biomass~(B/B[40])),
    y = "Catch"
  ) +
  
  # Annotations explaining the geometry
  annotate("text", x = B40_ref, y = (B40_ref)*F40_ref, 
           label = "Constant F", vjust = 1, hjust = -0.025) +
  annotate("text", x = B40_ref/1.5, y = 0.1, 
           label = "Ramping F", hjust = 1) +
  
  afscassess::theme_report(base_size = 16) -> p2


p1 + p2
