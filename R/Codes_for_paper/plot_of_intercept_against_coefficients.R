# PLOTS SHOWING RELATIONSHIP BETWEEN AUC AND PREVALENCE

# This script will only run after running the R script plots_multiple_regression

# Interpolate auc values for contour plot using coef_x1 and intercept
interp_result1 <- with(combinations_rearranged, 
                       akima::interp(x = intercept, y = coef_x1, z = auc, duplicate = "mean"))

# Interpolate auc values for contour plot using coef_x2 and intercept
interp_result2 <- with(combinations_rearranged, 
                       akima::interp(x = intercept, y = coef_x2, z = auc, duplicate = "mean"))

# Interpolate auc values for contour plot using coef_x3 and intercept
interp_result3 <- with(combinations_rearranged, 
                       akima::interp(x = intercept, y = coef_x3, z = auc, duplicate = "mean"))

# Determine midpoint for x and y based on contour data
x_mid <- mean(range(interp_result1$x))
y_mid <- mean(range(interp_result1$y))

# Set x and y ranges to be centered around their respective midpoints
x_contour_range <- c(x_mid - 4.5, x_mid + 4.5)
y_contour_range <- c(y_mid - 4.5, y_mid + 4.5)


custom_colorscale <- list(
  list(0, "blue"),         # Starting color at the minimum (0.5)
  list(1, "lightblue")      # Ending color at the maximum
)

#cat("\u03B2\u2080")
#cat("\u03B2\u2081")
#cat("\u03B2\u2082")
#cat("\u03B2\u2083")

# Create the 2D contour plot with the custom color scale
fig_contour1 <- plot_ly(x = interp_results$x, 
                        y = interp_results$y, 
                        z = interp_results$z, 
                        type = 'contour', 
                        colorscale = custom_colorscale,   # Use the custom color scale
                        zmin = 0.45,                       # Set minimum AUC value
                        zmax = max(interp_results$z, na.rm = TRUE)) %>%
  layout(title = 'Intercept, Coefficient of x3, and AUC',
         xaxis = list(title = 'Intercept β₀' , range = x_contour_range),         # Label for x-axis
         yaxis = list(title = 'Coefficient β₁' , range = y_contour_range)) %>%
  colorbar(title = 'AUC')  # Title for the color legend as "AUC"

fig_contour1


# Create the 2D contour plot with the custom color scale
fig_contour2 <- plot_ly(x = interp_result2$x, 
                        y = interp_result2$y, 
                        z = interp_result2$z, 
                        type = 'contour', 
                        colorscale = custom_colorscale,   # Use the custom color scale
                        zmin = 0.45,                       # Set minimum AUC value
                        zmax = max(interp_result2$z, na.rm = TRUE)) %>%
  layout(title = 'Intercept, Coefficient of x2, and AUC',
         xaxis = list(title = 'Intercept β₀' , range = x_contour_range),         # Label for x-axis
         yaxis = list(title = 'Coefficient β₂' , range = y_contour_range)) %>%
  colorbar(title = 'AUC')  # Title for the color legend as "AUC"

fig_contour2


# Create the 2D contour plot with the custom color scale
fig_contour3 <- plot_ly(x = interp_result3$x, 
                        y = interp_result3$y, 
                        z = interp_result3$z, 
                        type = 'contour', 
                        colorscale = custom_colorscale,   # Use the custom color scale
                        zmin = 0.45,                       # Set minimum AUC value
                        zmax = max(interp_result3$z, na.rm = TRUE)) %>%
  layout(title = 'Intercept, Coefficient of x3, and AUC',
         xaxis = list(title = 'Intercept β₀' , range = x_contour_range),         # Label for x-axis
         yaxis = list(title = 'Coefficient β₃' , range = y_contour_range)) %>%
  colorbar(title = 'AUC')  # Title for the color legend as "AUC"

fig_contour3


# Determine midpoint for x and y based on contour data
#x_mid <- mean(range(interp_result1$x))
#y_mid <- mean(range(interp_result1$y))

# Set x and y ranges to a fixed 9-unit span centered around their respective midpoints
x_contour_range <- c(x_mid - 4.5, x_mid + 4.5)
y_contour_range <- c(y_mid - 4.5, y_mid + 4.5)

# Filter high AUC points within the new constrained 9-unit range and sample up to 20 points
high_auc_points <- combinations_rearranged %>%
  filter(auc >= 0.45,
         coef_x1 >= x_contour_range[1], coef_x3 <= x_contour_range[2],
         intercept >= y_contour_range[1], intercept <= y_contour_range[2]) %>%
  sample_n(min(700, n()))  # Sample up to 20 points or fewer if less are available

################################################################################

# The design matrix x
x <- cbind(1, x1, x2, x3) 

# Initialize results data frame with prevalence
results <- high_auc_points %>%
  mutate(prevalence = numeric(nrow(high_auc_points)))  # Initialize the prevalence column

# Calculate prevalence for each row in high_auc_points using matrix multiplication
for (i in 1:nrow(high_auc_points)) {
  # Ensure the coefficients are treated as a numeric vector
  coeffs <- as.numeric(high_auc_points[i, 1:ncol(x)])
  
  # Calculate the probability
  prob <- 1 / (1 + exp(- (x %*% coeffs)))
  
  # Simulate binary outcomes based on the probability
  y <- rbinom(n, 1, prob)  # Generate binary outcomes
  
  # Calculate the prevalence as the mean of the simulated outcomes
  results$prevalence[i] <- mean(y)  # Store the prevalence
}

# Display the results with prevalence
head(results)

################################################################################

interp_results_prev1 <- with(results, 
                            akima::interp(x = intercept, y = coef_x1, z = prevalence, duplicate = "mean"))

interp_results_prev2 <- with(results, 
                            akima::interp(x = intercept, y = coef_x2, z = prevalence, duplicate = "mean"))

interp_results_prev3 <- with(results, 
                            akima::interp(x = intercept, y = coef_x3, z = prevalence, duplicate = "mean"))


# Determine midpoint for x and y based on contour data
x_mid_prev <- mean(range(interp_results_prev1$x))
y_mid_prev <- mean(range(interp_results_prev1$y))

# Set x and y ranges to a fixed 9-unit span centered around their respective midpoints
x_contour_range_prev <- c(x_mid_prev - 4, x_mid_prev + 4)
y_contour_range_prev <- c(y_mid_prev - 4, y_mid_prev + 4)


# Create a 2D contour plot with axis labels and colorbar title
fig_contour_prev1 <- plot_ly(x = interp_results_prev1$y, 
                            y = interp_results_prev1$x, 
                            z = interp_results_prev1$z, 
                            type = 'contour', 
                            colors = colorRamp(c("yellow", "red"))) %>%
  layout(title = 'Intercept, Coefficient of x1, and Prevalence',
         xaxis = list(title = 'Intercept β₀', range = x_contour_range_prev),         # Label for x-axis
         yaxis = list(title = 'Coefficient β₁', range = y_contour_range_prev)) %>%
  colorbar(title = 'Prevalence')  # Title for the color legend as "AUC"

fig_contour_prev1


# Create a 2D contour plot with axis labels and colorbar title
fig_contour_prev2 <- plot_ly(x = interp_results_prev2$y, 
                            y = interp_results_prev2$x, 
                            z = interp_results_prev2$z, 
                            type = 'contour', 
                            colors = colorRamp(c("yellow", "red"))) %>%
  layout(title = 'Intercept, Coefficient of x2, and Prevalence',
         xaxis = list(title = 'Intercept β₀', range = x_contour_range_prev),         # Label for x-axis
         yaxis = list(title = 'Coefficient β₂', range = y_contour_range_prev)) %>%
  colorbar(title = 'Prevalence')  # Title for the color legend as "AUC"

fig_contour_prev2

# Create a 2D contour plot with axis labels and colorbar title
fig_contour_prev3 <- plot_ly(x = interp_results_prev3$y, 
                            y = interp_results_prev3$x, 
                            z = interp_results_prev3$z, 
                            type = 'contour', 
                            colors = colorRamp(c("yellow", "red"))) %>%
  layout(title = 'Intercept, Coefficient of x3, and Prevalence',
         xaxis = list(title = 'Intercept β₀', range = x_contour_range_prev),         # Label for x-axis
         yaxis = list(title = 'Coefficient β₃', range = y_contour_range_prev)) %>%
  colorbar(title = 'Prevalence')  # Title for the color legend as "AUC"

fig_contour_prev3
