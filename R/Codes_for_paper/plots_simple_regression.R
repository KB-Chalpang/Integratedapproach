rm(list = ls())
# Load necessary libraries
library(dplyr)
library(plotly)
library(pROC)
library(readr)
# parameter values (i.e., intercept (β₀) and coefficient (β₁) values) that yield target 
#prevalence and target auc values. These are the parameter values we want our 
#algorithm to find.
parameter_values <- read_csv("Data/true_par_values.csv")
n <- 1000  # Number of combinations
set.seed(34532)
x1 <- rnorm(n)  # explanatory variable x1
intercepts <- runif(n, -5.5, 5.5)  # Random intercepts between -5 and 5
coef_x1 <- runif(n, -5.5, 5.5)  
combinations <- data.frame(intercept = intercepts, coef_x1 = coef_x1)
# Function to simulate outcomes and calculate AUC
calculate_auc <- function(intercept, coef_x1) {
  # Linear predictor
  linear_pred <- intercept + coef_x1 * x1
  # Convert linear predictor to probabilities using the logistic function
  probabilities <- 1 / (1 + exp(-linear_pred))
  # Generate binary outcomes based on probabilities
  y <- rbinom(n, 1, probabilities)
  # Check if 'y' has two levels; if not, return NA for AUC
  if (length(unique(y)) < 2) {
    return(NA)
  }
  # Calculate AUC, suppressing warnings
  auc <- suppressWarnings(roc(y, probabilities)$auc)
  return(auc)
}
# Step 4: Use apply to compute AUC for each combination
auc_values <- apply(combinations, 1, function(row) {
  calculate_auc(row['intercept'], row['coef_x1'])
})
# Add the AUC values to the combinations data frame
combinations$auc <- auc_values
# Filter the data frame to keep rows with AUC in the specified ranges
combinations_filtered <- combinations %>%
  filter((auc >= 0 & auc <= 0.65) | (auc > 0.68 & auc <= 0.79) | (auc > 0.8 & auc <= 1))
# Rearrange the rows from low AUC to high AUC
combinations_rearranged <- combinations_filtered %>%
  arrange(auc)
# Interpolate auc values for contour plot using coef_x1 and intercept
interp_results <- with(combinations_rearranged, 
                       akima::interp(x = intercept, y = coef_x1, z = auc, duplicate = "mean"))
# Determine midpoint for x and y based on contour data
x_mid <- mean(range(interp_results$x))
y_mid <- mean(range(interp_results$y))
# Set x and y ranges to a fixed 9-unit span centered around their respective midpoints
x_contour_range <- c(x_mid - 4.5, x_mid + 4.5)
y_contour_range <- c(y_mid - 4.5, y_mid + 4.5)
custom_colorscale <- list(
  list(0, "blue"),         # Starting color at the minimum (0.5)
  list(1, "lightblue")      # Ending color at the maximum
)
#cat("\u03B2\u2080")
#cat("\u03B2\u2081")
# Create the 2D contour plot with the custom color scale
fig_contour1 <- plot_ly(x = interp_results$x, 
                        y = interp_results$y, 
                        z = interp_results$z, 
                        type = 'contour', 
                        colorscale = custom_colorscale,   # Use the custom color scale
                        zmin = 0.45,                       # Set minimum AUC value
                        zmax = max(interp_results$z, na.rm = TRUE)) %>%
  layout(title = 'Intercept, Coefficient of x1, and AUC',
         xaxis = list(title = 'Intercept β₀' , range = x_contour_range),         # Label for x-axis
         yaxis = list(title = 'Coefficient β₁' , range = y_contour_range)) %>%
  colorbar(title = 'AUC')  # Title for the color legend as "AUC"
fig_contour1
###############################################################################
# Filter high AUC points within the new constrained 9-unit range and sample up to 20 points
high_auc_points <- combinations_rearranged %>%
  filter(auc >= 0.45,
         coef_x1 >= x_contour_range[1], coef_x1 <= x_contour_range[2],
         intercept >= y_contour_range[1], intercept <= y_contour_range[2]) %>%
  sample_n(min(400, n()))  # Sample up to 20 points or fewer if less are available
# Create the contour plot with updated x and y range
fig_contour2 <- plot_ly(x = interp_results$x, 
                        y = interp_results$y, 
                        z = interp_results$z, 
                        type = 'contour', 
                        colors = colorRamp(c("blue", "lightblue"))) %>%
  layout(title = 'Intercept, Coefficient of x, and AUC',
         xaxis = list(title = 'Intercept β₀', range = x_contour_range),  
         yaxis = list(title = 'Coefficient β₁', range = y_contour_range)) %>%
  colorbar(title = 'AUC')
# Overlay the filtered high AUC points on the contour plot
fig_contour2 <- fig_contour2 %>%
  add_trace(data = high_auc_points, 
            x = ~intercept,  
            y = ~coef_x1, 
            type = 'scatter',  
            mode = 'markers',   
            marker = list(color = 'red', size = 8),
            name = 'High AUC (>= 0.75)')
# Display the final plot
fig_contour2
################################################################################
# Create design matrix x
x <- cbind(1, x1)  # First column is 1's, second is the covariate
# Initialize results in data frame with prevalence
results <- high_auc_points %>%
  mutate(prevalence = numeric(nrow(high_auc_points)))  # Initialize the prevalence column
# Calculate prevalence for each row in high_auc_points using matrix multiplication
for (i in 1:nrow(high_auc_points)) {
  # Ensure the coefficients are treated as a numeric vector
  coeffs <- as.numeric(high_auc_points[i, 1:2])
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
# Interpolate auc values for contour plot using coef_x1 and intercept
interp_results_prev <- with(results, 
                            akima::interp(x = intercept, y = coef_x1, z = prevalence, duplicate = "mean"))
# Determine midpoint for x and y based on contour data
x_mid_prev <- mean(range(interp_results_prev$x))
y_mid_prev <- mean(range(interp_results_prev$y))
# Set x and y ranges to a fixed 9-unit span centered around their respective midpoints
x_contour_range_prev <- c(x_mid_prev - 3.8, x_mid_prev + 3.8)
y_contour_range_prev <- c(x_mid_prev - 4, y_mid_prev + 4)
# Create the contour plot with updated x and y range
fig_contour_prev <- plot_ly(x = interp_results_prev$x, 
                            y = interp_results_prev$y, 
                            z = interp_results_prev$z, 
                            type = 'contour', 
                            colors = colorRamp(c("yellow", "red"))) %>%
  layout(title = 'Intercept, Coefficient of x, and prevalence',
         xaxis = list(title = 'Intercet β₀', range = x_contour_range_prev),  
         yaxis = list(title = 'Coefficient β₁', range = y_contour_range_prev)) %>%
  colorbar(title = 'Prevalence')
fig_contour_prev
##### Intercept, coefficient, Prevalence, and AUC values  
# Determine the limits based on the scatter plot
x_limits <- range(high_auc_points$coef_x1)
y_limits <- range(high_auc_points$intercept)
# Create a 2D contour plot with axis labels and colorbar title
fig_contour3 <- plot_ly(x = interp_results$x, 
                        y = interp_results$y, 
                        z = interp_results$z, 
                        type = 'contour', 
                        colors = colorRamp(c("blue", "lightblue")), 
                        showscale = TRUE) %>%  # Show color scale for AUC
  layout(title = 'Intercept, Coefficient of x, AUC and prevalence',
         xaxis = list(title = 'Intercept β₀', range = x_limits),  # Set x-axis range
         yaxis = list(title = 'Coefficent β₁', range = y_limits)) %>%      # Set y-axis range
  colorbar(title = 'AUC',  # Title for the color legend as "AUC"
           thickness = 15, 
           len = 0.5,  # Length of the color bar
           x = 1.05,  # x position of the color bar
           y = 0.5)   # y position of the color bar
# Overlay the prevalence points
fig_contour3 <- fig_contour3 %>%
  add_trace(x = high_auc_points$coef_x1, 
            y = high_auc_points$intercept, 
            mode = 'markers', 
            type = 'scatter', 
            marker = list(size = 8, 
                          color = results$prevalence, 
                          colorscale = list(c(0, "yellow"), c(1, "red")), 
                          colorbar = list(title = 'Prevalence', 
                                          thickness = 15,  # Match thickness with AUC
                                          len = 0.5,      # Match length with AUC
                                          x = 1.05,       # Align with AUC
                                          y = 0.75)),    # Raise the position of prevalence
            name = 'Prevalence Points', 
            showlegend = TRUE)  # Show legend for prevalence points
# Display the plot
fig_contour3
##################################################################################
# Integrated approach  
fig_contour6 <- plot_ly(x = interp_results$x, 
                        y = interp_results$y, 
                        z = interp_results$z, 
                        type = 'contour', 
                        colors = colorRamp(c("blue", "lightblue")), 
                        showscale = TRUE) %>%  # Show color scale for AUC
  layout(title = 'Intercept, Coefficient of x, AUC, and prevalence',
         xaxis = list(title = ' Intercept β₀', range = x_limits),  # Set x-axis range
         yaxis = list(title = 'Coefficient β₁', range = y_limits)) %>%      # Set y-axis range
  colorbar(title = 'AUC',  # Title for the color legend as "AUC"
           thickness = 15, 
           len = 0.5,  # Length of the color bar
           x = 1.05,  # x position of the color bar
           y = 0.5) 
# Filter points where Prevalence is exactly 0.1, 0.5, or 0.8
#highlight_points1 <- parameter_values
# This just for symmetry of the target prevalence and AUC values, by negating 
# the coefficent of x1.
parameter_valuesCopy <- parameter_values
parameter_valuesCopy[, 3] <- -parameter_valuesCopy[,3]
combined_par_values <- rbind(parameter_values,parameter_valuesCopy)
# Add the target AUC values as diamond shapes
fig_contour6 <- fig_contour6 %>%
  add_trace(x = combined_par_values$beta1,
            y = combined_par_values$beta0, 
            mode = 'markers', 
            type = 'scatter', 
            marker = list(size = 15, 
                          color = "white", # Distinct color for visibility
                          symbol = "diamond", # Diamond shape for unique marking
                          line = list(color = "white", width = 1.5)), 
            name = 'Exact Prevalence Points',
            showlegend = F) %>%
  layout(legend = list(x = 1.1, y = 0.5)) # Adjust legend position
# Overlay the prevalence points
fig_contour6 <- fig_contour6 %>%
  add_trace(x = combined_par_values$beta1, 
            y = combined_par_values$beta0, 
            mode = 'markers', 
            type = 'scatter', 
            marker = list(size = 8, 
                          color = combined_par_values$Prevalence, 
                          colorscale = list(c(0, "yellow"), c(1, "red")), 
                          colorbar = list(title = 'Prevalence', 
                                          thickness = 15,  # Match thickness with AUC
                                          len = 0.5,      # Match length with AUC
                                          x = 1.05,       # Align with AUC
                                          y = 0.75)),    # Raise the position of prevalence
            name = 'Prevalence Points', 
            showlegend = F)  # Show legend for prevalence points
# Display the plot
fig_contour6
# Austin's 2-steps approach
two_steps1 <- read_csv("Data/results_from_two_steps_approach.csv")
#View(results_df)
two_steps2 <- two_steps1
two_steps2[, 3] <- -two_steps2[,3]
#highlight_points <- rbind(highlight_points1,highlight_points2)
combined_par_two_steps <- rbind(two_steps1, two_steps2)
fig_contour7 <- plot_ly(x = interp_results$x, 
                        y = interp_results$y, 
                        z = interp_results$z, 
                        type = 'contour', 
                        colors = colorRamp(c("blue", "lightblue")), 
                        showscale = TRUE) %>%  # Show color scale for AUC
  layout(title = 'Intercept, Coefficient of x, AUC, and prevalence',
         xaxis = list(title = 'Intercept β₀', range = x_limits),  # Set x-axis range
         yaxis = list(title = 'Coefficient β₁', range = y_limits)) %>%      # Set y-axis range
  colorbar(title = 'AUC',  # Title for the color legend as "AUC"
           thickness = 15, 
           len = 0.5,  # Length of the color bar
           x = 1.05,  # x position of the color bar
           y = 0.5) 
# Add filtered points to the contour plot as distinct markers
fig_contour7 <- fig_contour7 %>%
  add_trace(x = combined_par_values$beta1, 
            y = combined_par_values$beta0, 
            mode = 'markers', 
            type = 'scatter', 
            marker = list(size = 15, 
                          color = "white", # Distinct color for visibility
                          symbol = "diamond", # Diamond shape for unique marking
                          line = list(color = "white", width = 1.5)), 
            name = 'Exact Prevalence Points',
            showlegend = F) %>%
  layout(legend = list(x = 1.1, y = 0.5)) # Adjust legend position
# Overlay the prevalence points
fig_contour7 <- fig_contour7 %>%
  add_trace(x = combined_par_two_steps$beta1, 
            y = combined_par_two_steps$beta0, 
            mode = 'markers', 
            type = 'scatter', 
            marker = list(size = 8, 
                          color = combined_par_two_steps$Prevalence, 
                          colorscale = list(c(0, "yellow"), c(1, "red")), 
                          colorbar = list(title = 'Prevalence', 
                                          thickness = 15,  # Match thickness with AUC
                                          len = 0.5,      # Match length with AUC
                                          x = 1.05,       # Align with AUC
                                          y = 0.75)),    # Raise the position of prevalence
            name = 'Prevalence Points', 
            showlegend = F)  # Show legend for prevalence points
# Display the plot
fig_contour7
################################################################################
