rm(list = ls())

# Load necessary libraries
library(dplyr)
library(plotly)
library(pROC)

set.seed(45632)
n <- 3000  # Number of combinations
intercepts <- runif(n, -5.5, 5.5)  # Random intercepts between -5 and 5
coef_x1 <- runif(n, -5, 5)  
coef_x2 <- runif(n, -5, 5)
coef_x3 <- runif(n, -5, 5)

combinations <- data.frame(intercept = intercepts, coef_x1 = coef_x1, coef_x2 = coef_x2,
                           coef_x3 = coef_x3)

# Step 3: Simulate a binary response variable based on logistic regression
# For each combination, generate x1 and x2 values
x1 <- rnorm(n)  
x2 <- rnorm(n)  
x3 <- rbinom(n,1,0.5)
#x4 <- rnorm(n)
# Function to simulate outcomes and calculate AUC
calculate_auc <- function(intercept, coef_x1, coef_x2, coef_x3) {
  # Linear predictor
  linear_pred <- intercept + coef_x1*x1 + coef_x2*x2 + coef_x3*x3
  
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
  calculate_auc(row['intercept'], row['coef_x1'], row['coef_x2'], row['coef_x3'])
})

# Add the AUC values to the combinations data frame
combinations$auc <- auc_values
head(combinations)


# Filter the data frame to keep rows with AUC in the specified ranges
combinations_filtered <- combinations %>%
  filter((auc >= 0 & auc <= 0.65) | (auc > 0.68 & auc <= 0.78) | (auc > 0.8 & auc <= 1))

# Rearrange the rows from low AUC to high AUC
combinations_rearranged <- combinations_filtered %>%
  arrange(auc)

#combinations_rearranged <- combinations_rearranged[10:nrow(combinations_filtered),]

head(combinations_rearranged)

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
cat("\u03B2\u2082")
#cat("\u03B2\u2083")

# Create the 2D contour plot with the custom color scale
fig_contour1 <- plot_ly(x = interp_results$x, 
                        y = interp_results$y, 
                        z = interp_results$z, 
                        type = 'contour', 
                        colorscale = custom_colorscale,   # Use the custom color scale
                        zmin = 0.45,                       # Set minimum AUC value
                        zmax = max(interp_results$z, na.rm = TRUE)) %>%
  layout(title = 'Coefficient x1, Coefficient of x3, and AUC',
         xaxis = list(title = 'Coefficent β₃' , range = x_contour_range),         # Label for x-axis
         yaxis = list(title = 'Coefficient β₁' , range = y_contour_range)) %>%
  colorbar(title = 'AUC')  # Title for the color legend as "AUC"

fig_contour1


# Determine midpoint for x and y based on contour data
x_mid <- mean(range(interp_results$x))
y_mid <- mean(range(interp_results$y))

# Set x and y ranges to a fixed 9-unit span centered around their respective midpoints
x_contour_range <- c(x_mid - 4.5, x_mid + 4.5)
y_contour_range <- c(y_mid - 4.5, y_mid + 4.5)

# Filter high AUC points within the new constrained 9-unit range and sample up to 20 points
high_auc_points <- combinations_rearranged %>%
  filter(auc >= 0.45,
         coef_x1 >= x_contour_range[1], coef_x1 <= x_contour_range[2],
         intercept >= y_contour_range[1], intercept <= y_contour_range[2]) %>%
  sample_n(min(700, n()))  # Sample up to 20 points or fewer if less are available

# Create the contour plot with updated x and y range
#fig_contour2 <- plot_ly(x = interp_results$x, 
#                        y = interp_results$y, 
#                        z = interp_results$z, 
#                        type = 'contour', 
#                        colors = colorRamp(c("blue", "lightblue"))) %>%
#  layout(title = 'Intercept, Coefficient of x, and AUC',
#         xaxis = list(title = 'Coefficient of x', range = x_contour_range),  
#         yaxis = list(title = 'Intercept', range = y_contour_range)) %>%
#  colorbar(title = 'AUC')
#
# Overlay the filtered high AUC points on the contour plot
#fig_contour2 <- fig_contour2 %>%
#  add_trace(data = high_auc_points, 
#            x = ~coef_x1,  
#            y = ~intercept, 
#            type = 'scatter',  
#            mode = 'markers',   
#            marker = list(color = 'red', size = 8),
#            name = 'High AUC (>= 0.75)')

# Display the final plot
#fig_contour2

################################################################################

# Create an n x 2 matrix x
x <- cbind(1, x1, x2, x3)  # First column is 1's, second is the covariate

# Initialize results data frame with prevalence
results <- high_auc_points %>%
  mutate(prevalence = numeric(nrow(high_auc_points)))  # Initialize the prevalence column

# Calculate prevalence for each row in high_auc_points using matrix multiplication
#set.seed(123)  # For reproducibility
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

interp_results_prev <- with(results, 
                            akima::interp(x = intercept, y = coef_x1, z = prevalence, duplicate = "mean"))


# Determine midpoint for x and y based on contour data
x_mid_prev <- mean(range(interp_results_prev$x))
y_mid_prev <- mean(range(interp_results_prev$y))

# Set x and y ranges to a fixed 9-unit span centered around their respective midpoints
x_contour_range_prev <- c(x_mid_prev - 4, x_mid_prev + 4)
y_contour_range_prev <- c(y_mid_prev - 4, y_mid_prev + 4)


# Create a 2D contour plot with axis labels and colorbar title
fig_contour_prev <- plot_ly(x = interp_results_prev$y, 
                            y = interp_results_prev$x, 
                            z = interp_results_prev$z, 
                            type = 'contour', 
                            colors = colorRamp(c("yellow", "red"))) %>%
  layout(title = 'Intercept, Coefficient of x, and Prevalence',
         xaxis = list(title = 'Intercept β₀', range = x_contour_range_prev),         # Label for x-axis
         yaxis = list(title = 'Coefficient β₁', range = y_contour_range_prev)) %>%
  colorbar(title = 'Prevalence')  # Title for the color legend as "AUC"

fig_contour_prev


# Determine the limits based on your scatter plot
x_limits <- range(high_auc_points$coef_x1)
y_limits <- range(high_auc_points$intercept)


# Create a 2D contour plot with axis labels and colorbar title
fig_contour3 <- plot_ly(x = interp_results$x, 
                        y = interp_results$y, 
                        z = interp_results$z, 
                        type = 'contour', 
                        colors = colorRamp(c("blue", "lightblue")), 
                        showscale = TRUE) %>%  # Show color scale for AUC
  layout(title = 'Intercept, Coefficient of x, and AUC',
         xaxis = list(title = 'Coefficient of x', range = x_limits),  # Set x-axis range
         yaxis = list(title = 'Intercept', range = y_limits)) %>%      # Set y-axis range
  colorbar(title = 'AUC',  # Title for the color legend as "AUC"
           thickness = 15, 
           len = 0.5,  # Length of the color bar
           x = 1.05,  # x position of the color bar
           y = 0.5)   # y position of the color bar

# Overlay the prevalence points
fig_contour3 <- fig_contour3 %>%
  add_trace(x = results$coef_x1, 
            y = results$intercept, 
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

