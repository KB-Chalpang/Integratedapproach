rm(list = ls())

# Before running this code, you should first run the code "Integrated_algorithm"
# It should be noted that how fast convergence is achieved depends 
# on the tolerance level

# Initialization
outcome.prev <- 0
cstat.emp <- 0
iter <- 0

# Define the number of parameters to store (intercept, scalar, cstat, and slopes)
num_params <- length(B.outcome.fixed) + 3  # +3 for intercept, scalar, and c-stat
results_matrix <- matrix(NA, nrow = max.iter, ncol = num_params)

# Column names for the results matrix
colnames(results_matrix) <- c("Intercept", "Scalar", "Cstat", paste0("Slope_", 1:length(B.outcome.fixed)))

# Initialize vectors to store the differences for convergence plot
outcome_diff <- numeric(max.iter)   # Difference between outcome prevalence and target
cstat_diff <- numeric(max.iter)     # Difference between empirical and target c-statistics

# Joint optimization loop for both intercept and scalar
while ((abs(outcome.prev - target.prev.outcome) > tolerance.prev || abs(cstat.emp - target.cstat) > tolerance.cstat) && iter < max.iter) {
  iter <- iter + 1
  
  # Step 1: Adjust intercept to match prevalence
  int.mid <- (int.low + int.high) / 2
  scalar <- (scalar.low + scalar.high) / 2
  
  B.outcome <- scalar * B.outcome.fixed
  beta.outcome.modified <- c(int.mid, B.outcome)
  
  # Compute linear predictor and probabilities
  XB <- X %*% beta.outcome.modified
  p.outcome <- exp(XB) / (1 + exp(XB))
  Y <- rbinom(N, 1, p.outcome)
  
  # Step 2: Evaluate outcome prevalence
  outcome.prev <- mean(Y)
  
  # Step 3: Evaluate c-statistic (AUC)
  model <- lrm(Y ~ x1+x2+x3)
  cstat.emp <- model$stats["C"]
  
  # Update intervals for both intercept and scalar based on prevalence and AUC
  if (outcome.prev < target.prev.outcome) {
    int.low <- int.mid
  } else {
    int.high <- int.mid
  }
  
  if (cstat.emp < target.cstat) {
    scalar.low <- scalar
  } else {
    scalar.high <- scalar
  }
  
  # Store the intercept, scalar, c-stat, and slopes in the matrix
  results_matrix[iter, ] <- c(int.mid, scalar, cstat.emp, B.outcome)
  
  # Display current status
  cat("Iteration:", iter, "\n",
      "Prevalence: ", outcome.prev, "Target Prevalence: ", target.prev.outcome, "\n",
      "C-Statistic: ", cstat.emp, "Target C-Statistic: ", target.cstat, "\n",
      "Current intercept: ", int.mid, "\n",
      "Current scalar: ", scalar, "\n\n")
  
  # Store the differences for convergence tracking
  outcome_diff[iter] <- outcome.prev - target.prev.outcome
  cstat_diff[iter] <- cstat.emp - target.cstat
}

# Truncate the difference vectors to the actual number of iterations
outcome_diff <- outcome_diff[1:iter]
cstat_diff <- cstat_diff[1:iter]

# Plot the convergence of differences
par(mfrow = c(1, 2)) # Set up a 2-row plot

# Plot for outcome prevalence difference
plot(outcome_diff, type = "l", col = "blue", main = "Prevalence = 0.8",
     xlab = "Iteration", ylab = "Prevalence residuals")
abline(h = 0, col = "red")  # Add a horizontal line at y = 0


# Plot for c-statistic difference
plot(cstat_diff, type = "l", col = "blue", main = "AUC = 0.85",
     xlab = "Iteration", ylab = "AUC residuals")
abline(h = 0, col = "red")  # Add a horizontal line at y = 0


