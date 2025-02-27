rm(list = ls())

library(rms)
library(pROC)

# Parameters
N <- 2000
target.prev.outcome <- 0.1
target.cstat <- 0.85
tolerance.prev <- 0.001
tolerance.cstat <- 0.001
max.iter <- 10000

# Generate covariates
set.seed(4555632)
x1 <- rnorm(N, 0, 1)
x2 <- rnorm(N, 0, 1)
x3 <- rnorm(N, 0, 1)
X <- cbind(1, x1, x2, x3)

# Initial coefficients
B.outcome.fixed <- c(log(2), log(2), log(2))

# Bisection intervals
int.low <- -20
int.high <- 20
scalar.low <- 0
scalar.high <- 20

# Initialization
outcome.prev <- 0
cstat.emp <- 0
iter <- 0

# Define the number of parameters to store (intercept, scalar, cstat, and slopes)
num_params <- length(B.outcome.fixed) + 3  # +3 for intercept, scalar, and c-stat
results_matrix <- matrix(NA, nrow = max.iter, ncol = num_params)

# Column names for the results matrix
colnames(results_matrix) <- c("Intercept", "Scalar", "Cstat", paste0("Coefficient_", 1:length(B.outcome.fixed)))

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
  model <- lrm(Y ~ x1 + x2 + x3)
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
  
  # Stop if both criteria are met within tolerance
  if (abs(outcome.prev - target.prev.outcome) <= tolerance.prev && abs(cstat.emp - target.cstat) <= tolerance.cstat) {
    cat("Converged after", iter, "iterations.\n")
    break
  }
}

cat("Final intercept: ", int.mid, "\nFinal scalar: ", scalar, "\n")

# Remove empty rows from the matrix (in case convergence happened before max.iter)
results_matrix <- results_matrix[1:iter, ]

# View the stored results
head(results_matrix)
tail(results_matrix)

# check if the target AUC and target prevalence are met 
predicted_probabilities <- plogis(X%*%beta.outcome.modified)
# target AUC value
auc(Y, predicted_probabilities)
# target prevalence
mean(Y)
# The intercept and coefficients that yield the target prevalence and target AUC value
beta.outcome.modified



