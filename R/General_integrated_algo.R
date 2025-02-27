rm(list = ls())

# This code is a modified version of Austin PC. (2023) : 
#https://github.com/peter-austin/BMC_MRM-bisection-procedures-for-Monte-Carlo-simulations

# Install packages
#install.packages("rms")
#install.packages("pROC")
library(rms)
library(pROC)

N <- 10000
target.prev.outcome <- 0.3 # target prevalence
target.cstat <- 0.99
tolerance.prev <- 0.0001
tolerance.cstat <- 0.0001
max.iter <- 10000

# Generate covariates
set.seed(676869)
x1 <- rnorm(N, 0, 1)
x2 <- rnorm(N, 0, 1)
x3 <- rnorm(N, 0, 1)
x4 <- rnorm(N, 0, 1)
x5 <- rnorm(N, 0, 1)
x6 <- rbinom(N, 1, 0.5)
x7 <- rbinom(N, 1, 0.5)
x8 <- rbinom(N, 1, 0.5)
x9 <- rbinom(N, 1, 0.5)
x10 <- rbinom(N, 1, 0.5)
X <- cbind(1, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)

# Initial coefficients
B.outcome.fixed <- c(log(1.2), log(1.5), log(1.0), log(1.8), log(1.3),
                     log(1.4), log(1.7), log(1.9), log(1.6), log(1.55))

# Bisection intervals for prevalence and AUC
# Bisection intervals
# Depending on the target prevalence and AUC value, the intervals have to be set 
# accordingly to ensure quick convergence.

int.low <- -50
int.high <- 50
scalar.low <- 0
scalar.high <- 15

# Initialization
outcome.prev <- 0
cstat.emp <- 0
iter <- 0

# Joint optimization loop for both intercept and scalar
while ((abs(outcome.prev - target.prev.outcome) > tolerance.prev || abs(cstat.emp - target.cstat) > tolerance.cstat) && iter < max.iter) {
  iter <- iter + 1
  # Step 1: Adjust intercept to match prevalence
  int.mid <- (int.low + int.high) / 2
  scalar <- (scalar.low + scalar.high) / 2
  B.outcome <- scalar * B.outcome.fixed
  beta.outcome.modified <- c(int.mid, B.outcome)
  # Compute linear predictors and probabilities
  XB <- X %*% beta.outcome.modified
  p.outcome <- exp(XB) / (1 + exp(XB))
  Y <- rbinom(N, 1, p.outcome)
  
  # Step 2: Evaluate outcome prevalence
  outcome.prev <- mean(Y)
  # Step 3: Evaluate c-statistic (AUC)
  psm.pop <- lrm(Y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10)
  cstat.emp <- psm.pop$stats["C"]
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

#This code is better described as an iterative bisection adjustment method
#The key idea is to update the intercept to match the desired prevalence and 
#then update the coefficients to match the desired AUC

# check if the target AUC and target prevalence are met 
predicted_probabilities <- plogis(X%*%beta.outcome.modified)
# target AUC value
auc(Y, predicted_probabilities)
# target prevalence
mean(Y)
# The intercept and coefficients that yield the target prevalence and target AUC value
beta.outcome.modified
# The initial values of the coefficients
B.outcome.fixed
# It can be noticed that the coefficients that achieve our target prevalence
# and the target AUC values are significantly different. 
