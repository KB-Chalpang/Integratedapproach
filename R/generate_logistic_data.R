# Generate logistic-regression data using nested bisection (Algorithm 3).
# Edit the two targets below, then source this script or run it in RStudio.
# Required package: pROC.
# Outputs are kept in R: calibration, coefficients, simulated_data and trace.
# Convergence plots are displayed in the current graphics device.

# 1. Choose the targets -------------------------------------------------------

target_prevalence <- 0.8
target_auc <- 0.85

# 2. Optional settings --------------------------------------------------------

# Predictors may be continuous, binary or a mixture. The defaults below use
# two continuous predictors and one binary predictor. All slopes share a
# multiplier; their relative sizes are held fixed during calibration.
N <- 100000L
base_slopes <- c(x1 = log(2), x2 = log(1.2), x3 = log(2.8))
prevalence_tolerance <- 0.001
auc_tolerance <- 0.001
max_inner_iterations <- 1000L
max_outer_iterations <- 1000L
auc_replicates <- 5L
intercept_bounds <- c(-20, 20)
multiplier_bounds <- c(0, 20)
covariate_seed <- 4555632L
calibration_seed <- 17L
data_seed <- 2026L
show_convergence_plots <- TRUE

# 3. Check inputs and generate fixed covariates --------------------------------

if (!requireNamespace("pROC", quietly = TRUE)) {
  stop("Install pROC before running this script: install.packages('pROC')")
}
check_probability <- function(value, name) {
  if (length(value) != 1L || !is.finite(value) || value <= 0 || value >= 1) {
    stop(name, " must be a number strictly between 0 and 1.")
  }
}
check_probability(target_prevalence, "target_prevalence")
check_probability(target_auc, "target_auc")
stopifnot(N >= 2, N == floor(N), length(base_slopes) == 3L,
          all(is.finite(base_slopes)), prevalence_tolerance > 0,
          auc_tolerance > 0, max_inner_iterations >= 1,
          max_outer_iterations >= 1, auc_replicates >= 2,
          length(intercept_bounds) == 2L, all(is.finite(intercept_bounds)),
          intercept_bounds[1] < intercept_bounds[2],
          length(multiplier_bounds) == 2L, all(is.finite(multiplier_bounds)),
          multiplier_bounds[1] >= 0, multiplier_bounds[1] < multiplier_bounds[2])

# Change these three lines to use your desired predictor distributions.
# rnorm() generates continuous values; rbinom(..., size = 1) generates binary
# values. All three predictors can be continuous, binary, or mixed.
# The default draws are independent. Keep the column names x1, x2 and x3;
# base_slopes supplies their relative effects in that order.
# Recalibrate after changing the covariate distributions or their values,
# because the coefficients needed to achieve the targets can change.
set.seed(covariate_seed)
covariates <- data.frame(
  x1 = rnorm(N, mean = 0, sd = 1),   # Continuous predictor.
  x2 = rbinom(N, size = 1, prob = 0.5), # Binary predictor: 0 or 1.
  x3 = rnorm(N, mean = 0, sd = 1)    # Continuous predictor.
)
X <- cbind(intercept = 1, as.matrix(covariates))

# 4. Calibration helpers ------------------------------------------------------

# AUC requires both outcome classes. Keep pROC's default direction selection.
# Convert predictions to a vector to avoid the deprecated matrix-input warning.
calculate_auc <- function(y, probability) {
  if (length(unique(y)) < 2L) {
    stop("An outcome draw contains only one class. Increase N or use less extreme targets.")
  }
  as.numeric(pROC::auc(y, as.numeric(probability), quiet = TRUE))
}

# For a proposed slope multiplier, find the intercept giving mean(p) equal
# to the target prevalence on these fixed covariates.
find_intercept <- function(multiplier) {
  slopes <- multiplier * base_slopes
  lo <- intercept_bounds[1]
  hi <- intercept_bounds[2]
  expected_prevalence <- function(intercept) {
    mean(plogis(X %*% c(intercept, slopes)))
  }
  if (expected_prevalence(lo) > target_prevalence ||
      expected_prevalence(hi) < target_prevalence) {
    stop("The intercept bounds do not bracket the target. Widen intercept_bounds.")
  }
  for (iteration in seq_len(max_inner_iterations)) {
    intercept <- (lo + hi) / 2
    prevalence <- expected_prevalence(intercept)
    if (prevalence < target_prevalence) lo <- intercept else hi <- intercept
    if (abs(prevalence - target_prevalence) <= prevalence_tolerance) break
  }
  list(intercept = intercept, slopes = slopes, prevalence = prevalence)
}

# 5. Run Algorithm 3 ----------------------------------------------------------

# Each outer iteration proposes a multiplier and recalibrates the intercept.
# AUC is averaged over several simulated outcome vectors. It is a noisy
# estimate: bisection is not guaranteed to reach every requested AUC.
lo <- multiplier_bounds[1]
hi <- multiplier_bounds[2]
history <- vector("list", max_outer_iterations)
converged <- FALSE

for (iteration in seq_len(max_outer_iterations)) {
  multiplier <- (lo + hi) / 2
  inner <- find_intercept(multiplier)
  probability <- as.numeric(plogis(X %*% c(inner$intercept, inner$slopes)))

  auc_values <- numeric(auc_replicates)
  for (r in seq_len(auc_replicates)) {
    set.seed(iteration * 100 + r + calibration_seed)
    y <- rbinom(N, 1, probability)
    auc_values[r] <- calculate_auc(y, probability)
  }
  average_auc <- mean(auc_values)
  history[[iteration]] <- data.frame(
    iteration = iteration, multiplier = multiplier,
    intercept = inner$intercept,
    beta1 = unname(inner$slopes[1]),
    beta2 = unname(inner$slopes[2]),
    beta3 = unname(inner$slopes[3]),
    expected_prevalence = inner$prevalence,
    average_auc = average_auc, auc_sd = sd(auc_values)
  )

  converged <- abs(inner$prevalence - target_prevalence) <= prevalence_tolerance &&
    abs(average_auc - target_auc) <= auc_tolerance
  if (converged) break
  if (average_auc < target_auc) lo <- multiplier else hi <- multiplier
}
trace <- do.call(rbind, history[seq_len(iteration)])

# 6. Retain the coefficients and calibration status ----------------------------

# These are the final candidate values. Use them as calibrated coefficients
# only when calibration$converged is TRUE.
coefficients <- c(
  intercept = inner$intercept,
  beta1 = unname(inner$slopes[1]),
  beta2 = unname(inner$slopes[2]),
  beta3 = unname(inner$slopes[3])
)
calibration <- list(
  converged = converged,
  coefficients = coefficients,
  multiplier = multiplier,
  target_prevalence = target_prevalence,
  target_auc = target_auc,
  expected_prevalence = inner$prevalence,
  average_auc = average_auc,
  auc_sd = sd(auc_values),
  iterations = iteration,
  covariates = covariates,
  N = N,
  trace = trace
)

cat("\nCalibration targets and achieved values:\n")
print(data.frame(
  metric = c("Expected prevalence", "Average simulated AUC"),
  target = c(target_prevalence, target_auc),
  achieved = c(inner$prevalence, average_auc)
), row.names = FALSE)
cat("\nBoth tolerances met:", converged, "\n")
if (converged) {
  cat("Converged after", iteration, "outer iterations.\n")
} else {
  cat("Stopped after", iteration, "outer iterations without meeting both targets.\n")
}
# Each outer iteration proposes a slope multiplier and runs an inner
# intercept search. calibration$iterations records the outer iteration count.
cat("\nIntercept and slopes:\n")
print(coefficients)

# 7. Plot the calibration trace -----------------------------------------------

# These plots describe the calibration search, not an MCMC fit.
# Dashed lines mark prevalence/AUC targets and the final intercept value.
# The AUC band is mean +/- one SD across the five outcome draws (or the
# chosen auc_replicates); it is not a confidence interval.
plot_convergence <- function(calibration) {
  history <- calibration$trace
  step <- history$iteration
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 3, 0))

  plot(step, history$expected_prevalence, type = "o", pch = 16,
       cex = 0.5, col = "#0077BB",
       ylim = range(history$expected_prevalence, calibration$target_prevalence),
       xlab = "Outer iteration", ylab = "Expected prevalence",
       main = "Prevalence")
  abline(h = calibration$target_prevalence, lty = 2, col = "grey30")

  lower <- history$average_auc - history$auc_sd
  upper <- history$average_auc + history$auc_sd
  plot(step, history$average_auc, type = "n",
       ylim = range(lower, upper, calibration$target_auc),
       xlab = "Outer iteration", ylab = "Average AUC", main = "AUC (band: +/- 1 SD)")
  if (length(step) > 1L) {
    polygon(c(step, rev(step)), c(lower, rev(upper)),
            col = adjustcolor("#E63946", alpha.f = 0.15), border = NA)
  } else {
    segments(step, lower, step, upper, col = "#E63946")
  }
  lines(step, history$average_auc, type = "o", pch = 16, cex = 0.5,
        col = "#E63946")
  abline(h = calibration$target_auc, lty = 2, col = "grey30")

  plot(step, history$intercept, type = "o", pch = 16, cex = 0.5,
       xlab = "Outer iteration", ylab = "Intercept", main = "Intercept")
  abline(h = tail(history$intercept, 1), lty = 2, col = "grey50")

  slope_colours <- c("#E63946", "#D89000", "#0077BB")
  matplot(step, as.matrix(history[c("beta1", "beta2", "beta3")]),
          type = "o", pch = 16, cex = 0.4, lty = 1,
          col = slope_colours, xlab = "Outer iteration", ylab = "Scaled slope",
          main = "Slopes")
  legend("topright", legend = c("beta1", "beta2", "beta3"),
         col = slope_colours, lty = 1, bty = "n", cex = 0.8)

  status <- if (calibration$converged) "Targets reached" else "Targets not reached"
  mtext(sprintf("%s after %d outer iterations", status, calibration$iterations),
        outer = TRUE, line = 1, font = 2)
  invisible(NULL)
}

# Display the trace even if the iteration limit was reached without convergence.
# To display it again later, run plot_convergence(calibration).
if (show_convergence_plots) plot_convergence(calibration)

# 8. Generate a dataset after successful calibration --------------------------

# The target prevalence concerns mean(probability). A new binary outcome draw
# has sampling variation, so its observed prevalence and AUC need not equal
# the targets. Calibration applies to this covariate distribution and design.
simulated_data <- NULL
if (converged) {
  set.seed(data_seed)
  simulated_data <- data.frame(
    covariates,
    probability = probability,
    y = rbinom(N, 1, probability)
  )
  cat("\nGenerated", N, "observations in simulated_data.\n")
  cat("Observed prevalence:", mean(simulated_data$y), "\n")
  # A one-class final draw has undefined AUC but is still a valid random draw.
  observed_auc <- if (length(unique(simulated_data$y)) == 2L) {
    calculate_auc(simulated_data$y, simulated_data$probability)
  } else NA_real_
  cat("Observed AUC:", observed_auc, "\n")
} else {
  warning(
    "The search did not meet both tolerances. No dataset was generated. ",
    "The reported coefficients are final candidates, not a successful calibration. ",
    "Inspect trace and consider the search bounds, AUC variability and whether ",
    "the targets are attainable for these covariates and relative slopes."
  )
}

# The objects remain available after sourcing this script. For example:
# glm(y ~ x1 + x2 + x3, data = simulated_data, family = binomial())
# Run this example only after successful calibration. Its fitted coefficients
# are estimates and will differ from the generating coefficients above.
