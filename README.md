# Data generation for the logistic regression model

These files calibrate a logistic model to a target prevalence and AUC
(c-statistic). They provide a standalone data generator and a simulation
comparison of three calibration algorithms.

## Files

| File | Purpose |
|---|---|
| `generate_logistic_data.R` | Runs Algorithm 3 once, retains the calibrated coefficients and generates one dataset after successful calibration. |
| `prevalence_auc_simulation.Rmd` | Compares Algorithms 1, 2 and 3 across four scenarios, with summaries and plots. |
| `distribution_of_predicted_probabilities_all_algorithms.R` | Distribution of predicted probabilities for Algorithms 1, 2 and 3. |
| `distribution_of_predicted_probabilities_under_perfect_separation.R`| Distribution of predicted probabilities when AUC is 1. |
| `effect_of_priors_on_predicted_probabilities.R` | Effect of prior distribution on predicted probabilites. |

The model uses an intercept and three slopes. Calibration multiplies the
initial slopes by one common scaling factor, preserving their relative sizes.

## 1. Install the required packages

Run this once if the packages are not already installed:

```r
install.packages(c("pROC", "knitr", "ggplot2", "dplyr", "tidyr"))
```

The standalone R script requires only `pROC`; its convergence plots use base R.
The R Markdown document uses all five packages. No external data files are
required. Run the workflows separately, preferably in fresh R sessions,
because both create objects in the working environment.

## 2. Generate one dataset with Algorithm 3

Open `generate_logistic_data.R` and set the two targets:

```r
target_prevalence <- 0.8
target_auc <- 0.85
```

From the directory containing the file, run:

```r
source("generate_logistic_data.R")
```

The default design has 1,000 observations and mixed covariates: `x1` and `x3`
are standard normal, and `x2` is Bernoulli with probability 0.5. The initial
relative slopes are `log(2)`, `log(1.2)` and `log(2.8)`.

To change the predictors, edit the covariate-generation section. Predictors
may be continuous, binary or mixed. Keep the names `x1`, `x2`, `x3` and their
ordering consistent with `base_slopes`. Recalibrate after changing their
distributions or values.

### How calibration works

1. Propose a common slope multiplier.
2. Find the intercept that gives the target expected prevalence on the fixed covariates.
3. Generate five outcome vectors using those same coefficients and average their AUCs.
4. Update the multiplier and repeat until both tolerances are met or the iteration limit is reached.

The five outcome draws stabilise the AUC estimate. They use the same intercept
and slopes within an iteration; the intercept itself is not averaged.
This is one calibration run, with no 100-run performance study.

Defaults are 0.0001 for both tolerances, up to 1,000 iterations in each search,
intercept bounds of [-20, 20], and multiplier bounds of [0, 20]. These settings
can be edited near the top of the script.

### Results

The script prints the achieved calibration values, intercept, slopes,
convergence status and number of outer iterations. It retains:

| Object | Contents |
|---|---|
| `coefficients` | Final intercept and slopes `beta1`, `beta2`, `beta3`. |
| `calibration` | Targets, achieved values, convergence flag, multiplier, iteration count, covariates and trace. |
| `simulated_data` | Covariates, outcome probabilities and a binary outcome `y`; `NULL` if calibration did not converge. |
| `trace` | Intercept, slopes, multiplier, expected prevalence and AUC statistics at each outer iteration. |

Inspect the results with:

```r
calibration$converged
calibration$iterations
coefficients
head(simulated_data)
```

Only treat the returned coefficients as successfully calibrated when
`calibration$converged` is `TRUE`. Otherwise, they are the final search
candidates, and no dataset is generated.

After successful calibration, a logistic model can be fitted to the data:

```r
if (calibration$converged) {
  fitted_model <- glm(
    y ~ x1 + x2 + x3,
    data = simulated_data,
    family = binomial()
  )
  summary(fitted_model)
}
```

Fitted coefficients are estimates and can differ from the generating values.

### Convergence plots

Four plots appear automatically: expected prevalence, average AUC, intercept
and scaled slopes. Display them again with:

```r
plot_convergence(calibration)
```

Set `show_convergence_plots <- FALSE` to disable automatic plotting. Dashed
lines show the prevalence and AUC targets, and the final intercept. The AUC
band is the mean plus or minus one SD across the outcome draws, not a
confidence interval. These plots describe calibration, not MCMC convergence.

## 3. Compare the three algorithms

Open `prevalence_auc_simulation.Rmd` in RStudio and run its code chunks in
order. No knitting or HTML generation is required.

- **Algorithm 1:** adjusts the intercept to target prevalence, with the slope multiplier fixed at 1.
- **Algorithm 2:** adjusts the slope multiplier to target AUC, with the intercept fixed at zero.
- **Algorithm 3:** uses nested bisection to target both quantities.

The document runs all four scenarios:

| Covariates | Target prevalence | Target AUC |
|---|---:|---:|
| Binary | 0.1 | 0.85 |
| Binary | 0.8 | 0.85 |
| Continuous | 0.1 | 0.85 |
| Continuous | 0.8 | 0.85 |

Binary predictors use Bernoulli probabilities 0.5, 0.5 and 0.7. Continuous
predictors are standard normal. Each design has 5,000 observations. Covariates
are fixed across runs, and the two prevalence scenarios reuse the same design
within each covariate type.

The default `n.sim <- 100` gives 100 runs per algorithm per scenario. Use
`n.sim <- 2` for a small initial check. This setting is separate from the five
outcome draws used within each AUC-calibration iteration.

The document retains:

- `results_by_scenario`: a named list of the four scenarios' results.
- `results`: combined results with scenario labels and replicate IDs.
- `summary_results`: prevalence/AUC bias and SD by scenario and method.
- `conv`: the nested-search trace for one selected scenario.

Four comparison plots display the scenario results. Convergence plots are
provided only for `trace_scenario`, which defaults to `"binary_prev08"`.
Change that setting to another scenario name in the settings table if needed.
Algorithm 2 does not use the prevalence target, so its results are identical
across the two prevalence targets for a given covariate design and seed sequence.

## 4. Interpret the targets and outputs

Expected prevalence is the mean of the model probabilities. Calibration AUC
is the average AUC across simulated outcome draws. The final dataset is a new
random draw, so its realised prevalence and AUC need not equal those calibration
values exactly.

AUC calibration remains stochastic, and the requested targets may not be
attainable for every covariate design, set of relative slopes or search interval.
Check convergence and achieved values before using the coefficients. A stable
coefficient trace alone does not establish that both targets were reached.

The scripts retain pROC's automatic AUC direction selection and pass predictions
as numeric vectors to avoid the deprecated matrix-predictor warning.

## 5. Reproducibility and saving

The standalone script separates covariate, calibration and final-data seeds.
The comparison document uses fixed covariate and replicate seeds. Reproducing
results requires the same settings, covariates and compatible R/package versions.

Results remain in the R session. To save the standalone output explicitly:

```r
saveRDS(
  list(calibration = calibration, data = simulated_data),
  file = "logistic_data_generation.rds"
)
```


