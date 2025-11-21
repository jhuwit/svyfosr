# Survey-weighted Functional on Scalar Regression

A high-level wrapper for running survey-weighted function-on-scalar
regression with bootstrap inference.

## Usage

``` r
svyfui(
  formula,
  data,
  weights = NULL,
  family = gaussian(),
  boot_type = "weighted",
  num_boots = 500,
  nknots_min = NULL,
  nknots_min_fpca = NULL,
  seed = 2025,
  conf_level_pw = 0.95,
  conf_level_joint = 0.95,
  verbose = TRUE,
  parallel = FALSE,
  n_cores = NULL,
  ...
)
```

## Arguments

- formula:

  Formula with functional outcome on predictors, e.g. `Y ~ X1 + X2`.

- data:

  A data frame with functional outcome columns, predictors, weights,
  etc.

- weights:

  Optional bare column name for weights or external weight vector

- family:

  Outcome distribution family (e.g., "gaussian", "binomial").

- boot_type:

  Bootstrap method: "BRR", "Rao-Wu-Yue-Beaumont", "weighted",
  "unweighted".

- num_boots:

  Number of bootstrap replicates.

- nknots_min:

  Minimum number of knots for smoothing (optional).

- nknots_min_fpca:

  Minimum number of knots for FPCA (optional).

- seed:

  Random seed for reproducibility.

- conf_level_pw:

  Confidence level for pointwise confidence intervals (default 0.95).

- conf_level_joint:

  Confidence level for joint confidence intervals (default 0.95).

- verbose:

  Whether to print messages about progress

- parallel:

  Whether to run bootstrap in parallel (default FALSE).

- n_cores:

  Number of cores to use if parallel = TRUE

- ...:

  Additional arguments passed to helpers.

## Value

A list with components:

- betaHat:

  Smoothed coefficient functions

- cis:

  Bootstrap confidence intervals

- boots:

  Raw bootstrap draws of coefficients

- tidy_df:

  Tidy data frame for plotting with columns `l` (functional domain),
  `beta_hat` (estimate), `lower_pw` (pointwise lower CI), `upper_pw`
  (pointwise upper CI), `lower_joint` (joint lower CI), `upper_joint`
  (joint upper CI), and `var_name` (variable name from regression)

## Examples

``` r
# \donttest{
   fit <- svyfui(Y ~ X, data = svyfosr::sample_df, weights = weight)
#> Estimating coefficients
#> Smoothing coefficients
#> Bootstrapping
#> Obtaining pointwise and joint confidence intervals
#> Completed!
   plot(fit)

# }
```
