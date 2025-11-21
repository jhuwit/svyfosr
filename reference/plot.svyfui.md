# Plot method for `svyfui` objects

Plot method for `svyfui` objects

## Usage

``` r
# S3 method for class 'svyfui'
plot(
  x,
  include_joint_ci = TRUE,
  joint_fill = "#56B4E9FF",
  pw_fill = "#0072B2FF",
  joint_alpha = 0.3,
  pw_alpha = 0.3,
  ...
)
```

## Arguments

- x:

  An object of class
  [`svyfui`](https://jhuwit.github.io/svyfosr/reference/svyfui.md),
  usually returned by
  [`svyfui`](https://jhuwit.github.io/svyfosr/reference/svyfui.md).

- include_joint_ci:

  Logical; if TRUE, includes joint confidence intervals in the plot
  (default TRUE).

- joint_fill:

  Color for joint confidence interval fill (default "#56B4E9FF").

- pw_fill:

  Color for pointwise confidence interval fill (default "#0072B2FF").

- joint_alpha:

  Transparency level for joint confidence interval fill (default 0.3).

- pw_alpha:

  Transparency level for pointwise confidence interval fill (default
  0.3).

- ...:

  Additional arguments passed to plotting functions.

## Value

A `ggplot` object (recommended).

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
