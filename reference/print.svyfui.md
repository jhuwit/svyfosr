# Print an `svyfui` object

Print an `svyfui` object

## Usage

``` r
# S3 method for class 'svyfui'
print(x, ...)
```

## Arguments

- x:

  an object of class
  [svyfui](https://jhuwit.github.io/svyfosr/reference/svyfui.md)

- ...:

  additional arguments to pass to other
  [print](https://rdrr.io/r/base/print.html) methods

## Value

The object, invisibly

## Examples

``` r
# \donttest{
   fit <- svyfui(Y ~ X, data = svyfosr::sample_df, weights = weight)
#> Estimating coefficients
#> Smoothing coefficients
#> Bootstrapping
#> Obtaining pointwise and joint confidence intervals
#> Completed!
   print(fit)
#> Survey FUI object
#> -----------------
#> Number of bootstrap replicates: 500 
#> # A tibble: 100 × 7
#>        l beta_hat lower_pw upper_pw lower_joint upper_joint var_name   
#>    <int>    <dbl>    <dbl>    <dbl>       <dbl>       <dbl> <chr>      
#>  1     1   -0.179   -0.195   -0.163      -0.202      -0.156 (Intercept)
#>  2     2   -0.191   -0.205   -0.178      -0.211      -0.172 (Intercept)
#>  3     3   -0.203   -0.215   -0.192      -0.220      -0.187 (Intercept)
#>  4     4   -0.214   -0.225   -0.204      -0.230      -0.199 (Intercept)
#>  5     5   -0.224   -0.236   -0.212      -0.241      -0.207 (Intercept)
#>  6     6   -0.231   -0.243   -0.220      -0.248      -0.215 (Intercept)
#>  7     7   -0.237   -0.248   -0.225      -0.253      -0.220 (Intercept)
#>  8     8   -0.240   -0.253   -0.228      -0.258      -0.223 (Intercept)
#>  9     9   -0.242   -0.254   -0.229      -0.260      -0.224 (Intercept)
#> 10    10   -0.242   -0.254   -0.229      -0.260      -0.224 (Intercept)
#> # ℹ 90 more rows
# }
```
