# Example functional survey data

A dataset included in svyfosr to illustrate usage.

## Usage

``` r
sample_df
```

## Format

A data frame with 975 rows and 8 columns

- id:

  Numeric. Subject identifier

- X:

  Numeric. X variable (sample covariate)

- strata:

  Numeric. Strata variable for survey design, ranges from 1 to 10

- psu:

  Character. PSU (cluster) variable for survey design. Two unique values
  per stratum

- weight:

  Numeric. Survey weight, inverse of sampling probability

- p_stage1:

  Numeric. Stage 1 sampling probability. Ranges from 0 to 1

- p_stage2:

  Numeric. Stage 2 sampling probability. Ranges from 0 to 1

- Y:

  Matrix with N = 975 rows and L = 50 columns. Functional response
  variable observed on a regular grid

## Source

Simulated data
