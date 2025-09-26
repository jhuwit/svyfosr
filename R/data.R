#' Example functional survey data
#'
#' A dataset included in svyfosr to illustrate usage.
#'
#' @format A data frame with 975 rows and 8 columns
#' \describe{
#'   \item{id}{Numeric. Subject identifier}
#'   \item{X}{Numeric. X variable (sample covariate)}
#'   \item{strata}{Numeric. Strata variable for survey design, ranges from 1 to 10}
#'   \item{psu}{Character. PSU (cluster) variable for survey design. Two unique values per stratum}
#'   \item{weight}{Numeric. Survey weight, inverse of sampling probability}
#'   \item{p_stage1}{Numeric. Stage 1 sampling probability. Ranges from 0 to 1}
#'   \item{p_stage2}{Numeric. Stage 2 sampling probability. Ranges from 0 to 1}
#'   \item{Y}{Matrix with N = 975 rows and L = 50 columns. Functional response variable observed on a regular grid}
#' }
#' @source Simulated data
"sample_df"
