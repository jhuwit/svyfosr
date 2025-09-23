#' Survey-weighted Functional on Scalar Regression
#'
#' @description A high-level wrapper for running survey-weighted
#'   function-on-scalar regression with bootstrap inference.
#'
#' @param data A data frame with functional outcome columns, predictors, weights, etc.
#' @param family Outcome distribution family (e.g., "gaussian", "binomial").
#' @param estimator_type Type of estimator: "weighted" or "unweighted".
#' @param boot_type Bootstrap method: "BRR", "Rao-Wu-Yue-Beaumont", "weighted", "unweighted".
#' @param num_boots Number of bootstrap replicates.
#' @param nknots_min Minimum number of knots for smoothing (optional).
#' @param seed Random seed for reproducibility.
#' @param model_formula Formula with functional outcome on predictors, e.g. \code{Y ~ X1 + X2}.
#' @param ... Additional arguments passed to helpers.
#'
#' @importFrom stats model.matrix as.formula
#' @return A list with components:
#'   \item{betaHat}{Smoothed coefficient functions}
#'   \item{cis}{Bootstrap confidence intervals}
#'   \item{boots}{Raw bootstrap draws of coefficients}
#'
#' @export
svyfui <- function(data,
                   family = "gaussian",
                   boot_type,
                   weighted = TRUE,
                   num_boots = 500,
                   nknots_min = NULL,
                   seed = 2025,
                   model_formula = as.formula("Y ~ X1 + X2"),
                   ...) {
  # step 0: extract matrices from data
  out_index <- grep(paste0("^", model_formula[2]), names(data))
  Y_mat <- as.matrix(data[, out_index])
  # add intercept to X
  X_base <- stats::model.matrix(stats::as.formula(paste0("~", model_formula[3])), data = data)
  if(!(("weight") %in% colnames(data)) & weighted){
    stop("data must contain a 'weight' column for weighted estimation")
  }
  # step 1: get betatilde
  betaTilde <- get_betatilde(X_mat = X_bas,
                             Y_mat = Y_mat,
                             w = data$weights,
                             family = family)

  # step 2: Smooth to get betaHat
  betaHat <- get_betahat(betaTilde = betaTilde, nknots_min = nknots_min)

  # step 3: Run bootstrap
  boots <- run_boots_fast(
    data = data,
    boot_type = boot_type,
    family = family,
    num_boots = num_boots,
    seed = seed,
    X_base = X_base,
    Y_mat = Y_mat,
    L = ncol(betaHat),
    ...
  )
  ### start here
  # Step 4: Compute confidence intervals
  cis <- get_cis(betaTilde_boot = boots,
                 betaHat = betaHat,
                 L = ncol(betaHat),
                 nknots_min = nknots_min)

  out <- list(betaHat = betaHat,
              boots = boots,
              cis = cis)
  class(out) <- "svyfui"
  return(out)
}
