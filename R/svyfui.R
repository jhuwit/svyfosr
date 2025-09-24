#' Survey-weighted Functional on Scalar Regression
#'
#' @description A high-level wrapper for running survey-weighted
#'   function-on-scalar regression with bootstrap inference.
#'
#' @param model_formula Formula with functional outcome on predictors, e.g. \code{Y ~ X1 + X2}.
#' @param data A data frame with functional outcome columns, predictors, weights, etc.
#' @param weights Optional column name for weights or external weight vector
#' @param family Outcome distribution family (e.g., "gaussian", "binomial").
#' @param estimator_type Type of estimator: "weighted" or "unweighted".
#' @param boot_type Bootstrap method: "BRR", "Rao-Wu-Yue-Beaumont", "weighted", "unweighted".
#' @param num_boots Number of bootstrap replicates.
#' @param nknots_min Minimum number of knots for smoothing (optional).
#' @param seed Random seed for reproducibility.
#' @param ... Additional arguments passed to helpers.
#'
#' @importFrom stats model.matrix as.formula
#' @importFrom rlang is_symbol as_name
#' @return A list with components:
#'   \item{betaHat}{Smoothed coefficient functions}
#'   \item{cis}{Bootstrap confidence intervals}
#'   \item{boots}{Raw bootstrap draws of coefficients}
#'
#' @export
svyfui <- function(formula,
                   data,
                   weights = NULL,
                   family = gaussian(),
                   boot_type = "weighted",
                   num_boots = 500,
                   nknots_min = NULL,
                   seed = 2025,
                   ...) {

  # deal with weights (copied from glm fn)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action"),
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  Y_mat <- as.matrix(model.response(mf))  # lhs as matrix
  X_base <- stats::model.matrix(attr(mf, "terms"), mf)

  # --- weights ---
  w <- as.vector(model.weights(mf))
  if (is.null(w)) w <- rep(1, nrow(mf))  # default to uniform weights


  # step 1: get betatilde
  betaTilde <- get_betatilde(X_mat = X_base,
                             Y_mat = Y_mat,
                             w = w,
                             family = family)

  # step 2: Smooth to get betaHat
  betaHat <- get_betahat(betaTilde = betaTilde, nknots_min = nknots_min)

  # step 3: Run bootstrap
  boots <- run_boots(
    data = data,
    boot_type = boot_type,
    weights  = w,
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
