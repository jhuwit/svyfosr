#' Estimate unsmoothed coefficient functions (beta tilde)
#'
#' @description
#' Computes raw estimates of coefficient functions across the functional domain,
#' before smoothing. Typically used inside \code{\link{svyfui}}
#' @param X_mat Design matrix of predictors (including intercept if desired)
#' @param Y_mat Matrix of functional outcomes (rows = subjects, columns = functional points)
#' @param w Optional vector of weights for weighted regression
#' @param family Error distribution family (e.g., "gaussian", "binomial", "poisson")
#' @return A matrix of coefficient estimates \eqn{B} with
#'   predictors in rows and time/functional points in columns.
#'
#' @keywords internal
#' @noRd
get_betatilde = function(X_mat, Y_mat, w = NULL, family = "gaussian"){
  if (is.character(family)) {
    family <- match.fun(family)()
  }

  if (family$family == "gaussian") {
    coef_mat = lm_wls_multi(X_mat, Y_mat, w)
  } else {
    coef_mat = glm_batch_multiY(X = X_mat, Y = Y_mat, w = w, family = family, return_se = FALSE,
                                add_intercept = FALSE)$coef
  }
  # Obtain betaTilde, fixed effects estimates
  colnames(coef_mat) <- seq_len(ncol(Y_mat))

  return(coef_mat)
}

#' Smooth coefficient estimates
#'
#' @description
#' Applies basis smoothing (e.g., B-splines) to unsmoothed coefficient
#' functions to produce final \eqn{\hat{\beta}(s)}.
#'
#' @param betaTilde A matrix of unsmoothed coefficient estimates.
#' @param nknots_min Minimum number of knots for the smoother (optional).
#' @importFrom mgcv gam
#' @return A smoothed matrix of coefficient estimates.
#'
#' @keywords internal
#' @noRd
get_betahat = function(betaTilde, nknots_min = NULL){
  L <- ncol(betaTilde)
  argvals <- seq_len(L)
  nknots <- if (is.null(nknots_min)) round(L / 2) else min(round(L / 2), nknots_min)
  betaHat <- t(apply(betaTilde, 1, function(x) mgcv::gam(x ~ s(argvals, bs = "tp",
                                                               k = (nknots + 1)),
                                                         method = "GCV.Cp")$fitted.values))
  rownames(betaHat) <- rownames(betaTilde)
  colnames(betaHat) <- argvals
  return(betaHat)
}
