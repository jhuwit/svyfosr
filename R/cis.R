#' Obtaining CIs for functional regression coefficients
#'
#' @description
#' Computes confidence intervals from bootstrap coefficient array
#' Implemented in \code{\link{svyfui}}
#' @param betaTilde_boot An array of bootstrap coefficient estimates of dimension p x L x B
#' @param betaHat A matrix of smoothed coefficient estimates (p x L)
#' @param smooth_for_ci Whether to smooth bootstrap estimates for CI calculation (default TRUE)
#' @param smooth_for_variance Whether to smooth bootstrap estimates for variance calculation (default TRUE)
#' @param L Number of functional points
#' @param nknots_min Minimum number of knots for smoothing (optional)
#' @param nknots_min_cov Minimum number of knots for covariance smoothing (default 35)
#' @param mult_fac Multiplicative factor for variance inflation (default 1.2)
#' @importFrom refund fpca.face
#' @importFrom mgcv gam
#' @importFrom stats rnorm quantile
#' @importFrom Rfast colVars
#' @return a list with components:
#' betaHat, the smoothed coefficient estimate matrix (p x L)
#' betaHat.var the variance estimate array (L x L x p)
#' qn the critical values for joint CIs for each coefficient (length p)
#'
#'
#' @keywords internal
#' @noRd
get_cis = function(betaTilde_boot,
                   betaHat,
                   smooth_for_ci = TRUE,
                   smooth_for_variance = TRUE,
                   L,
                   nknots_min = NULL,
                   nknots_min_cov = 35,
                   mult_fac = 1.2) {
  argvals = 1:L
  B <- dim(betaTilde_boot)[3]
  nknots <- if (is.null(nknots_min)) round(L / 2) else min(round(L / 2), nknots_min)

  nknots_cov <- ifelse(is.null(nknots_min_cov), 35, nknots_min_cov)
  nknots_fpca <- min(round(L / 2), 35)
  betaHat_boot <- array(NA, dim = c(nrow(betaHat), ncol(betaHat), ncol(betaTilde_boot[1, , ])))
  betaHat.var <- array(NA, dim = c(L, L, nrow(betaHat)))
  # smooth


  for (b in 1:B) {
    betaHat_boot[, , b] <- t(apply(betaTilde_boot[, , b], 1, function(x)
      mgcv::gam(x ~ s(
        argvals, bs = "tp", k = (nknots + 1)
      ), method = "GCV.Cp")$fitted.values))
  }

  for (r in 1:nrow(betaHat)) {
    if (smooth_for_variance) {
      betaHat.var[, , r] <- mult_fac * var(t(betaHat_boot[r, , ]))
    } else{
      betaHat.var[, , r] <- mult_fac * var(t(betaTilde_boot[r, , ]))
    }
  }

  N <- 10000
  qn <- rep(0, length = nrow(betaHat))
  set.seed(456)
  for (i in 1:length(qn)) {
    if (smooth_for_ci) {
      est_bs <- t(betaHat_boot[i, , ])
    } else{
      est_bs <- t(betaTilde_boot[i, , ])
    }
    fit_fpca <- suppressWarnings(refund::fpca.face(est_bs, knots = nknots_fpca)) # suppress sqrt(Eigen$values) NaNs
    ## extract estimated eigenfunctions/eigenvalues
    phi <- fit_fpca$efunctions
    lambda <- fit_fpca$evalues
    K <- length(fit_fpca$evalues)

    ## simulate random coefficients
    theta <- matrix(stats::rnorm(N * K), nrow = N, ncol = K) # generate independent standard normals
    if (K == 1) {
      theta <- theta * sqrt(lambda) # scale to have appropriate variance
      X_new <- tcrossprod(theta, phi) # simulate new functions
    } else{
      theta <- theta %*% diag(sqrt(lambda)) # scale to have appropriate variance
      X_new <- tcrossprod(theta, phi) # simulate new functions
    }
    x_sample <- X_new + t(fit_fpca$mu %o% rep(1, N)) # add back in the mean function
    Sigma_sd <- Rfast::colVars(x_sample, std = TRUE, na.rm = FALSE) # standard deviation: apply(x_sample, 2, sd)
    x_mean <- colMeans(est_bs)
    # x_sample: N x p matrix
    # x_mean: length-p vector
    # Sigma_sd: length-p vector

    # Vectorized computation
    z <- sweep(x_sample, 2, x_mean, "-")     # center
    z <- sweep(z, 2, Sigma_sd, "/")          # standardize
    un <- apply(abs(z), 1, max)              # row-wise max of absolute values

    qn[i] <- stats::quantile(un, 0.95)
  }

  return(list(
    betaHat = betaHat,
    betaHat.var = betaHat.var,
    qn = qn
  ))

}
