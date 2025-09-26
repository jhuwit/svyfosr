#' Batched WLS
#'
#' @description
#' Computes fast batched (optionally weighted) LM with matrix Y and matrix X using QR decomposition.
#' Typically used inside \code{\link{svyfui}}
#' @param X Design matrix of predictors (including intercept if desired)
#' @param Y Matrix of functional outcomes (rows = subjects, columns = functional points)
#' @param w Optional vector of weights for weighted regression
#'
#' @importFrom assertthat assert_that
#' @return A matrix of coefficient estimates \eqn{B} with
#'   predictors in rows and time/functional points in columns.
#'
#' @keywords internal
#' @noRd
lm_wls_multi <- function(X, Y, w = NULL) {
  assertthat::assert_that(is.matrix(X), is.matrix(Y), msg = "X and Y must be matrices")
  n <- nrow(X); p <- ncol(X); B <- ncol(Y)
  assertthat::assert_that(nrow(Y) == n, msg = "X and Y must have the same number of rows")

  # ---- process weights ----
  if (is.null(w)) w = rep(1, n)   # default: equal weights

  assertthat::assert_that(is.numeric(w), length(w) == n, msg = "w must be numeric with length n")

  # ---- weighted WLS ----
  W_half <- sqrt(w)
  Xw <- X * W_half
  Yw <- Y * W_half
  coef_mat <- qr.coef(qr(Xw), Yw)

  coef_mat
}

#' Batched GLM
#' @description
#' Computes fast batched (optionally weighted) GLM with matrix Y and matrix X using IRLS
#' Typically used inside \code{\link{svyfui}}
#' @param X Design matrix of predictors (including intercept if desired)
#' @param Y Matrix of functional outcomes (rows = subjects, columns = functional points)
#' @param w Optional vector of weights for weighted regression
#' @param family A GLM family object (e.g., from \code{binomial()}, \code{poisson()}), can also be character
#' @param add_intercept Whether to add an intercept column to X (default FALSE)
#' @param start Optional starting value for coefficients (length p or p x B)
#' @param maxit Maximum number of IRLS iterations (default 50)
#' @param tol Convergence tolerance (default 1e-8)
#' @param ridge Ridge penalty to ensure numerical stability (default 1e-8)
#' @param verbose Whether to print iteration info (default FALSE)
#'
#' @importFrom assertthat assert_that
#' @return A matrix of coefficient estimates \eqn{B} with
#'   predictors in rows and time/functional points in columns.
#'
#' @keywords internal
#' @noRd
glm_batch_multiY <- function(
    X, Y, w = NULL, data = NULL,       # data optional if w is a column name
    family,
    add_intercept = FALSE,
    start = NULL,
    maxit = 50, tol = 1e-8, ridge = 1e-8,
    verbose = FALSE
) {
  assertthat::assert_that(is.matrix(X), is.matrix(Y), msg = "X and Y must be matrices")
  n <- nrow(X); p <- ncol(X); B <- ncol(Y)
  assertthat::assert_that(nrow(Y) == n, msg = "X and Y must have the same number of rows")

  # ---- process weights ----
  if (is.null(w)) w <- rep(1, n)   # default: equal weights


  assertthat::assert_that(is.numeric(w), length(w) == n, msg = "w must be numeric with length n")

  if (add_intercept) X <- cbind(Intercept = 1, X)
  p <- ncol(X)

  # ---- family functions ----
  linkinv    <- family$linkinv
  mu_eta_fun <- family$mu.eta
  variance   <- family$variance
  dev_resids <- family$dev.resids
  eps <- .Machine$double.eps^(2/3)

  # ---- IRLS init ----
  Beta <- if (!is.null(start)) {
    stopifnot(length(start) == p)
    matrix(start, p, B)
  } else matrix(0, p, B)

  Xt <- t(X)
  converged <- rep(FALSE, B)

  for (it in seq_len(maxit)) {
    Eta <- X %*% Beta
    Eta <- pmin(pmax(Eta, -35), 35)
    Mu  <- linkinv(Eta)

    mu_eta <- mu_eta_fun(Eta); mu_eta[abs(mu_eta) < eps] <- eps
    VarMu  <- variance(Mu); VarMu[VarMu < eps] <- eps
    if (is.vector(VarMu) && length(VarMu) == n * B) VarMu <- matrix(VarMu, n, B)

    Ww <- (mu_eta^2 / VarMu) * matrix(w, n, B)
    Z  <- Eta + (Y - Mu) / mu_eta

    RHS <- Xt %*% (Ww * Z)
    step_max <- rep(0, B)
    active <- which(!converged)
    if (!length(active)) break

    s_all <- (mu_eta^2 / VarMu)
    for (b in active) {
      wb <- w * s_all[, b]
      if (!any(is.finite(wb)) || sum(wb) < eps) next
      Xw <- X * wb
      H  <- Xt %*% Xw
      diag(H) <- diag(H) + ridge
      beta_new <- tryCatch({
        Rchol <- chol(H)
        backsolve(Rchol, forwardsolve(t(Rchol), RHS[, b]))
      }, error = function(e) solve(H, RHS[, b], tol = 1e-12))
      step_max[b] <- max(abs(beta_new - Beta[, b]))
      Beta[, b]   <- beta_new
    }

    converged[(!converged) & (step_max < tol)] <- TRUE
    if (all(converged)) break
  }

  out <- list(coef = Beta, iter = it, converged = converged)
  return(out)
}

