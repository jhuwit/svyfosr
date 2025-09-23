#' Batched WLS
#'
#' @description
#' Computes fast batched (optionally weighted) LM with matrix Y and matrix X using QR decomposition.
#' Typically used inside \code{\link{svyfui}}
#' @param X Design matrix of predictors (including intercept if desired)
#' @param Y Matrix of functional outcomes (rows = subjects, columns = functional points)
#' @param w Optional vector of weights for weighted regression
#'
#' @return A matrix of coefficient estimates \eqn{B} with
#'   predictors in rows and time/functional points in columns.
#'
#' @keywords internal
#' @noRd
lm_wls_multi <- function(X, Y, w = NULL) {
  stopifnot(is.matrix(X), is.matrix(Y))
  n <- nrow(X); p <- ncol(X); B <- ncol(Y)
  stopifnot(nrow(Y) == n)

  # X: (n × p), Y: (n × L), w: (n)
  if (is.null(w)) {
    coef_mat <- qr.coef(qr(X), Y)
  } else {
    W_half <- sqrt(w)
    Xw <- X * W_half
    Yw <- Y * W_half
    # Solve weighted least squares for all outcomes at once
    coef_mat <- qr.coef(qr(Xw), Yw)
  }
  # returns (p × L) coefficient matrix
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
#' @param offset Optional offset vector or matrix (length n or n x B)
#' @param start Optional starting value for coefficients (length p or p x B)
#' @param maxit Maximum number of IRLS iterations (default 50)
#' @param tol Convergence tolerance (default 1e-8)
#' @param ridge Ridge penalty to ensure numerical stability (default 1e-8)
#' @param return_se Whether to compute and return standard errors (default FALSE)
#' @param estimate_phi Whether to estimate dispersion parameter φ for non-Poisson/Binomial families (default FALSE)
#' @param verbose Whether to print iteration info (default FALSE)
#'
#' @return A matrix of coefficient estimates \eqn{B} with
#'   predictors in rows and time/functional points in columns.
#'
#' @keywords internal
#' @noRd
glm_batch_multiY <- function(
    X, Y, w = NULL,
    family,
    add_intercept = FALSE,
    offset = NULL,              # length-n or n x B (broadcasted if needed)
    start = NULL,               # optional p[+1]-vector to warm start all columns
    maxit = 50, tol = 1e-8, ridge = 1e-8,
    return_se = FALSE, estimate_phi = FALSE, verbose = FALSE
) {
  stopifnot(is.matrix(X), is.matrix(Y))
  n <- nrow(X); p <- ncol(X); B <- ncol(Y)
  stopifnot(nrow(Y) == n)

  if (is.null(w)) {
    w <- rep(1, n)   # no weights → all equal to 1
  }
  stopifnot(length(w) == n, is.numeric(w))

  if (add_intercept) X <- cbind(Intercept = 1, X)
  p <- ncol(X)

  # offsets: allow NULL (zeros), length-n (shared), n x 1 (shared), or n x B (per column)
  if (is.null(offset)) {
    offset <- matrix(0, n, B)
  } else if (is.vector(offset) && length(offset) == n) {
    offset <- matrix(offset, n, B)
  } else if (is.matrix(offset) && nrow(offset) == n && ncol(offset) == 1) {
    offset <- matrix(offset[,1], n, B)
  } else if (is.matrix(offset) && nrow(offset) == n && ncol(offset) == B) {
    # as provided
  } else {
    stop("offset must be NULL, length-n, n×1, or n×B")
  }

  # family bits
  linkinv    <- family$linkinv
  mu_eta_fun <- family$mu.eta
  variance   <- family$variance
  dev_resids <- family$dev.resids

  eps <- .Machine$double.eps^(2/3)

  # warm start
  if (!is.null(start)) {
    stopifnot(length(start) == p)
    Beta <- matrix(start, p, B)
  } else {
    Beta <- matrix(0, p, B)
  }

  # Precompute X' once
  Xt <- t(X)

  converged <- rep(FALSE, B)
  it <- 0L

  for (it in seq_len(maxit)) {
    # ETA & MU (n x B)
    Eta <- X %*% Beta
    Eta <- Eta + offset
    Eta <- pmin(pmax(Eta, -35), 35)
    Mu  <- linkinv(Eta)

    # dμ/dη and Var(μ) (n x B)
    mu_eta <- mu_eta_fun(Eta); mu_eta[abs(mu_eta) < eps] <- eps
    VarMu  <- variance(Mu);    VarMu[VarMu  < eps] <- eps
    if (is.vector(VarMu) && length(VarMu) == n * B) VarMu <- matrix(VarMu, n, B)  # expand vector to n × B

    # IRLS working weights and responses
    # Ww = w * (dμ/dη)^2 / Var(μ)  (n x B), z = η + (y - μ) / (dμ/dη)
    Ww <- (mu_eta^2 / VarMu)
    Ww <- Ww * matrix(w, n, B)      # broadcast w across columns
    Z  <- Eta + (Y - Mu) / mu_eta

    # RHS for all fits at once: X' (Ww * Z)  (p x B)
    RHS <- Xt %*% (Ww * Z)

    step_max <- rep(0, B)
    active <- which(!converged)
    if (!length(active)) break

    # For each column, H_b = X' diag(w * s_b) X where s_b = (dμ/dη)^2 / Var(μ)
    s_all <- (mu_eta^2 / VarMu)      # n x B (without prior weights)
    for (b in active) {
      wb <- w * s_all[, b]           # length-n
      if (!any(is.finite(wb)) || sum(wb) < eps) next

      # Efficient: crossprod(X, wb * X) without forming diag
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

    newly_conv <- (!converged) & (step_max < tol)
    converged[newly_conv] <- TRUE

    if (verbose) {
      cat(sprintf("Iter %d: active=%d, max step=%.3e\n",
                  it, length(active), if (length(active)) max(step_max[active]) else 0))
    }
    if (all(converged)) break
  }

  out <- list(coef = Beta, iter = it, converged = converged)

  if (return_se) {
    # recompute at solution
    Eta <- X %*% Beta + offset
    Eta <- pmin(pmax(Eta, -35), 35)
    Mu  <- linkinv(Eta)
    mu_eta <- mu_eta_fun(Eta); mu_eta[abs(mu_eta) < eps] <- eps
    VarMu  <- variance(Mu);    VarMu[VarMu  < eps] <- eps
    if (is.vector(VarMu) && length(VarMu) == n * B) VarMu <- matrix(VarMu, n, B)  # expand vector to n × B
    s_all  <- (mu_eta^2 / VarMu)

    # dispersion (φ) per column if needed
    fam_name <- tolower(family$family)
    phi <- rep(1, B)
    needs_phi <- estimate_phi && !grepl("binomial|poisson", fam_name)
    if (needs_phi) {
      for (b in seq_len(B)) {
        wb <- w
        wb[!is.finite(wb) | wb < 0] <- 0
        phi[b] <- sum(dev_resids(Y[, b], Mu[, b], wb)) / max(n - p, 1L)
      }
    }

    SE <- matrix(NA_real_, p, B, dimnames = list(colnames(X), colnames(Y)))
    vcov_list <- vector("list", B)
    for (b in seq_len(B)) {
      wb <- w * s_all[, b]
      if (!any(is.finite(wb)) || sum(wb) < eps) next
      H <- Xt %*% (X * wb)
      diag(H) <- diag(H) + ridge
      invH <- tryCatch({
        Rchol <- chol(H); chol2inv(Rchol)
      }, error = function(e) solve(H, tol = 1e-12))
      vc <- invH * phi[b]
      vcov_list[[b]] <- vc
      SE[, b] <- sqrt(pmax(diag(vc), 0))
    }
    out$se <- SE
    out$vcov <- vcov_list
    out$phi <- phi
  }

  out
}
