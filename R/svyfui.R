#' Survey-weighted Functional on Scalar Regression
#'
#' @description A high-level wrapper for running survey-weighted
#'   function-on-scalar regression with bootstrap inference.
#'
#' @param formula Formula with functional outcome on predictors, e.g. \code{Y ~ X1 + X2}.
#' @param data A data frame with functional outcome columns, predictors, weights, etc.
#' @param weights Optional bare column name for weights or external weight vector
#' @param family Outcome distribution family (e.g., "gaussian", "binomial").
#' @param boot_type Bootstrap method: "BRR", "Rao-Wu-Yue-Beaumont", "weighted", "unweighted".
#' @param num_boots Number of bootstrap replicates.
#' @param nknots_min Minimum number of knots for smoothing (optional).
#' @param seed Random seed for reproducibility.
#' @param conf_level_joint Confidence level for joint confidence intervals (default 0.95).
#' @param conf_level_pw Confidence level for pointwise confidence intervals (default 0.95).
#' @param verbose Whether to print messages about progress
#' @param ... Additional arguments passed to helpers.
#'
#' @importFrom stats model.matrix as.formula model.frame qnorm model.response model.weights var gaussian
#' @importFrom tibble tibble
#' @importFrom dplyr mutate
#' @importFrom purrr map_dfr
#' @return A list with components:
#'   \item{betaHat}{Smoothed coefficient functions}
#'   \item{cis}{Bootstrap confidence intervals}
#'   \item{boots}{Raw bootstrap draws of coefficients}
#'
#' @export
#' @examples
#' \donttest{
#'    fit <- svyfui(Y ~ X, data = svyfosr::sample_df, weights = weight)
#'    plot(fit)
#' }
svyfui <- function(formula,
                   data,
                   weights = NULL,
                   family = gaussian(),
                   boot_type = "weighted",
                   num_boots = 500,
                   nknots_min = NULL,
                   seed = 2025,
                   conf_level_pw = 0.95,
                   conf_level_joint = 0.95,
                   verbose = TRUE,
                   ...) {
  beta_hat = l = lower_joint = lower_pw = upper_joint = upper_pw = NULL
  psu = strata = weight = . = NULL
  rm(list = c("beta_hat", "l", "lower_joint", "lower_pw", "upper_joint", "upper_pw",
              "psu", "strata", "weight", "."))

  # deal with weights (copied from glm fn)
  # first check character stuff - should prob move helper fn outside
  check_var(weights)

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

  if (verbose) message("Estimating coefficients")
  # step 1: get betatilde
  betaTilde <- get_betatilde(X_mat = X_base,
                             Y_mat = Y_mat,
                             w = w,
                             family = family)

  if (verbose) message("Smoothing coefficients")
  # step 2: Smooth to get betaHat
  betaHat <- get_betahat(betaTilde = betaTilde, nknots_min = nknots_min)

  if (verbose) message("Bootstrapping")
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
    ...
  )

  if (verbose) message("Obtaining pointwise and joint confidence intervals")
  cis <- get_cis(betaTilde_boot = boots,
                 betaHat = betaHat,
                 L = ncol(betaHat),
                 nknots_min = nknots_min,
                 conf_level = conf_level_joint)

  ## step 5: make into tidy data frame
  num_var = nrow(betaHat)
  m = stats::qnorm((1 + conf_level_pw) / 2) # pw ci multiplier

  plt_df = purrr::map_dfr(
    .x = 1:num_var,
    .f = function(r) {
      tibble::tibble(
        l = 1:length(betaHat[r, ]),
        beta_hat = betaHat[r, ],
        lower_pw = betaHat[r, ] - m * sqrt(diag(cis$betaHat.var[ , , r])),
        upper_pw = betaHat[r, ] + m * sqrt(diag(cis$betaHat.var[ , ,r ])),
        lower_joint = betaHat[r, ] - cis$qn[r] * sqrt(diag(cis$betaHat.var[, , r])),
        upper_joint = betaHat[r, ] + cis$qn[r] * sqrt(diag(cis$betaHat.var[, , r]))
      ) %>%
        dplyr::mutate(var_name = betaHat %>% rownames() %>% .[r])
    }
  )

  out <- list(betaHat = betaHat,
              boots = boots,
              cis = cis,
              tidy_df = plt_df)
  class(out) <- "svyfui"

  if(verbose) message("Completed!")
  return(out)
}
