#' Helper fn to sample data
#'
#' @description
#' Samples data for BRR
#' Implemented in \code{\link{run_boots}}
#' @param df A data frame with columns 'strata' and 'psu'
#' @importFrom dplyr group_by summarize
#'
#' @return A data frame with 1 psu selected per stratu,
#'
#' @keywords internal
#' @noRd
sample_data <- function(df) {
  strata = psu = NULL
  rm(list = c("strata", "psu"))
  df %>% dplyr::group_by(strata) %>%
    dplyr::summarize(psu = sample(psu, size = 1), .groups = "drop")
}

#' Bootstrapping
#'
#' @description
#' Computes bootstrap replicates for different bootstrap methods.
#' Implemented in \code{\link{svyfui}}
#' @param data A data frame with functional outcome columns, predictors, weights, etc.
#' @param X_base Design matrix of predictors (n x (p + 1)) including intercept
#' @param Y_mat Matrix of functional outcomes (n x L, each column is a time point)
#' @param weights External weight vector (length n)
#' @param boot_type Bootstrap method: "BRR", "RWYB", "weighted", "unweighted"
#' @param family Error distribution family (e.g., "gaussian", "binomial", "poisson")
#' @param num_boots Number of bootstrap replicates (default 500)
#' @param seed Random seed for reproducibility (default 2025)
#' @param samp_method_by_stage Sampling method by stage for RWYB bootstrap (e.g., c("PPSWR", "Poisson"))
#' @param parellel Logical whether bootstraping should be done in parallel
#' @param n_cores number of cores to be used for parallelization
#' @importFrom survey svydesign
#' @importFrom svrep make_rwyb_bootstrap_weights
#' @importFrom dplyr mutate group_by summarize row_number
#' @importFrom stats as.formula
#' @importFrom parallelly availableCores
#' @importFrom future plan multisession
#' @importFrom future.apply future_lapply
#' @importFrom assertthat assert_that
#' @return An array of bootstrap coefficient estimates of dimension p x L x num_boots
#'
#' @keywords internal
#' @noRd
run_boots <- function(data, X_base, Y_mat, weights = NULL, boot_type,
                      family = "gaussian", num_boots = 500,
                      seed = 2025, samp_method_by_stage = NULL,
                      parallel = FALSE, n_cores = NULL) {

  strata = psu = weight = NULL
  rm(list = c("strata", "psu", "weight"))

  if (!all(samp_method_by_stage %in% c("SRSWR", "SRSWOR", "PPSWR",
                                       "PPSWOR", "POISSON"))) {
    stop("Each element of `samp_method_by_stage` must be one of the following: \"SRSWR\", \"SRSWOR\", \"PPSWR\", \"PPSWOR\", or \"Poisson\"")
  }

  if (!(boot_type %in% c("BRR", "RWYB", "weighted", "unweighted"))) {
    stop("boot_type must be one of 'BRR', 'RWYB', 'weighted', 'unweighted'")
  }

  # ignore case for col names
  colnames(data) <- tolower(colnames(data))
  # add row ID for sampling
  data <- data %>% dplyr::mutate(row_id = dplyr::row_number())
  # handle if family input as character
  if (is.character(family)) family <- get(family, mode = "function", envir = parent.frame())()

  # ----- coefficient function to get coefficients given X, Y, w -----
  coef_fun <- function(X_tmp, Y_tmp, w_tmp = NULL) {
    if (family$family == "gaussian") {
      lm_wls_multi(X_tmp, Y_tmp, w = w_tmp)
    } else {
      glm_batch_multiY(X = X_tmp, Y = Y_tmp, w = w_tmp, family = family,
                       return_se = FALSE, add_intercept = FALSE)$coef
    }
  }

  set.seed(seed)
  # initalize array to store matrices
  betaTilde_boot <- NULL
  n = nrow(data)

  # parallel options
  if (parallel) {
    if (is.null(n_cores)) {
      n_cores <- parallelly::availableCores() - 1
      message("using ", n_cores, " cores for parallelization")
    }
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = n_cores)
    lapply_fn <- function(X, FUN) {
      future.apply::future_lapply(X, FUN, future.seed = TRUE)
    }
  } else {
    lapply_fn <- lapply
  }

  # ----- generate bootstrap indices or weights -----
  if (boot_type %in% c("weighted", "unweighted")) {
    # Determine the weight vector
    if (is.null(weights) & boot_type == "weighted") {
      stop("Weighted bootstrap requires a weight vector")
    } else if (is.null(weights) & boot_type == "unweighted") {
      weights <- rep(1, n)
    }
    assertthat::assert_that(is.numeric(weights), length(weights) == n, msg = "weights must be numeric with length `n`")
    assertthat::assert_that(all(weights > 0), sum(is.na(weights)) == 0, msg = "weights must be positive and non-missing")

    # create probability vector for sampling
    probs <- if (boot_type == "weighted") weights / sum(weights) else NULL

    # generate bootstrap indices
    boot_indices <- replicate(
      num_boots,
      sample(seq_len(nrow(data)), size = nrow(data), replace = TRUE, prob = probs),
      simplify = FALSE
    )

    coefs <- lapply_fn(boot_indices, function(idx) {
      X_tmp <- X_base[idx, , drop = FALSE]
      Y_tmp <- Y_mat[idx, , drop = FALSE]
      w_tmp <- if (!is.null(probs)) weights[idx] else NULL   # use weights if weighted otherwise unweighted
      coef_fun(X_tmp, Y_tmp, w_tmp)
    })

    betaTilde_boot <- simplify2array(coefs)
  }

  else if (boot_type == "BRR") {
    # Validate BRR setup
    if (!all(c("strata", "psu", "weight") %in% colnames(data))) stop("BRR requires strata, psu, weight")
    strata_counts <- data %>% dplyr::group_by(strata) %>%
      dplyr::summarize(n_psu = length(unique(psu)), .groups = "drop")
    if (any(strata_counts$n_psu != 2)) warning("Each strata should have exactly 2 PSUs for BRR")

    indices <- replicate(num_boots, sample_data(data), simplify = FALSE)

    coefs <- lapply_fn(indices, function(idx_df) {
      dat_tmp <- data %>% dplyr::inner_join(idx_df, by = c("psu", "strata")) %>%
        dplyr::mutate(weight = weight * 2)
      X_tmp <- X_base[dat_tmp$row_id, , drop = FALSE]
      Y_tmp <- Y_mat[dat_tmp$row_id, , drop = FALSE]
      coef_fun(X_tmp, Y_tmp, dat_tmp$weight)
    })
    betaTilde_boot <- simplify2array(coefs)
  }

  else if (boot_type == "RWYB") {
    if (!all(c("strata", "psu", "weight", "p_stage1", "p_stage2") %in% colnames(data))) {
      stop("RWYB bootstrap requires strata, psu, weight, p_stage1, p_stage2")
    }
    if (is.null(samp_method_by_stage)) stop("Specify samp_method_by_stage")
    if (!("id" %in% colnames(data))) {
      data <- data %>% dplyr::mutate(id = dplyr::row_number())
      message("Creating 'id' column for data")
    }

    svy_design <- survey::svydesign(ids = ~ psu + id, strata = ~ strata,
                                    weights = ~ weight, data = data, nest = TRUE)
    rwyb_wts <- svrep::make_rwyb_bootstrap_weights(
      num_replicates = num_boots,
      samp_unit_ids = svy_design$cluster,
      strata_ids = svy_design$strata,
      samp_unit_sel_probs = matrix(c(data$p_stage1, data$p_stage2), ncol = 2),
      samp_method_by_stage = samp_method_by_stage,
      allow_final_stage_singletons = TRUE,
      output = "weights"
    )

    coefs <- lapply_fn(seq_len(num_boots), function(r) {
      coef_fun(X_base, Y_mat, rwyb_wts[, r])
    })
    betaTilde_boot <- simplify2array(coefs)
  }

  betaTilde_boot
}


