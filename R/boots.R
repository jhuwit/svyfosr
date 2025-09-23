#' Bootstrapping
#'
#' @description
#' Computes bootstrap replicates for different bootstrap methods.
#' Implemented in \code{\link{svyfui}}
#' @param data A data frame with functional outcome columns, predictors, weights, etc.
#' @param boot_type Bootstrap method: "BRR", "Rao-Wu-Yue-Beaumont", "weighted", "unweighted"
#' @param family Error distribution family (e.g., "gaussian", "binomial", "poisson")
#' @param num_boots Number of bootstrap replicates (default 500)
#' @param seed Random seed for reproducibility (default 2025)
#' @param samp_stages Sampling method by stage for Rao-Wu-Yue-Beaumont bootstrap (e.g., c("PPSWR", "Poisson"))
#' @importFrom survey svydesign
#' @importFrom svrep make_rwyb_bootstrap_weights
#' @importFrom dplyr mutate group_by summarize row_number
#' @importFrom stats as.formula
#' @return An array of bootstrap coefficient estimates of dimension p x L x num_boots
#'
#' @keywords internal
#' @noRd
run_boots = function(data, X_base, Y_mat, boot_type, family = "gaussian",
                             num_boots = 500, seed = 2025, L,
                             samp_stages = NULL) {
    if(!(boot_type %in% c("BRR", "Rao-Wu-Yue-Beaumont", "weighted", "unweighted"))){
      stop("please specify boot_type as one of 'BRR', 'Rao-Wu-Yue-Beaumont', 'weighted', 'unweighted'")
    }
    colnames(data) = tolower(colnames(data)) # make colnames lowercase
    if (boot_type == "BRR") {
      # check strata and psu are in colnames
      if(!(all(c("strata", "psu", "weight")) %in% colnames(data))) stop("data must contain 'weight', 'strata' and 'psu' columns for BRR bootstrap")
      # check 2 psu per strata
      strata_counts = data %>%
        dplyr::group_by(strata) %>%
        dplyr::summarize(n_psu = length(unique(psu)))
      if(unique(strata_counts$n_psu) != 2) stop("each strata must contain at least 2 psu for BRR bootstrap")
    }
    if (boot_type == "weighted" & !("weight" %in% colnames(data))) stop("'weight' must be column in data for weighted bootstrap")
    if (boot_type == "Rao-Wu-Yue-Beaumont") {
      # check strata, psu, weight, p_stage1, p_stage2 are in colnames
      if(!(all(c("strata", "psu", "weight", "p_stage1", "p_stage2")) %in% colnames(data))) stop("data must contain 'strata', 'psu', 'weight', 'p_stage1', and 'p_stage2' columns for Rao-Wu-Yue-Beaumont bootstrap")
      if(is.null(samp_stages) | length(samp_stages) != 2) stop("please specify samp_stages as a vector of length 2 for Rao-Wu-Yue-Beaumont bootstrap")
      if(!("id" %in% colnames(data))){
        data = data %>% dplyr::mutate(id = dplyr::row_number())
        message("creating id column for data")
      }
  }
    argvals = 1:L
    data = data %>%
      dplyr::mutate(row_id = dplyr::row_number())

    if (is.character(family)) {
      family <- get(family, mode = "function", envir = parent.frame())()
    }
    if (boot_type == "BRR") {
      sample_data = function(df) {
        df %>%
          dplyr::group_by(strata) %>%
          dplyr::summarize(psu = sample(psu, size = 1))
      }
      set.seed(seed)
      indices = replicate(num_boots, sample_data(data), simplify = FALSE)
    } else if (boot_type == "weighted") {
      set.seed(seed)
      boot_indices <- replicate(
        num_boots,
        sample(
          seq_len(nrow(data)),
          size = nrow(data),
          replace = TRUE,
          prob = data$weight / sum(data$weight)
        ),
        simplify = FALSE
      )
    } else if (boot_type == "unweighted") {
      set.seed(seed)
      boot_indices <- replicate(num_boots,
                                sample(
                                  seq_len(nrow(data)),
                                  size = nrow(data),
                                  replace = TRUE
                                ),
                                simplify = FALSE)
    }
    if (boot_type == "BRR") {
      # Weighted least squares for Gaussian case
      coefs = lapply(1:num_boots, function(r) {
        dat_tmp = data %>%
          dplyr::inner_join(indices[[r]], by = c("psu", "strata")) %>%
          mutate(weight = weight * 2)

        X_tmp = X_base[dat_tmp$row_id, , drop = FALSE]
        Y_tmp = Y_mat[dat_tmp$row_id, , drop = FALSE]  # all L outcomes

        if (family$family == "gaussian") {
          # Fit all L responses at once
          coef_mat = lm_wls_multi(X_tmp, Y_tmp, dat_tmp$weight)
        } else {
          coef_mat = glm_batch_multiY(X = X_tmp, y_mat = Y_tmp, w = dat_tmp$weight, family = family, return_se = FALSE,
                                      add_intercept = FALSE)$coef
        }
        coef_mat  # (p × L)
      })

      # Convert list of matrices → array (p × L × num_boots)
      betaTilde_boot = simplify2array(coefs)
    } else if (boot_type == "unweighted") {
      coefs = lapply(1:num_boots, function(r) {
        dat_tmp <- data[boot_indices[[r]], ]
        X_tmp = X_base[dat_tmp$row_id, , drop = FALSE]
        Y_tmp = Y_mat[dat_tmp$row_id, , drop = FALSE]  # all L outcomes

        if (family$family == "gaussian") {
          # Fit all L responses at once
          coef_mat = lm_wls_multi(X_tmp, Y_tmp)
        } else {
          coef_mat = glm_batch_multiY(X = X_tmp, y_mat = Y_tmp, family = family, return_se = FALSE,
                                      add_intercept = FALSE)$coef
        }
        coef_mat  # (p × L)
      })
      betaTilde_boot = simplify2array(coefs)
    } else if (boot_type == "weighted") {
      coefs = lapply(1:num_boots, function(r) {
        dat_tmp <- data[boot_indices[[r]], ]
        X_tmp = X_base[dat_tmp$row_id, , drop = FALSE]
        Y_tmp = Y_mat[dat_tmp$row_id, , drop = FALSE]  # all L outcomes
        if (family$family == "gaussian") {
          # Fit all L responses at once
          coef_mat = lm_wls_multi(X_tmp, Y_tmp, dat_tmp$weight)
        } else {
          coef_mat = glm_batch_multiY(X = X_tmp, y_mat = Y_tmp, w = dat_tmp$weight, family = family, return_se = FALSE,
                                      add_intercept = FALSE)$coef
        }
        coef_mat  # (p × L)
      })
      # Convert list of matrices → array (p × L × num_boots)
      betaTilde_boot = simplify2array(coefs)
    } else if (boot_type == "Rao-Wu-Yue-Beaumont") {
      svy_design <- survey::svydesign(
        ids = ~ psu + id,
        strata = ~ strata,
        weights = ~ weight,
        data = data,
        nest = TRUE
      )
      set.seed(seed)
      rwyb_wts = svrep::make_rwyb_bootstrap_weights(
        num_replicates = num_boots,
        samp_unit_ids = svy_design$cluster,
        strata_ids = svy_design$strata,
        samp_unit_sel_probs = matrix(
          c(data$p_stage1, data$p_stage2),
          byrow = FALSE,
          ncol = 2
        ),
        samp_method_by_stage = samp_stages,
        allow_final_stage_singletons = TRUE,
        output = "weights"
      )

      coefs = lapply(1:num_boots, function(r) {
        wts_temp = rwyb_wts[, r]

        if (family$family == "gaussian") {
          # Fit all L responses at once
          coef_mat = lm_wls_multi(X_base, Y_mat, wts_temp)
        } else {
          coef_mat = glm_batch_multiY(X = X_base, y_mat = Y_mat, w = wts_temp, family = family, return_se = FALSE,
                                      add_intercept = FALSE)$coef
        }
        coef_mat  # (p × L)
      })
      # Convert list of matrices → array (p × L × num_boots)
      betaTilde_boot = simplify2array(coefs)
    }
    return(betaTilde_boot)

}
