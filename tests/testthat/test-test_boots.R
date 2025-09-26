test_that("weighted bootstrap returns array with expected dimensions", {
  set.seed(123)
  n <- 20
  L <- 50
  n_boot <- 100
  dat <- dplyr::tibble(
    Y = matrix(rnorm(n * L), nrow = n, ncol = L),
    X1 = rnorm(n),
    X2 = rnorm(n),
    weights = runif(n, 0.5, 1.5)
  )

  X_base = stats::model.matrix(~ X1 + X2, data = dat)
  Y_mat = dat$Y
  boots <- run_boots(
    data = dat,
    X_base,
    Y_mat,
    weights = dat$weights,
    boot_type = "weighted",
    family = gaussian(),
    num_boots = n_boot,
    seed = 1012
  )
  expect_true(is.array(boots))
  expect_true(dim(boots)[1] == ncol(X_base))
  expect_true(dim(boots)[2] == L)
  expect_true(dim(boots)[3] == 100)
})

test_that("unweighted bootstrap returns array with expected dimensions", {
  set.seed(123)
  n <- 20
  L <- 50
  n_boot <- 100
  dat <- dplyr::tibble(
    Y = matrix(rnorm(n * L), nrow = n, ncol = L),
    X1 = rnorm(n),
    X2 = rnorm(n),
    weights = runif(n, 0.5, 1.5)
  )

  X_base = stats::model.matrix(~ X1 + X2, data = dat)
  Y_mat = dat$Y
  boots <- run_boots(
    data = dat,
    X_base,
    Y_mat,
    weights = NULL,
    boot_type = "unweighted",
    family = gaussian(),
    num_boots = n_boot,
    seed = 1012
  )
  expect_true(is.array(boots))
  expect_true(dim(boots)[1] == ncol(X_base))
  expect_true(dim(boots)[2] == L)
  expect_true(dim(boots)[3] == 100)
})

test_that("BRR bootstrap returns array with expected dimensions", {
  set.seed(123)
  n <- 200
  L <- 50
  n_boot <- 100
  dat <- dplyr::tibble(
    Y = matrix(rnorm(n * L), nrow = n, ncol = L),
    X1 = rnorm(n),
    X2 = rnorm(n),
    weight = runif(n, 0.5, 1.5),
    strata = rep(1:(n/20), each = 20),
    psu = rep(1:2, n / 2)
  )

  X_base = stats::model.matrix(~ X1 + X2, data = dat)
  Y_mat = dat$Y
  boots <- run_boots(
    data = dat,
    X_base,
    Y_mat,
    weights = weights,
    boot_type = "BRR",
    family = gaussian(),
    num_boots = n_boot,
    seed = 1012
  )
  expect_true(is.array(boots))
  expect_true(dim(boots)[1] == ncol(X_base))
  expect_true(dim(boots)[2] == L)
  expect_true(dim(boots)[3] == 100)
})

test_that("RWYB bootstrap returns array with expected dimensions", {
  set.seed(123)
  n <- 200
  L <- 50
  n_boot <- 100
  dat <- dplyr::tibble(
    Y = matrix(rnorm(n * L), nrow = n, ncol = L),
    X1 = rnorm(n),
    X2 = rnorm(n),
    weight = runif(n, 0.5, 1.5),
    strata = rep(1:(n/20), each = 20),
    psu = rep(1:2, n / 2),
    p_stage1 = runif(n, 0.01, 0.2),
    p_stage2 = runif(n, 0.5, 0.9)
  )

  X_base = stats::model.matrix(~ X1 + X2, data = dat)
  Y_mat = dat$Y
  boots <- run_boots(
    data = dat,
    X_base,
    Y_mat,
    weights = weights,
    boot_type = "RWYB",
    samp_method_by_stage= c("PPSWOR", "PPSWR"),
    family = gaussian(),
    num_boots = n_boot,
    seed = 1012
  )
  expect_true(is.array(boots))
  expect_true(dim(boots)[1] == ncol(X_base))
  expect_true(dim(boots)[2] == L)
  expect_true(dim(boots)[3] == 100)
})

test_that("run_boots errors with bad inputs", {
  set.seed(123)
  n <- 20
  L <- 50
  n_boot <- 10
  dat <- dplyr::tibble(
    Y = matrix(rnorm(n * L), nrow = n, ncol = L),
    X1 = rnorm(n),
    X2 = rnorm(n),
    weight = runif(n, 0.5, 1.5),
    psu = rep(1:2, n / 2),
    strata = rep(1:(n/4), each = 4)
  )

  X_base = stats::model.matrix(~ X1 + X2, data = dat)
  Y_mat = dat$Y

  expect_error(run_boots(
    data = dat,
    X_base,
    Y_mat,
    weights = dat$weight,
    boot_type = "invalid",
    family = gaussian(),
    num_boots = n_boot,
    seed = 1012
  ), "boot_type must be one of")

  expect_error(run_boots(
    data = dat,
    X_base,
    Y_mat,
    weights = dat$weight,
    boot_type = "RWYB",
    family = gaussian(),
    num_boots = n_boot,
    samp_method_by_stage = c("SSR", "PSR"),
    seed = 1012
  ), "must be one of the following")

  expect_error(run_boots(
    data = dat,
    X_base,
    Y_mat,
    weights = NULL,
    boot_type = "weighted",
    family = gaussian(),
    num_boots = n_boot,
    seed = 1012
  ), "requires a weight vector")

  expect_error(run_boots(
    data = dat %>% dplyr::select(-strata, -psu),
    X_base,
    Y_mat,
    weights = NULL,
    boot_type = "BRR",
    family = gaussian(),
    num_boots = n_boot,
    seed = 1012
  ), "requires strata, psu")

  expect_error(run_boots(
    data = dat,
    X_base,
    Y_mat,
    weights = dat$weight,
    boot_type = "RWYB",
    samp_method_by_stage = c("PPSWR", "PPSWR"),
    family = gaussian(),
    num_boots = n_boot,
    seed = 1012
  ), "requires strata, psu, weight, p_stage1")

  expect_error(run_boots(
    data = dat %>% dplyr::mutate(p_stage1 = 1, p_stage2 = 1),
    X_base,
    Y_mat,
    weights = dat$weight,
    boot_type = "RWYB",
    family = gaussian(),
    num_boots = n_boot,
    seed = 1012
  ), "Specify samp")

  dat_psu <- dplyr::tibble(
    Y = matrix(rnorm(n * L), nrow = n, ncol = L),
    X1 = rnorm(n),
    X2 = rnorm(n),
    weight = runif(n, 0.5, 1.5),
    psu = c(rep(1:3, 6), 1:2),
    strata = rep(1:(n/4), each = 4)
  )

  expect_warning(run_boots(
    data = dat_psu,
    X_base,
    Y_mat,
    weights = dat$weight,
    boot_type = "BRR",
    family = gaussian(),
    num_boots = n_boot,
    seed = 1012
  ), "exactly 2")

})

test_that("resolve_family works as expected", {
  fam1 <- resolve_family("gaussian")
  expect_true(inherits(fam1, "family"))
  expect_equal(fam1$family, "gaussian")

  fam2 <- resolve_family(stats::poisson())
  expect_true(inherits(fam2, "family"))
  expect_equal(fam2$family, "poisson")

  expect_error(resolve_family(123), "'family' must be a family object or a character name of a family")
})
