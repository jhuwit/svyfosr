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
    seed = 1012,
    L = L
  )
  expect_true(is.array(boots))
  expect_true(dim(boots) == c(ncol(X_base), L, 100))
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
    seed = 1012,
    L = L
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
    seed = 1012,
    L = L
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
    samp_stages= c("PPSWOR", "PPSWR"),
    family = gaussian(),
    num_boots = n_boot,
    seed = 1012,
    L = L
  )
  expect_true(is.array(boots))
  expect_true(dim(boots)[1] == ncol(X_base))
  expect_true(dim(boots)[2] == L)
  expect_true(dim(boots)[3] == 100)
})
