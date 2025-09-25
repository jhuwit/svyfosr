test_that("weighted bootstrap returns array with expected dimensions in parallel", {
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
    parallel = TRUE
  )
  expect_true(is.array(boots))
  expect_true(dim(boots)[1] == ncol(X_base))
  expect_true(dim(boots)[2] == L)
  expect_true(dim(boots)[3] == 100)
})
