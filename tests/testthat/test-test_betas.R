test_that("get_betatilde returns a matrix with expected dimensions", {
  # toy data
  set.seed(123)
  n <- 20
  L <- 50
  dat <- dplyr::tibble(
    Y = matrix(rnorm(n * L), nrow = n, ncol = L),
    X1 = rnorm(n),
    X2 = rnorm(n),
    weights = runif(n, 0.5, 1.5)
  )

  X_mat = dat %>% dplyr::select(starts_with("X")) %>% as.matrix()
  # run helper
  beta_tilde <- get_betatilde(
    data = dat,
    family = "gaussian",
    Y_mat = dat$Y,
    X_mat = X_mat,
    w = dat$weights
  )

  expect_true(is.matrix(beta_tilde))
  expect_equal(ncol(beta_tilde), L)
  expect_equal(nrow(beta_tilde), 2) # two predictors (X1, X2)
})


test_that("get_betahat smooths correctly", {
  # simulate raw coefficients
  L <- 40
  beta_tilde <- matrix(rnorm(2 * L), nrow = 2, ncol = L)

  # run smoother
  beta_hat <- get_betahat(beta_tilde, nknots_min = 5)

  expect_true(is.matrix(beta_hat))
  expect_equal(dim(beta_hat), dim(beta_tilde))

  # check smoothing doesn't introduce NAs
  expect_false(anyNA(beta_hat))
})
