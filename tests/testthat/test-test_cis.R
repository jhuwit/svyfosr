test_that("CIs return list with expected dimensions", {
  ## create boots array
  set.seed(123)
  n <- 20
  L <- 50
  p <- 2
  n_boot <- 100
  beta_hat = matrix(rnorm(p * L), nrow = p, ncol = L)
  beta_tilde_boot <- array(rnorm(n * L * n_boot), dim = c(p, L, n_boot))
  ci_res = get_cis(beta_tilde_boot,
                   beta_hat,
                   L = L)
  expect_true(is.list(ci_res))
  expect_true(is.array(ci_res[[1]]))
  expect_equal(dim(ci_res[[1]]), c(L, L, p))
  expect_true(is.vector(ci_res[[2]]))
  expect_true(length(ci_res[[2]]) == p)

  ci_res = get_cis(beta_tilde_boot,
                   beta_hat,
                   L = L,
                   smooth_for_ci = FALSE,
                   smooth_for_variance = FALSE)
  expect_true(is.list(ci_res))
  expect_true(is.array(ci_res[[1]]))
  expect_equal(dim(ci_res[[1]]), c(L, L, p))
  expect_true(is.vector(ci_res[[2]]))
  expect_true(length(ci_res[[2]]) == p)

})
