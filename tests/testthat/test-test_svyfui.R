test_that("svyfui works", {
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


  test_res = svyfui(Y ~ X1 + X2, data = dat, weights = weight, family = gaussian(),
                    boot_type = "RWYB", samp_method_by_stage = c("PPSWOR", "PPSWR"),
                    num_boots = n_boot, seed = 1012, nknots_min = 5, conf_level_pw = 0.95, conf_level_joint = 0.95)
  expect_true(class(test_res) == "svyfui")
  expect_true(is.list(test_res))
  expect_true(length(test_res) == 4)
  expect_true(is.matrix(test_res$betaHat))
  expect_equal(dim(test_res$betaHat), c(3, L))
  ## add some more
})
