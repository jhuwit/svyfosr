test_that("glm_batch_multiY basic gaussian regression works", {
  set.seed(123)
  n <- 50; p <- 2; B <- 5
  X <- cbind(rnorm(n), rnorm(n))
  Y <- matrix(rnorm(n * B), n, B)
  X_mat <- stats::model.matrix(~ X)

  fit <- glm_batch_multiY(X, Y, family = gaussian(), add_intercept = FALSE)

  expect_type(fit, "list")
  expect_true(all(dim(fit$coef) == c(p, B)))
  expect_true(length(fit$converged) == B)
  expect_true(fit$iter <= 50)
})

test_that("glm_batch_multiY works with intercept and weights", {
  n <- 40; p <- 1; B <- 3
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * B), n, B)
  w <- runif(n, 0.5, 2)

  fit <- glm_batch_multiY(X, Y, w = w, family = gaussian(), add_intercept = TRUE)

  expect_equal(nrow(fit$coef), p + 1) # intercept added
  expect_true(all(fit$converged))
})

test_that("glm_batch_multiY handles binomial outcomes", {
  set.seed(42)
  n <- 60; p <- 1; B <- 2
  X <- matrix(rnorm(n * p), n, p)
  prob <- plogis(X) # true probs
  Y <- matrix(rbinom(n * B, size = 1, prob = rep(prob, B)), n, B)

  fit <- glm_batch_multiY(X, Y, family = binomial())

  expect_equal(dim(fit$coef), c(p, B))
  expect_true(all(fit$converged))
})

test_that("glm_batch_multiY errors with bad inputs", {
  X <- matrix(rnorm(20), 10, 2)
  Y <- matrix(rnorm(15), 5, 3)

  expect_error(glm_batch_multiY(X, Y, family = gaussian()), "same number of rows")

  Y2 <- matrix(rnorm(20), 10, 2)
  expect_error(glm_batch_multiY(X, Y2, w = 1:5, family = gaussian()), "length n")

  expect_error(glm_batch_multiY(X, Y2, offset = matrix(1, 5, 2), family = gaussian()))
})

test_that("glm_batch_multiY respects offsets", {
  n <- 30; p <- 1; B <- 2
  X <- matrix(0, n, p)  # no predictors
  Y <- matrix(rnorm(n * B), n, B)

  off <- matrix(1, n, B)
  fit_no_offset <- glm_batch_multiY(X, Y, family = gaussian())
  fit_with_offset <- glm_batch_multiY(X, Y, family = gaussian(), offset = off)

  expect_false(isTRUE(all.equal(fit_no_offset$coef, fit_with_offset$coef)))
})
