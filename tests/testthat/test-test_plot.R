test_that("plot works", {
  ## create fake plot object
  tidy_df <- tibble::tibble(
    l = seq(0, 1, length.out = 10),
    beta_hat = rnorm(10),
    lower_pw = rnorm(10, -1),
    upper_pw = rnorm(10, 1),
    lower_joint = rnorm(10, -1.5),
    upper_joint = rnorm(10, 1.5),
    var_name = rep(c("X1", "X2"), each = 5)
  )

  fit <- list(betaHat = NULL, boots = NULL, cis = NULL, tidy_df = tidy_df)
  class(fit) <- "svyfui"

  # --- test with joint CI ---
  p1 <- plot(fit, include_joint_ci = TRUE)
  expect_s3_class(p1, "ggplot")

  # --- test without joint CI ---
  p2 <- plot(fit, include_joint_ci = FALSE)
  expect_s3_class(p2, "ggplot")

  # --- test that both calls produce different legends ---
  expect_false(identical(ggplot2::ggplot_build(p1)$plot$scales,
                         ggplot2::ggplot_build(p2)$plot$scales))
})
