#' Plot method for svyfui objects
#'
#' @param x An object of class \code{svyfui}, usually returned by \code{\link{svyfui}}.
#' @param include_joint_ci Logical; if TRUE, includes joint confidence intervals in the plot (default TRUE).
#' @param joint_fill Color for joint confidence interval fill (default "#56B4E9FF").
#' @param pw_fill Color for pointwise confidence interval fill (default "#0072B2FF").
#' @param joint_alpha Transparency level for joint confidence interval fill (default 0.3).
#' @param pw_alpha Transparency level for pointwise confidence interval fill (default 0.3).
#' @param ... Additional arguments passed to plotting functions.
#'
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line labs theme_classic facet_wrap scale_fill_manual element_text theme theme_sub_axis theme_sub_legend
#' @importFrom assertthat assert_that
#' @return A \code{ggplot} object (recommended).
#' @export
#'
#' @examples
#' \dontrun{
#' fit <- svyfui(Y ~ X1 + X2, data = dat, weights = "weights")
#' plot(fit)
#' }
plot.svyfui <- function(x, include_joint_ci = TRUE,
                        joint_fill = "#56B4E9FF", pw_fill = "#0072B2FF",
                        joint_alpha = 0.3, pw_alpha = 0.3, ...) {

  beta_hat = l = lower_joint = lower_pw = upper_joint = upper_pw = NULL
  rm(list = c("beta_hat", "l", "lower_joint", "lower_pw", "upper_joint", "upper_pw"))


  assertthat::assert_that(inherits(x, "svyfui"), msg = "plot.svyfui() requires an object of class 'svyfui'.")
  plt_df = x$tidy_df

  if(include_joint_ci) {
    p = plt_df %>%
      ggplot2::ggplot(ggplot2::aes(x = l, y = beta_hat)) +
      ggplot2::facet_wrap(.~var_name, scales = "free_y") +
      ggplot2::geom_ribbon(ggplot2::aes(x = l, ymin = lower_pw, ymax = upper_pw, fill = "Pointwise"), alpha = pw_alpha) +
      ggplot2::geom_ribbon(ggplot2::aes(x = l, ymin = lower_joint, ymax = upper_joint, fill = "Joint"), alpha = joint_alpha) +
      ggplot2::geom_line(linewidth = 1.1)  +
      ggplot2::labs(x = "Functional Domain", y = "Coefficient Estimate") +
      ggplot2::theme_classic() +
      ggplot2::theme_sub_legend(position = "bottom",
                                title = ggplot2::element_text(size = 14),
                                text = ggplot2::element_text(size = 12)) +
      ggplot2::theme_sub_axis(title = ggplot2::element_text(size = 14),
                              text = ggplot2::element_text(size = 12)) +
      ggplot2::theme(strip.text = ggplot2::element_text(size = 12)) +
      ggplot2::scale_fill_manual(values = c("Joint" = joint_fill, "Pointwise" =  pw_fill), name = "Confidence Interval")

  } else {
    p = plt_df %>%
      ggplot2::ggplot(ggplot2::aes(x = l, y = beta_hat)) +
      ggplot2::facet_wrap(.~var_name, scales = "free_y") +
      ggplot2::geom_ribbon(ggplot2::aes(x = l, ymin = lower_pw, ymax = upper_pw, fill = "Pointwise Confidence Interval"), alpha = pw_alpha) +
      ggplot2::geom_line(linewidth = 1.1)  +
      ggplot2::labs(x = "Functional Domain", y = "Coefficient Estimate") +
      ggplot2::theme_classic() +
      ggplot2::theme_sub_legend(position = "bottom",
                                title = ggplot2::element_blank(),
                                text = ggplot2::element_text(size = 12)) +
      ggplot2::theme_sub_axis(title = ggplot2::element_text(size = 14),
                              text = ggplot2::element_text(size = 12)) +
      ggplot2::theme(strip.text = ggplot2::element_text(size = 12)) +
      ggplot2::scale_fill_manual(values = c("Pointwise Confidence Interval" =  pw_fill))

  }

  return(p)

}
