#' Print an `svyfui` object
#'
#' @param x an object of class [svyfui]
#' @param ... additional arguments to pass to other [print] methods
#'
#' @returns The object, invisibly
#' @export
#'
#' @examples
#' \donttest{
#'    fit <- svyfui(Y ~ X, data = svyfosr::sample_df, weights = weight)
#'    print(fit)
#' }
print.svyfui = function(x, ...) {
  cat("Survey FUI object\n")
  cat("-----------------\n")
  n_boot = dim(x$boots)[3]
  cat("Number of bootstrap replicates:", n_boot, "\n")
  print(x$tidy_df)
  invisible(x)
}
