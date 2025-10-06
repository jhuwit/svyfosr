#' Turn family from character into family object
#'
#' @description Turns character into family object if needed, and checks that it's a family object
#' @param family A GLM family object (e.g., from \code{binomial()}, \code{poisson()}), can also be character
#'
#' @return A family object
#'
#' @keywords internal
#' @noRd
resolve_family <- function(family) {
  # If a string is provided, convert it to a family object
  if (is.character(family)) {
    fam_fun <- match.fun(family)      # safely get function by name
    family <- fam_fun()               # call to get family object
  }

  # Check that it's actually a family object
  if (!inherits(family, "family")) {
    stop("'family' must be a family object or a character name of a family")
  }

  return(family)
}

#' Check that variable name is bare or external vector
#'
#' @description Checks input is a bare variable name or numeric vector
#' @param x A bare variable name or numeric vector
#'
#' @return NULL, throws error if not bare variable name or numeric vector
#'
#' @keywords internal
#' @noRd
check_var <- function(x) {
  expr <- substitute(x)

  if (is.character(expr)) {
    stop("Quoted variable names are not allowed. Use bare variable names instead, e.g. x = var, not x = 'var'.")
  }

  if (!is.symbol(expr)) {
    stop("Input must be a bare variable name, `$` expression, or external vector")
  }
}
