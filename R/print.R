#' Print method for wdsmatch objects
#'
#' @param x A \code{wdsmatch} object returned by \code{\link{wdsmatchATE}}
#'   or \code{\link{wdsmatchATT}}.
#' @param digits Number of significant digits.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns the input object \code{x}. Called for its
#'   side effect of printing a formatted summary to the console, including
#'   the point estimate, standard error, confidence interval, number of
#'   matches, and sample sizes.
#' @export
print.wdsmatch <- function(x, digits = 4, ...) {
  cat("\nWeighted Double Score Matching\n")
  cat("------------------------------\n")
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat(sprintf("Estimate:  %.*f\n", digits, x$estimate))
  if (!is.na(x$se)) {
    cat(sprintf("Std. Err:  %.*f\n", digits, x$se))
    cat(sprintf("95%% CI:    [%.*f, %.*f]\n",
                digits, x$ci[1], digits, x$ci[2]))
  }
  cat(sprintf("M:         %d\n", x$M))
  cat(sprintf("N:         %d (treated: %d, control: %d)\n",
              x$n, x$n.treated, x$n.control))
  invisible(x)
}

#' Summary method for wdsmatch objects
#'
#' @param object A \code{wdsmatch} object returned by \code{\link{wdsmatchATE}}
#'   or \code{\link{wdsmatchATT}}.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns the input object. Called for its side effect
#'   of printing a formatted summary to the console.
#' @export
summary.wdsmatch <- function(object, ...) {
  print.wdsmatch(object, ...)
}
