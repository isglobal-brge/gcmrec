#' Summary of a gcmrec model fit
#'
#' Computes hazard ratios (exponentiated covariate coefficients) with
#' confidence intervals for a fitted [gcmrec()] model. Unlike older
#' versions of the package, the summary is returned as an object (of class
#' `summary.gcmrec`) with its own print method, following the standard S3
#' idiom.
#'
#' @param object An object of class `gcmrec`.
#' @param level Confidence level for the intervals (default 0.95).
#' @param ... Ignored; kept for compatibility with the generic.
#'
#' @return An object of class `summary.gcmrec`, a list with components:
#'   \item{call}{the call that produced the fit}
#'   \item{conf.int}{matrix with the hazard ratio and its confidence
#'     limits for each covariate}
#'   \item{level}{the confidence level used}
#'   \item{loglik, n, nk}{fit summaries copied from the model object}
#'
#' @seealso [gcmrec()], [print.gcmrec()]
#'
#' @examples
#' data(readmission)
#' mod <- gcmrec(Survr(id, time, event) ~ as.factor(dukes),
#'               data = readmission, s = 3000)
#' summary(mod)
#' @export
summary.gcmrec <- function(object, level = 0.95, ...) {
  if (!inherits(object, "gcmrec")) {
    stop("Invalid data")
  }

  z <- abs(qnorm((1 - level) / 2))

  # Covariate coefficients: with rho = alpha^k the first coefficient is
  # alpha and is excluded from the hazard ratio table
  if (object$rho.type == 2) {
    co <- object$coef[-1]
    se <- sqrt(diag(object$var))[-1]
  } else {
    co <- object$coef
    se <- sqrt(diag(object$var))
  }
  se <- se[seq_along(co)]

  conf.int <- cbind(
    "exp(coef)" = exp(co),
    lower = exp(co - z * se),
    upper = exp(co + z * se)
  )
  rownames(conf.int) <- names(co)

  out <- list(
    call = object$call,
    conf.int = conf.int,
    level = level,
    loglik = object$loglik,
    n = object$n,
    nk = object$nk
  )
  class(out) <- "summary.gcmrec"
  out
}

#' Print a gcmrec summary
#'
#' @param x An object of class `summary.gcmrec`.
#' @param digits Number of decimal digits in the table.
#' @param lab Label of the hazard ratio column.
#' @param ... Ignored; kept for compatibility with the generic.
#'
#' @return Invisibly, `x`.
#' @seealso [summary.gcmrec()]
#' @export
print.summary.gcmrec <- function(x, digits = 2, lab = "hr", ...) {
  ci <- x$conf.int
  colnames(ci) <- c(lab, paste0(x$level * 100, "%"), "C.I.")
  formatted <- formatC(ifelse(ci > 999.99, Inf, ci), digits = digits,
                       width = 6, format = "f")

  labels <- rownames(ci)
  mx <- max(nchar(labels)) + 1
  cat(paste(rep(" ", mx), collapse = ""),
      paste("   ", colnames(ci)), "\n")
  for (i in seq_len(nrow(ci))) {
    cat(formatC(labels[i], width = mx), formatted[i, 1],
        "(", formatted[i, 2], "-", formatted[i, 3], ") \n")
  }
  invisible(x)
}
