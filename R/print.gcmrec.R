#' Print a gcmrec model fit
#'
#' Prints the coefficient table (with standard errors, z statistics and
#' p-values), the \eqn{\alpha} parameter of the \eqn{\rho} function, the
#' frailty parameter when present, and fit summaries.
#'
#' @param x An object of class `gcmrec`.
#' @param digits Number of significant digits to print.
#' @param ... Ignored; kept for compatibility with the generic.
#'
#' @return Invisibly, `x`.
#'
#' @seealso [gcmrec()], [summary.gcmrec()]
#' @export
print.gcmrec <- function(x, digits = max(options()$digits - 4, 3), ...) {
  if (!is.null(cl <- x$call)) {
    cat("Call:\n")
    dput(cl)
    cat("\n")
  }
  if (!is.null(x$fail)) {
    cat(" gmcrec failed.", x$fail, "\n")
    return(invisible(x))
  }
  savedig <- options(digits = digits)
  on.exit(options(savedig))

  label.rho <- if (x$rho.type == 1) "Identity" else "Alpha to k"

  # With the identity rho, alpha is fixed at 1 and has no standard error
  coef.full <- if (x$rho.type == 1) c(1, x$coef) else x$coef
  coef <- coef.full[-1]

  se <- sqrt(diag(x$var))
  if (x$rho.type == 1) {
    se <- c(NA, se)
  }
  se <- if (is.null(x$Xi)) {
    se[-1]
  } else {
    n.Xi <- nrow(x$var)
    se[-c(1, n.Xi)]
  }
  if (is.null(coef) || is.null(se)) {
    stop("Input is not valid")
  }

  se.label <- if (x$se.type == 2) "se(coef) Jacknife" else "se(coef)"
  tmp <- cbind(coef, exp(coef), se, coef / se,
               signif(1 - pchisq((coef / se)^2, 1), digits - 1))
  dimnames(tmp) <- list(names(coef),
                        c("coef", "exp(coef)", se.label, "z", "p"))
  cat("\n")
  prmatrix(tmp)

  cat("\n")
  cat("  General class model parameter estimates", "\n")
  cat("    rho function: ", label.rho, "\n")

  alpha <- coef.full[1]
  alpha.se <- if (x$rho.type == 2) sqrt(x$var[1, 1]) else "--"
  alpha.label <- if (x$se.type == 2) "alpha (s.e. Jacknife): " else
    "alpha (s.e.): "
  cat("      ", alpha.label, alpha, " (", alpha.se, ")", "\n", sep = "")

  if (!is.null(x$Xi)) {
    xi.se <- sqrt(x$var[n.Xi, n.Xi])
    cat("    Frailty parameter, Xi (s.e. Jacknife): ", x$Xi, " (", xi.se,
        ")", "\n")
  }
  cat(" \n")

  # With frailties, the likelihood is conditional on the frailties
  # estimated in the last EM iteration
  loglik.label <- if (is.null(x$Xi)) "log-likelihood=" else
    "Partial log-likelihood (given frailties)="
  iter.label <- if (is.null(x$Xi)) "Newton-Raphson" else "EM steps"
  cat(paste("  ", loglik.label, round(x$loglik, 2), sep = ""), "\n")
  cat("  n=", x$n, "\n")
  cat("  n times=", x$nk, "\n")
  cat("  number of iterations: ", x$kiter, " ", iter.label, "\n")
  cat("\n")

  invisible(x)
}
