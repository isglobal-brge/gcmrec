#' Likelihood ratio tests for gcmrec models
#'
#' Compares nested `gcmrec` models with a likelihood ratio test. The
#' typical use is to ask whether there is real event accumulation (is
#' \eqn{\alpha} different from 1?) or whether a covariate is worth
#' keeping.
#'
#' Given a single model, the test compares it against the same model with
#' \eqn{\rho \equiv 1} (that is, \eqn{\alpha = 1}: occurrences carry no
#' effect), refitting it internally. Given two or more models, they are
#' compared in sequence, each against the previous one, in increasing
#' order of complexity.
#'
#' @section What can be compared:
#' A likelihood ratio test is only meaningful for **nested models fitted
#' to the same data**, so the following are checked and refused:
#'
#' \describe{
#'   \item{Different data}{all the models must have the same number of
#'     subjects and records.}
#'   \item{Different effective age}{models must share the effective age
#'     scale, that is, the same `s`, `typeEffage`, `effageData` and
#'     `cancer`. A perfect repair fit and a minimal repair one are *not*
#'     nested and cannot be tested against each other.}
#'   \item{Non-nested covariates}{the covariates of each model must be a
#'     subset of those of the next one.}
#'   \item{Frailties}{the log-likelihood of a frailty fit is conditional
#'     on the estimated frailties and is not on the same scale as that of
#'     a fit without them, so the two cannot be compared. Judge the
#'     frailties from \eqn{\xi} and its jackknife standard error instead.}
#' }
#'
#' Comparing two frailty models with each other is allowed but only
#' approximate, since each is conditional on its own estimated frailties;
#' a warning is issued.
#'
#' @param object An object of class `gcmrec`.
#' @param ... Further `gcmrec` objects to compare, in increasing order of
#'   complexity.
#' @param test Ignored; kept for compatibility with the generic.
#'
#' @return A data frame of class `anova.gcmrec` with the log-likelihood,
#'   degrees of freedom, chi-square statistic and p-value of each
#'   comparison.
#'
#' @seealso [logLik.gcmrec()], [gcmrec()]
#'
#' @examples
#' data(readmission)
#' mod <- gcmrec(Survr(id, time, event) ~ as.factor(dukes),
#'               data = readmission, s = 3000)
#'
#' # Is there event accumulation? (alpha = 1 against alpha free)
#' anova(mod)
#'
#' # Is sex worth adding to the model?
#' mod2 <- gcmrec(Survr(id, time, event) ~ as.factor(dukes) + sex,
#'                data = readmission, s = 3000)
#' anova(mod, mod2)
#' @export
anova.gcmrec <- function(object, ..., test = "Chisq") {
  others <- list(...)
  if (!all(vapply(others, inherits, logical(1), "gcmrec"))) {
    stop("all the objects compared must be gcmrec fits")
  }

  if (length(others) == 0) {
    # Test of alpha = 1: refit the model with the identity rho
    if (object$rho.type == 1) {
      stop("the model already has rho = Identity (alpha = 1); ",
           "supply a second model to compare")
    }
    reduced <- update_rho_identity(object)
    models <- list(reduced, object)
    labels <- c("rho = Identity (alpha = 1)", "rho = alpha^k")
  } else {
    models <- c(list(object), others)
    labels <- c(deparse(substitute(object)),
                vapply(substitute(list(...))[-1], deparse, character(1)))
  }

  check_comparable(models)

  loglik <- vapply(models, function(m) m$loglik, numeric(1))
  npar <- vapply(models, model_npar, numeric(1))

  out <- data.frame(
    model = labels,
    npar = npar,
    logLik = loglik,
    stringsAsFactors = FALSE
  )
  out$Df <- c(NA, diff(npar))
  out$Chisq <- c(NA, 2 * diff(loglik))
  out$`Pr(>Chisq)` <- c(NA, pchisq(out$Chisq[-1], out$Df[-1],
                                   lower.tail = FALSE))

  class(out) <- c("anova.gcmrec", "data.frame")
  out
}

# Checks that the models can be compared with a likelihood ratio test:
# same data, same effective age scale, nested covariates, and
# homogeneous with respect to frailties.
check_comparable <- function(models) {
  stop_anova <- function(...) stop(..., call. = FALSE)

  # 1. Same data
  n <- vapply(models, function(m) m$n, numeric(1))
  nk <- vapply(models, function(m) m$nk, numeric(1))
  if (length(unique(n)) > 1 || length(unique(nk)) > 1) {
    stop_anova("the models are not fitted to the same data: they have ",
               "different numbers\n  of subjects (",
               paste(unique(n), collapse = ", "), ") or records (",
               paste(unique(nk), collapse = ", "), ")")
  }

  # 2. Same effective age scale (same s, typeEffage, effageData,
  #    cancer): if it differs, the models are not nested
  diseff <- models[[1]]$diseff
  same.scale <- vapply(models[-1], function(m) {
    isTRUE(all.equal(m$diseff, diseff))
  }, logical(1))
  if (!all(same.scale)) {
    stop_anova("the models do not share the same effective age scale ",
               "(different s,\n  typeEffage, effageData or cancer), so ",
               "they are not nested")
  }

  # 3. Homogeneous frailties
  has.frailty <- vapply(models, function(m) !is.null(m$Xi), logical(1))
  if (length(unique(has.frailty)) > 1) {
    stop_anova("cannot compare models with and without frailties: their ",
               "log-likelihoods\n  are not on the same scale. See ",
               "?anova.gcmrec.")
  }
  if (all(has.frailty)) {
    warning("each frailty fit is conditional on its own estimated ",
            "frailties, so this\n  likelihood ratio test is only ",
            "approximate", call. = FALSE)
  }

  # 4. Nested covariates and increasing complexity
  covariates <- lapply(models, model_covariates)
  for (i in seq_len(length(models) - 1)) {
    if (!all(covariates[[i]] %in% covariates[[i + 1]])) {
      extra <- setdiff(covariates[[i]], covariates[[i + 1]])
      stop_anova("the models are not nested: model ", i, " has terms ",
                 "absent from model ", i + 1, "\n  (",
                 paste(extra, collapse = ", "), ")")
    }
  }
  npar <- vapply(models, model_npar, numeric(1))
  if (any(diff(npar) <= 0)) {
    stop_anova("the models must be given in increasing order of ",
               "complexity, and each\n  must add at least one parameter")
  }

  invisible(TRUE)
}

# Covariate names of a fit (alpha is not a covariate)
model_covariates <- function(m) {
  if (m$rho.type == 2) names(m$coef)[-1] else names(m$coef)
}

# Number of free parameters of a fit (coefficients plus xi)
model_npar <- function(m) {
  length(m$coef) + as.integer(!is.null(m$Xi))
}

# Refits a model forcing rho = Identity (alpha = 1).
#
# The prepared data that the fit stores in `gcmrec.data` are reused, so
# the refit depends neither on the environment where gcmrec() was called
# nor on the original data still being available.
update_rho_identity <- function(object) {
  d <- object$gcmrec.data
  if (is.null(d)) {
    stop("the fit does not carry the data needed to refit it (it was ",
         "produced with\n  keep.data = FALSE). Fit the reduced model ",
         "yourself and compare both:",
         "\n    m0 <- update(fit, rhoFunc = 'Identity'); anova(m0, fit)",
         call. = FALSE)
  }

  ans <- .gcmrec_fit(d$s, as.integer(d$dataOK$k), d$dataOK$tau,
                     d$dataOK$caltimes, d$dataOK$gaptimes,
                     d$dataOK$censored, d$dataOK$intercepts,
                     d$dataOK$slopes, d$dataOK$lastperrep,
                     d$dataOK$effagebegin, d$dataOK$effage, d$cov,
                     as.double(d$offset), d$alphaSeed, d$betaSeed,
                     rep(1, d$dataOK$n), 1L, d$tol, as.integer(d$maxit))
  if (ans$search != 1) {
    warning("the reduced model (alpha = 1) did not converge in ", d$maxit,
            " iterations", call. = FALSE)
  }

  reduced <- object
  reduced$loglik <- ans$loglik
  reduced$coef <- as.vector(ans$estim)
  names(reduced$coef) <- model_covariates(object)
  reduced$var <- ans$info
  reduced$rho.type <- 1
  reduced$kiter <- ans$kiter
  reduced$search <- ans$search
  reduced
}

#' @param x An object of class `anova.gcmrec`.
#' @param digits Number of significant digits.
#' @param ... Ignored; kept for compatibility with the generic.
#' @return For `print`, invisibly `x`.
#' @rdname anova.gcmrec
#' @export
print.anova.gcmrec <- function(x, digits = 4, ...) {
  cat("Likelihood ratio test for gcmrec models\n\n")
  out <- as.data.frame(x)
  rownames(out) <- format(out$model)
  out$model <- NULL
  printCoefmat(as.matrix(out), digits = digits, has.Pvalue = TRUE,
               P.values = TRUE, na.print = "", cs.ind = NULL,
               tst.ind = which(names(out) == "Chisq"))
  invisible(x)
}
