#' Extract coefficients from a gcmrec model
#'
#' @param object An object of class `gcmrec`.
#' @param ... Ignored; kept for compatibility with the generic.
#'
#' @return Named vector of estimated parameters. With
#'   `rhoFunc = "alpha to k"` the first element is \eqn{\alpha}.
#' @seealso [gcmrec()]
#' @export
coef.gcmrec <- function(object, ...) {
  object$coef
}

#' Extract the covariance matrix from a gcmrec model
#'
#' @param object An object of class `gcmrec`.
#' @param ... Ignored; kept for compatibility with the generic.
#'
#' @return Covariance matrix of the estimated parameters (inverse of the
#'   information matrix, or the jackknife estimate if the model was fitted
#'   with `se = "Jacknife"`).
#' @seealso [gcmrec()]
#' @export
vcov.gcmrec <- function(object, ...) {
  object$var
}

#' Extract the log-likelihood from a gcmrec model
#'
#' @param object An object of class `gcmrec`.
#' @param ... Ignored; kept for compatibility with the generic.
#'
#' @return An object of class `logLik`. For frailty models it is the
#'   partial log-likelihood conditional on the estimated frailties.
#' @seealso [gcmrec()]
#' @export
logLik.gcmrec <- function(object, ...) {
  structure(object$loglik,
            df = length(object$coef) + as.integer(!is.null(object$Xi)),
            nobs = object$n,
            class = "logLik")
}
