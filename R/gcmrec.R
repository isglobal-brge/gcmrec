#' Fit the general class of models for recurrent event data
#'
#' Fits the general semiparametric model for recurrent events proposed by
#' Peña and Hollander (2004). This class of models incorporates an effective
#' age function which encodes the changes that occur after each event
#' occurrence (such as the impact of an intervention), allows for the
#' modelling of the impact of accumulating event occurrences on the unit,
#' admits a link function in which the effect of possibly time-dependent
#' covariates is incorporated, and allows the incorporation of unobservable
#' frailty components which induce dependence among the inter-event times
#' for each unit.
#'
#' Estimation with frailties is implemented using an
#' expectation-maximization (EM) algorithm. The marginal likelihood for
#' \eqn{\xi} is maximized either by Newton-Raphson or by Brent's method,
#' re-parameterizing \eqn{\xi^{*} = \log(\xi)} to avoid negative estimates.
#' Iteration terminates when successive values of \eqn{\xi / (1 + \xi)}
#' differ by no more than `tol`. Estimation under the frailty model can be
#' slow for large data sets.
#'
#' @param formula A formula object with a [Survr()] object as response (on
#'   the left of the `~` operator) and covariate terms on the right.
#'   Covariates are required.
#' @param data A data frame containing the variables named in the formula,
#'   including `id`, `time` and `event` (see [Survr()]). Alternatively, a
#'   list with elements `n` and `subject` in the legacy format described in
#'   [List.to.Dataframe()]; it is converted internally via
#'   [as_gcmrec_data()].
#' @param effageData Optional list with effective age information per
#'   subject (elements `intercepts`, `slopes`, `lastperrep`, `perrepind`,
#'   `effagebegin`, `effage`). If `NULL`, effective ages are generated under
#'   the perfect or minimal repair model according to `typeEffage`.
#' @param s A selected calendar time.
#' @param Frailty Logical. If `TRUE`, the model with gamma frailties is
#'   fitted via EM.
#' @param alphaSeed,betaSeed,xiSeed Initial values for \eqn{\alpha},
#'   \eqn{\beta} and \eqn{\xi}.
#' @param tol Tolerance of the maximization procedures.
#' @param maxit Maximum number of iterations of the maximization procedures.
#' @param rhoFunc Function \eqn{\rho(k; \alpha)} for the effect of
#'   accumulating event occurrences: `"alpha to k"` for
#'   \eqn{\rho(k; \alpha) = \alpha^k} (default) or `"Identity"` for
#'   \eqn{\rho(k; \alpha) = 1}. Partial matching is allowed.
#' @param typeEffage Effective age model: `"perfect"` (perfect repair,
#'   default) or `"minimal"` (minimal repair). Partial matching is allowed.
#' @param maxXi Maximization method for the marginal likelihood of
#'   \eqn{\xi}: `"Newton-Raphson"` (default) or `"Brent"`.
#' @param se Standard error estimation: `"Information matrix"` (default) or
#'   `"Jacknife"` (leave-one-out; can be time consuming).
#' @param cancer Optional character vector with the response achieved after
#'   each treatment for the cancer model of González et al. (2005), coded
#'   `"CR"`, `"PR"` or `"SD"`. See the [lymphoma] data set.
#' @param keep.data Keep the prepared data in the fitted object? It lets
#'   [anova()] refit the model without the original data, at the cost of a
#'   larger object. Set to `FALSE` for very large data sets.
#'
#' @return An object of class `gcmrec` with components including `coef`
#'   (estimated \eqn{\alpha} and \eqn{\beta}), `var` (covariance matrix),
#'   `loglik`, `Xi` and `frailties` (frailty model only), `Lambda` and
#'   `Survival` (baseline estimates at the distinct effective ages
#'   `diseff`). Methods are provided for [print()], [summary()], [plot()],
#'   [lines()], [coef()], [vcov()] and [logLik()].
#'
#' @references
#' Peña, E. and Hollander, M. (2004). Models for recurrent events in
#' reliability and survival analysis. In *Mathematical Reliability: An
#' Expository Perspective* (Chapter 6, pp. 105--123). Kluwer Academic
#' Publishers.
#'
#' González, J.R., Peña, E. and Slate, E. (2005). Modelling intervention
#' effects after cancer relapses. *Statistics in Medicine*, 24(24),
#' 3959--3975.
#'
#' Gail, M., Santner, T. and Brown, C. (1980). An analysis of comparative
#' carcinogenesis experiments based on multiple times to tumor.
#' *Biometrics*, 36, 255--266.
#'
#' @seealso [Survr()], [addCenTime()], [graph.caltimes()]
#'
#' @examples
#' data(readmission)
#'
#' # Perfect repair model with rho(k; alpha) = alpha^k
#' mod <- gcmrec(Survr(id, time, event) ~ as.factor(dukes),
#'               data = readmission, s = 3000)
#' mod
#'
#' # Minimal repair model
#' mod.min <- gcmrec(Survr(id, time, event) ~ as.factor(dukes),
#'                   data = readmission, s = 3000, typeEffage = "min")
#'
#' # Legacy list format converted internally
#' data(hydraulic)
#' gcmrec(Survr(id, time, event) ~ covar.1 + covar.2,
#'        data = hydraulic, s = 4753)
#'
#' \donttest{
#' # Model with frailties (EM algorithm)
#' mod.fra <- gcmrec(Survr(id, time, event) ~ as.factor(dukes),
#'                   data = readmission, s = 3000, Frailty = TRUE)
#'
#' # Cancer model with effective age given by treatment response
#' data(lymphoma)
#' mod.can <- gcmrec(Survr(id, time, event) ~ as.factor(distrib),
#'                   data = lymphoma, s = 1000, cancer = lymphoma$effage)
#' }
#' @export
gcmrec <- function(formula, data, effageData = NULL, s, Frailty = FALSE,
                   alphaSeed, betaSeed, xiSeed, tol = 10^(-6), maxit = 100,
                   rhoFunc = "alpha to k", typeEffage = "perfect",
                   maxXi = "Newton-Raphson", se = "Information matrix",
                   cancer = NULL, keep.data = TRUE) {
  # --- Argument validation -----------------------------------------------
  rho.type <- charmatch(rhoFunc, c("Identity", "alpha to k"), nomatch = 0)
  if (rho.type == 0) {
    stop("estimator must be 'alpha to k' or 'Identity' ")
  }
  effage.type <- charmatch(typeEffage, c("perfect", "minimal"))
  if (effage.type == 0) {
    stop("typeEffage must be perfect or minimal")
  }
  if (effage.type == 2) {
    effage.type <- 0
  }
  maxXi.type <- charmatch(maxXi, c("Newton-Raphson", "Brent"))
  if (maxXi.type == 0) {
    stop("maxXi must be Newton-Raphson or Brent")
  }
  se.type <- charmatch(se, c("Information matrix", "Jacknife"))
  if (se.type == 0) {
    stop("se must be Information matrix or Jacknife")
  }

  call <- match.call()
  if ((mode(call[[2]]) == "call" && call[[2]][[1]] == as.name("Survr")) ||
      inherits(formula, "Survr")) {
    stop("formula.default(object): invalid formula")
  }

  # --- Model frame construction ------------------------------------------
  m <- match.call(expand.dots = FALSE)
  m$s <- m$alphaSeed <- m$betaSeed <- m$xiSeed <- m$tol <- m$maxit <-
    m$rhoFunc <- m$Frailty <- m$effageData <- m$typeEffage <- m$maxXi <-
    m$se <- m$cancer <- m$keep.data <- m$... <- NULL
  Terms <- terms(formula, "strata")
  ord <- attr(Terms, "order")
  if (length(ord) && any(ord != 1)) {
    stop("Interaction terms are not valid for this function")
  }
  m$formula <- Terms
  # Normalize the input format (data.frame, legacy list, ...)
  if (!missing(data)) {
    m$data <- as_gcmrec_data(data)
  }
  m[[1]] <- as.name("model.frame")
  m <- eval(m, sys.parent())

  Y <- model.extract(m, "response")
  if (!is.Survr(Y)) {
    stop("Response must be a survival recurrent object")
  }

  offset <- attr(Terms, "offset")
  tt <- length(offset)
  offset <- if (tt == 0) {
    rep(0, nrow(Y))
  } else if (tt == 1) {
    m[[offset]]
  } else {
    ff <- m[[offset[1]]]
    for (i in 2:tt) ff <- ff + m[[offset[i]]]
    ff
  }

  mt <- attr(m, "terms")
  X <- if (!is.empty.model(mt)) model.matrix(mt, m)
  if (ncol(X) == 1) {
    stop("model need some covariates")
  }
  cov <- X[, -1]

  # --- Per-subject data preparation ---------------------------------------
  dataOK <- if (is.null(effageData)) {
    formatData(Y[, 1], Y[, 2], Y[, 3], cov, effage.type, cancer)
  } else {
    formatData.effage(Y[, 1], Y[, 2], Y[, 3], cov, effageData)
  }

  nvar <- ifelse(!is.null(ncol(cov)), ncol(cov), 1)
  if (missing(alphaSeed)) alphaSeed <- 1
  if (missing(betaSeed)) betaSeed <- rep(0, nvar)
  if (missing(xiSeed)) xiSeed <- 1

  # Covariate matrix, nvar x nk, for the C++ core (with a single
  # covariate formatData returns a vector)
  covMat <- matrix(dataOK$covariate, nrow = nvar)
  diseff <- sort(unique(dataOK$effage))

  # --- Estimation ----------------------------------------------------------
  if (!Frailty) {
    ans <- .gcmrec_fit(s, as.integer(dataOK$k), dataOK$tau,
                       dataOK$caltimes, dataOK$gaptimes, dataOK$censored,
                       dataOK$intercepts, dataOK$slopes, dataOK$lastperrep,
                       dataOK$effagebegin, dataOK$effage, covMat,
                       as.double(offset), alphaSeed, betaSeed,
                       rep(1, dataOK$n), rho.type, tol, as.integer(maxit))
    if (se.type == 2) {
      # Seeds for the leave-one-out refits: the estimates from the full
      # fit (with the identity rho, alpha is not estimated)
      if (rho.type == 2) {
        aSeed <- ans$estim[1]
        bSeed <- ans$estim[-1]
      } else {
        aSeed <- 1
        bSeed <- ans$estim
      }
      seJack <- .gcmrec_jackknife(s, as.integer(dataOK$k), dataOK$tau,
                                  dataOK$caltimes, dataOK$gaptimes,
                                  dataOK$censored, dataOK$intercepts,
                                  dataOK$slopes, dataOK$lastperrep,
                                  dataOK$effagebegin, dataOK$effage, covMat,
                                  as.double(offset), aSeed, bSeed, rho.type,
                                  tol, as.integer(maxit), resolve_threads())
      JackEstCov <- jackknife_cov(seJack)
    }
    if (ans$search != 1) {
      warning("Algorithm did not converge in: ", maxit, " iterations")
    }
    fit <- list(loglik = ans$loglik)
    fit$coef <- as.vector(ans$estim)
    fit$var <- if (se.type == 1) ans$info else JackEstCov
    fit$n <- dataOK$n
    fit$nk <- sum(dataOK$k)
    fit$kiter <- ans$kiter
    fit$Xi <- NULL
    fit$frailties <- rep(1, dataOK$n)
    fit$search <- ans$search
  } else {
    ans <- .gcmrec_frailty_fit(s, as.integer(dataOK$k), dataOK$tau,
                               dataOK$caltimes, dataOK$gaptimes,
                               dataOK$censored, dataOK$intercepts,
                               dataOK$slopes, dataOK$lastperrep,
                               dataOK$effagebegin, dataOK$effage, covMat,
                               as.double(offset), diseff, alphaSeed,
                               betaSeed, xiSeed, rep(1, dataOK$n), rho.type,
                               tol, as.integer(maxit), maxXi.type)
    ans$estim <- as.vector(ans$estim)
    if (se.type == 2) {
      seJack <- .gcmrec_jackknife_frailty(s, as.integer(dataOK$k),
                                          dataOK$tau, dataOK$caltimes,
                                          dataOK$gaptimes, dataOK$censored,
                                          dataOK$intercepts, dataOK$slopes,
                                          dataOK$lastperrep,
                                          dataOK$effagebegin, dataOK$effage,
                                          covMat, as.double(offset), diseff,
                                          ans$estim[1],
                                          ans$estim[2:(nvar + 1)],
                                          ans$estim[nvar + 2], rho.type,
                                          tol, as.integer(maxit),
                                          maxXi.type, resolve_threads())
      JackEstCov <- jackknife_cov(seJack)
    }
    fit <- list(loglik = ans$loglik)
    fit$search <- ans$search
    if (fit$search != 1) {
      warning("Algorithm did not converge after: ", maxit, " iterations")
    }
    fit$coef <- if (rho.type == 1) {
      ans$estim[2:(nvar + 1)]
    } else {
      ans$estim[1:(nvar + 1)]
    }
    fit$var <- if (se.type == 2) {
      JackEstCov
    } else if (rho.type == 1) {
      # Without the jackknife there is no variance estimate in the
      # frailty model (historical behaviour)
      matrix(NA, nvar + 1, nvar + 1)
    } else {
      matrix(NA, nvar + 2, nvar + 2)
    }
    fit$n <- dataOK$n
    fit$nk <- sum(dataOK$k)
    fit$kiter <- ans$kiter
    fit$Xi <- ans$estim[nvar + 2]
    fit$frailties <- ans$estim[(nvar + 3):(fit$n + nvar + 2)]
  }

  fit$method <- "Newton-Raphson"
  fit$rho.type <- rho.type

  # --- Baseline functions --------------------------------------------------
  if (fit$search == 1) {
    # With the identity rho, alpha is not part of coef: alpha = 1 is
    # passed (irrelevant for rho = 1) together with the full beta. The
    # original Fortran version sliced coef as if alpha were always
    # present, which left beta incomplete in this case.
    if (rho.type == 2) {
      alphaHat <- fit$coef[1]
      betaHat <- fit$coef[-1]
    } else {
      alphaHat <- 1
      betaHat <- fit$coef
    }
    EstLambSurv <- .gcmrec_baseline(s, as.integer(dataOK$k), dataOK$tau,
                                    dataOK$caltimes, dataOK$gaptimes,
                                    dataOK$censored, dataOK$intercepts,
                                    dataOK$slopes, dataOK$lastperrep,
                                    dataOK$effagebegin, dataOK$effage,
                                    covMat, as.double(offset), diseff,
                                    alphaHat, betaHat, fit$frailties,
                                    rho.type)
    fit$Lambda <- as.vector(EstLambSurv$Lambda)
    fit$Survival <- as.vector(EstLambSurv$Survival)
    fit$varLambda <- as.vector(EstLambSurv$varLambda)
  }

  fit$se.type <- se.type
  fit$diseff <- diseff
  fit$terms <- Terms
  # Factor levels from the fit, so that predict() can rebuild the same
  # design matrix from new data
  fit$xlevels <- .getXlevels(mt, m)
  # Prepared data and maximization options: they allow the model to be
  # refitted (for instance in anova()) without re-evaluating the call or
  # depending on the environment where it was fitted
  if (keep.data) {
    fit$gcmrec.data <- list(dataOK = dataOK, cov = covMat, offset = offset,
                            s = s, tol = tol, maxit = maxit,
                            alphaSeed = alphaSeed, betaSeed = betaSeed,
                            xiSeed = xiSeed, maxXi.type = maxXi.type,
                            se.type = se.type)
  }
  fit$call <- call
  if (rho.type == 2) {
    names(fit$coef) <- c("alpha", colnames(X)[-1])
  } else {
    names(fit$coef) <- colnames(X)[-1]
  }
  # The covariance matrix covers the coefficients and, in the frailty
  # model estimated by jackknife, also xi
  var.names <- c(names(fit$coef), "xi")[seq_len(nrow(fit$var))]
  if (length(var.names) == nrow(fit$var)) {
    dimnames(fit$var) <- list(var.names, var.names)
  }

  class(fit) <- "gcmrec"
  fit
}

# Jackknife covariance matrix from the matrix of leave-one-out
# estimates (one row per left-out subject)
jackknife_cov <- function(estimates) {
  n <- nrow(estimates)
  centered <- sweep(estimates, 2, colMeans(estimates))
  (n - 1) / n * crossprod(centered)
}
