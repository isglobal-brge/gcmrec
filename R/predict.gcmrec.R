#' Predictions from a gcmrec model
#'
#' Computes the linear predictor, the relative risk or the whole survivor
#' and cumulative hazard curves for given covariate profiles.
#'
#' The model is multiplicative, so a covariate profile \eqn{x} scales the
#' baseline by \eqn{\psi(x) = \exp(x'\beta)}: the survivor function of the
#' profile is \eqn{S_0(t)^{\psi(x)}} and its cumulative hazard is
#' \eqn{\psi(x)\Lambda_0(t)}. Note that the curves are functions of the
#' **effective age**, not of calendar time, and that they describe the gap
#' to the next event.
#'
#' @param object An object of class `gcmrec`.
#' @param newdata Data frame of covariate profiles, one row per profile,
#'   with the covariates used in the model. If `NULL`, the profiles of the
#'   data used to fit are taken.
#' @param type Type of prediction: `"risk"` (relative risk
#'   \eqn{\exp(x'\beta)}, the default), `"lp"` (linear predictor
#'   \eqn{x'\beta}), `"survival"` or `"hazard"` (the whole curve for each
#'   profile).
#' @param ... Ignored; kept for compatibility with the generic.
#'
#' @return For `"risk"` and `"lp"`, a numeric vector with one value per
#'   profile. For `"survival"` and `"hazard"`, a data frame with columns
#'   `profile`, `time` and `estimate`, ready to be plotted (see
#'   [plotPredict()]).
#'
#' @seealso [plotPredict()], [plot.gcmrec()]
#'
#' @examples
#' data(readmission)
#' mod <- gcmrec(Survr(id, time, event) ~ as.factor(dukes) + sex,
#'               data = readmission, s = 3000)
#'
#' # Relative risk of two patient profiles
#' profiles <- data.frame(dukes = c(1, 3), sex = c(1, 1))
#' predict(mod, profiles, type = "risk")
#'
#' # Whole survivor curve of each profile
#' head(predict(mod, profiles, type = "survival"))
#' @export
predict.gcmrec <- function(object, newdata = NULL,
                           type = c("risk", "lp", "survival", "hazard"),
                           ...) {
  type <- match.arg(type)

  beta <- if (object$rho.type == 2) object$coef[-1] else object$coef
  X <- predict_model_matrix(object, newdata)
  lp <- as.vector(X %*% beta)
  names(lp) <- rownames(X)

  if (type == "lp") {
    return(lp)
  }
  if (type == "risk") {
    return(exp(lp))
  }

  if (is.null(object$Survival)) {
    stop("baseline functions not available (the fit did not converge)")
  }
  psi <- exp(lp)
  curves <- lapply(seq_along(psi), function(i) {
    estimate <- if (type == "survival") {
      object$Survival^psi[i]
    } else {
      psi[i] * object$Lambda
    }
    data.frame(profile = names(psi)[i], time = object$diseff,
               estimate = estimate, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, curves)
  out$profile <- factor(out$profile, levels = names(psi))
  attr(out, "type") <- type
  out
}

# Builds the design matrix for the profiles to be predicted, reusing
# the model terms so that transformations and contrasts match those of
# the fit.
predict_model_matrix <- function(object, newdata) {
  beta.names <- if (object$rho.type == 2) names(object$coef)[-1] else
    names(object$coef)

  if (is.null(newdata)) {
    stop("'newdata' is required: supply the covariate profiles to predict")
  }
  if (!is.data.frame(newdata)) {
    stop("'newdata' must be a data frame with one row per profile")
  }

  # Model terms without the Survr response. The levels stored in the
  # fit (xlevels) guarantee that factors generate the same columns even
  # if newdata does not contain all the levels.
  terms.cov <- delete.response(object$terms)
  X <- tryCatch(
    model.matrix(terms.cov,
                 model.frame(terms.cov, newdata, na.action = na.pass,
                             xlev = object$xlevels)),
    error = function(e) {
      stop("'newdata' does not contain the covariates of the model: ",
           conditionMessage(e), call. = FALSE)
    }
  )
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]

  if (!identical(colnames(X), beta.names)) {
    missing.cols <- setdiff(beta.names, colnames(X))
    if (length(missing.cols) > 0) {
      stop("'newdata' does not produce these model terms: ",
           paste(missing.cols, collapse = ", "),
           ".\n  Factors need the same levels as in the fitted data; ",
           "consider passing them as factors with matching levels.")
    }
    X <- X[, beta.names, drop = FALSE]
  }
  if (is.null(rownames(X))) {
    rownames(X) <- paste("profile", seq_len(nrow(X)))
  }
  X
}

#' Plot predicted curves for covariate profiles
#'
#' Draws the survivor or cumulative hazard curve of each covariate
#' profile, which is the way to show what the model implies for concrete
#' subjects (for instance, the three Dukes stages side by side).
#'
#' @param object An object of class `gcmrec`.
#' @param newdata Data frame of covariate profiles, one row per profile.
#'   Its row names, if any, are used as labels in the legend.
#' @param type `"survival"` (default) or `"hazard"`.
#' @param labels Optional character vector of labels for the profiles,
#'   one per row of `newdata`.
#'
#' @return A `ggplot` object.
#'
#' @seealso [predict.gcmrec()]
#'
#' @examples
#' data(readmission)
#' mod <- gcmrec(Survr(id, time, event) ~ as.factor(dukes) + sex,
#'               data = readmission, s = 3000)
#'
#' profiles <- data.frame(dukes = c(1, 2, 3), sex = 1)
#' plotPredict(mod, profiles,
#'             labels = c("Dukes A-B", "Dukes C", "Dukes D"))
#' @export
plotPredict <- function(object, newdata, type = c("survival", "hazard"),
                        labels = NULL) {
  type <- match.arg(type)
  if (!is.null(labels)) {
    if (length(labels) != nrow(newdata)) {
      stop("'labels' must have one value per row of 'newdata'")
    }
    rownames(newdata) <- labels
  }

  df <- predict(object, newdata, type = type)
  # Step curves, one per profile
  steps <- lapply(levels(df$profile), function(p) {
    sub <- df[df$profile == p, ]
    st <- dostep(sub$time, sub$estimate)
    data.frame(profile = p, time = st$x, estimate = st$y,
               stringsAsFactors = FALSE)
  })
  st <- do.call(rbind, steps)
  st$profile <- factor(st$profile, levels = levels(df$profile))

  ylab <- if (type == "survival") "Predicted survivor function" else
    "Predicted cumulative hazard"

  ggplot2::ggplot(st, ggplot2::aes(x = .data$time, y = .data$estimate,
                                   colour = .data$profile)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::scale_colour_manual(values = gcmrec_palette(nlevels(st$profile)),
                                 name = NULL) +
    ggplot2::labs(x = "Effective age", y = ylab) +
    theme_gcmrec()
}
