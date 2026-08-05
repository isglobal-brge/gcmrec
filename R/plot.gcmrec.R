#' Plot the baseline functions of a gcmrec model
#'
#' Plots the estimated baseline survivor function or the baseline
#' cumulative hazard of a fitted [gcmrec()] model as a step function over
#' the effective age, with a pointwise confidence band.
#'
#' The band is built from the Aalen variance of the cumulative hazard,
#' transformed to the log scale for the hazard and to the log-log scale
#' for the survivor function, so that the limits stay positive and inside
#' \eqn{[0, 1]} respectively (the default transformation of
#' [survival::survfit()]).
#'
#' @param x An object of class `gcmrec`.
#' @param type.plot Curve to plot: `"survival"` (default) or `"hazard"`.
#'   Partial matching is allowed.
#' @param conf.int Draw the confidence band?
#' @param level Confidence level of the band.
#' @param ... Ignored; kept for compatibility with the generic.
#'
#' @return A `ggplot` object, which can be further styled with the usual
#'   `ggplot2` syntax.
#'
#' @seealso [plotBaseline()] to compare several models, [plotPredict()]
#'   for covariate-specific curves.
#'
#' @examples
#' data(readmission)
#' mod <- gcmrec(Survr(id, time, event) ~ as.factor(dukes),
#'               data = readmission, s = 3000)
#' plot(mod)
#' plot(mod, type.plot = "hazard", level = 0.99)
#'
#' # It is a ggplot, so it can be restyled
#' plot(mod) + ggplot2::labs(title = "Colorectal cancer readmissions")
#' @export
plot.gcmrec <- function(x, type.plot = "survival", conf.int = TRUE,
                        level = 0.95, ...) {
  df <- baseline_data(x, type.plot, conf.int = conf.int, level = level)
  ylab <- attr(df, "ylab")
  st <- step_frame(df)

  p <- ggplot2::ggplot(st, ggplot2::aes(x = .data$time))
  if (!is.null(st$lower)) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
      fill = gcmrec_palette(1), alpha = 0.15)
  }
  p +
    ggplot2::geom_line(ggplot2::aes(y = .data$estimate),
                       colour = gcmrec_palette(1), linewidth = 0.8) +
    ggplot2::labs(x = "Effective age", y = ylab,
                  caption = if (!is.null(st$lower))
                    paste0(level * 100, "% pointwise confidence band")) +
    theme_gcmrec()
}

#' Compare the baseline functions of several models
#'
#' Draws the baseline survivor or cumulative hazard functions of two or
#' more fitted models on the same axes, which is the natural way to see
#' what an assumption (for instance perfect versus minimal repair)
#' changes.
#'
#' @param models A named list of `gcmrec` objects. The names are used in
#'   the legend.
#' @param type.plot Curve to plot: `"survival"` (default) or `"hazard"`.
#' @param conf.int Draw confidence bands? Off by default, since
#'   overlapping bands are hard to read.
#' @param level Confidence level of the bands.
#'
#' @return A `ggplot` object.
#'
#' @seealso [plot.gcmrec()]
#'
#' @examples
#' data(readmission)
#' mod.per <- gcmrec(Survr(id, time, event) ~ as.factor(dukes),
#'                   data = readmission, s = 3000, typeEffage = "perfect")
#' mod.min <- gcmrec(Survr(id, time, event) ~ as.factor(dukes),
#'                   data = readmission, s = 3000, typeEffage = "minimal")
#' plotBaseline(list(perfect = mod.per, minimal = mod.min))
#' @export
plotBaseline <- function(models, type.plot = "survival", conf.int = FALSE,
                         level = 0.95) {
  if (inherits(models, "gcmrec")) {
    models <- list(models)
  }
  if (!all(vapply(models, inherits, logical(1), "gcmrec"))) {
    stop("'models' must be a list of gcmrec objects")
  }
  if (is.null(names(models))) {
    names(models) <- paste("model", seq_along(models))
  }

  curves <- lapply(names(models), function(nm) {
    df <- baseline_data(models[[nm]], type.plot, conf.int = conf.int,
                        level = level)
    st <- step_frame(df)
    st$model <- nm
    st
  })
  ylab <- attr(baseline_data(models[[1]], type.plot, conf.int = FALSE),
               "ylab")
  st <- do.call(rbind, curves)
  st$model <- factor(st$model, levels = names(models))

  p <- ggplot2::ggplot(st, ggplot2::aes(x = .data$time,
                                        colour = .data$model))
  if (conf.int && !is.null(st$lower)) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$lower, ymax = .data$upper,
                   fill = .data$model),
      alpha = 0.12, colour = NA) +
      ggplot2::scale_fill_manual(values = gcmrec_palette(length(models)),
                                 name = NULL)
  }
  p +
    ggplot2::geom_line(ggplot2::aes(y = .data$estimate), linewidth = 0.8) +
    ggplot2::scale_colour_manual(values = gcmrec_palette(length(models)),
                                 name = NULL) +
    ggplot2::labs(x = "Effective age", y = ylab) +
    theme_gcmrec()
}
