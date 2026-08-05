#' Add gcmrec baseline curves to a base graphics plot
#'
#' @description
#' `r lifecycle_badge_deprecated()`
#'
#' Adds the baseline curve of a fitted model to an open **base graphics**
#' plot. Since version 2.0 [plot.gcmrec()] returns a `ggplot` object, so
#' `lines()` can no longer add a curve to it. Use [plotBaseline()] with a
#' list of models instead, which draws all of them on the same axes with a
#' legend.
#'
#' @param x An object of class `gcmrec`.
#' @param type.plot Curve to add: `"survival"` (default) or `"hazard"`.
#' @param ... Further graphical parameters passed to [graphics::lines()].
#'
#' @return Invisibly, `x`.
#'
#' @seealso [plotBaseline()], the recommended replacement.
#'
#' @examples
#' data(readmission)
#' mod.per <- gcmrec(Survr(id, time, event) ~ as.factor(dukes),
#'                   data = readmission, s = 3000, typeEffage = "perfect")
#' mod.min <- gcmrec(Survr(id, time, event) ~ as.factor(dukes),
#'                   data = readmission, s = 3000, typeEffage = "minimal")
#'
#' # Recommended: both curves in one call
#' plotBaseline(list(perfect = mod.per, minimal = mod.min))
#' @export
lines.gcmrec <- function(x, type.plot = "survival", ...) {
  .Deprecated("plotBaseline",
              msg = paste("'lines.gcmrec' is deprecated: plot.gcmrec() now",
                          "returns a ggplot object.\n  Use",
                          "plotBaseline(list(a = fit1, b = fit2)) to",
                          "compare models."))
  df <- baseline_data(x, type.plot, conf.int = FALSE)
  st <- dostep(df$time, df$estimate)
  lines(st, ...)
  invisible(x)
}

# Lifecycle badge used in the documentation (avoids depending on the
# lifecycle package just for a badge)
lifecycle_badge_deprecated <- function() {
  "**\\[Deprecated\\]**"
}
