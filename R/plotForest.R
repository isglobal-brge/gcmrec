#' Forest plot of the hazard ratios
#'
#' Draws the exponentiated covariate coefficients of a fitted [gcmrec()]
#' model with their confidence intervals, on a logarithmic axis with a
#' reference line at 1. This is the figure that usually accompanies the
#' table returned by [summary.gcmrec()].
#'
#' @param x An object of class `gcmrec`.
#' @param level Confidence level of the intervals.
#' @param labels Optional named character vector to relabel the
#'   covariates: `c(sex = "Female vs male")`. Names are matched against
#'   the coefficient names.
#' @param sort Order the covariates by hazard ratio instead of keeping the
#'   order of the model?
#' @param title,subtitle Plot title and subtitle; `NULL` for none.
#'
#' @return A `ggplot` object.
#'
#' @seealso [summary.gcmrec()]
#'
#' @examples
#' data(readmission)
#' mod <- gcmrec(Survr(id, time, event) ~ as.factor(dukes) + sex,
#'               data = readmission, s = 3000)
#' plotForest(mod)
#'
#' # Readable labels for a report
#' plotForest(mod, labels = c("as.factor(dukes)2" = "Dukes C vs A-B",
#'                            "as.factor(dukes)3" = "Dukes D vs A-B",
#'                            "sex" = "Female vs male"))
#' @export
plotForest <- function(x, level = 0.95, labels = NULL, sort = FALSE,
                       title = "Hazard ratios", subtitle = NULL) {
  if (!inherits(x, "gcmrec")) {
    stop("'x' must be a gcmrec object")
  }
  ci <- summary(x, level = level)$conf.int
  df <- data.frame(
    term = rownames(ci),
    estimate = ci[, 1],
    lower = ci[, 2],
    upper = ci[, 3],
    stringsAsFactors = FALSE
  )

  if (!is.null(labels)) {
    hit <- match(df$term, names(labels))
    df$term[!is.na(hit)] <- unname(labels[hit[!is.na(hit)]])
  }
  df$term <- factor(df$term,
                    levels = if (sort) df$term[order(df$estimate)] else
                      rev(df$term))

  if (is.null(subtitle)) {
    subtitle <- paste0(level * 100, "% confidence intervals")
  }

  ggplot2::ggplot(df, ggplot2::aes(x = .data$estimate, y = .data$term)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed",
                        colour = "grey55") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$lower,
                                        xmax = .data$upper),
                           orientation = "y", width = 0.14,
                           colour = "grey35") +
    ggplot2::geom_point(size = 2.8, colour = gcmrec_palette(1)) +
    ggplot2::scale_x_continuous(trans = "log",
                                breaks = scales::breaks_log(n = 6)) +
    ggplot2::labs(x = "Hazard ratio (log scale)", y = NULL,
                  title = title, subtitle = subtitle) +
    theme_gcmrec() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_line(
                     linewidth = 0.3, colour = "grey88"))
}
