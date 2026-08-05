#' Event chart of recurrent events in calendar time
#'
#' Draws one row per subject with a segment covering its follow-up, a
#' point at each recurrence and a marker at the end of follow-up. It is
#' the first look you should take at recurrent event data: it shows at a
#' glance how many subjects recur, how often, and whether the events
#' cluster early or late.
#'
#' @param data A data frame with `id`, `time` and `event` variables (and
#'   possibly covariates), or a legacy gcmrec list (see
#'   [as_gcmrec_data()]).
#' @param var Optional variable used to colour the subjects, of the same
#'   length as `nrow(data)`. Categorical variables colour by group;
#'   continuous ones use a gradient.
#' @param effageData Optional effective age information; recurrences with
#'   `perrepind == 0` (imperfect repair) are drawn as hollow triangles.
#' @param sortevents How to order the subjects along the vertical axis:
#'   `"followup"` (default) by length of follow-up, `"events"` by number
#'   of recurrences, or `"none"` to keep the order of the data. In earlier
#'   versions this argument was accepted but ignored.
#' @param width Size of the plotting symbols.
#' @param lines Draw the follow-up segment of each subject?
#' @param max.subjects Maximum number of subjects to draw. With more, a
#'   random sample of this size is taken (with a message), since the chart
#'   becomes unreadable beyond a few hundred rows. `Inf` draws all.
#' @param ... Ignored; kept for compatibility with earlier versions.
#'
#' @return A `ggplot` object.
#'
#' @seealso [mcf()] for the cumulative view of the same data.
#'
#' @examples
#' data(readmission)
#' graph.caltimes(readmission[readmission$id %in% 1:40, ])
#'
#' # Coloured by a covariate, ordered by number of events
#' graph.caltimes(readmission[readmission$id %in% 1:40, ],
#'                var = readmission$sex[readmission$id %in% 1:40],
#'                sortevents = "events")
#' @export
graph.caltimes <- function(data, var = NULL, effageData = NULL,
                           sortevents = c("followup", "events", "none"),
                           width = 2, lines = TRUE, max.subjects = 300,
                           ...) {
  data <- as_gcmrec_data(data)
  # Compatibility: in version 1.x sortevents was a logical
  if (is.logical(sortevents)) {
    sortevents <- if (isTRUE(sortevents)) "followup" else "none"
  }
  sortevents <- match.arg(sortevents)

  id.unique <- unique(data$id)
  if (length(id.unique) > max.subjects) {
    message("Showing a random sample of ", max.subjects, " of the ",
            length(id.unique), " subjects; raise 'max.subjects' to draw ",
            "more.")
    id.unique <- sort(sample(id.unique, max.subjects))
    keep <- data$id %in% id.unique
    if (!is.null(var)) var <- var[keep]
    data <- data[keep, , drop = FALSE]
  }

  # One record per subject: total follow-up and number of recurrences
  subjects <- data.frame(
    id = id.unique,
    tau = vapply(id.unique, function(i) sum(data$time[data$id == i]),
                 numeric(1)),
    nevents = vapply(id.unique,
                     function(i) sum(data$event[data$id == i] == 1),
                     numeric(1))
  )
  subjects$row <- switch(
    sortevents,
    followup = rank(subjects$tau, ties.method = "first"),
    events   = rank(subjects$nevents + subjects$tau / (max(subjects$tau) + 1),
                    ties.method = "first"),
    none     = seq_len(nrow(subjects))
  )

  # One record per recurrence, in calendar time
  events <- do.call(rbind, lapply(id.unique, function(i) {
    sel <- data$id == i & data$event == 1
    if (!any(sel)) return(NULL)
    data.frame(id = i, caltime = cumsum(data$time[sel]))
  }))

  # Imperfect repair: hollow triangles
  if (!is.null(effageData)) {
    perfect <- unlist(lapply(effageData$subject, function(x) x[["perrepind"]]),
                      use.names = FALSE)
    if (!is.null(events) && length(perfect) >= nrow(events)) {
      events$repair <- ifelse(perfect[seq_len(nrow(events))] == 0,
                              "imperfect", "perfect")
    }
  }
  if (!is.null(events) && is.null(events$repair)) {
    events$repair <- "perfect"
  }

  # Colour variable: one value per subject (the first of each). Without
  # `var` a constant group is used and the legend is hidden, so the
  # aesthetic mapping is the same in both cases.
  subjects$group <- if (is.null(var)) {
    factor("all")
  } else {
    var[match(subjects$id, data$id)]
  }
  continuous.group <- is.numeric(subjects$group) &&
    length(unique(subjects$group)) > 5
  if (!continuous.group) {
    subjects$group <- factor(subjects$group)
  }
  if (!is.null(events)) {
    events$row <- subjects$row[match(events$id, subjects$id)]
    events$group <- subjects$group[match(events$id, subjects$id)]
  }

  p <- ggplot2::ggplot()
  if (lines) {
    p <- p + ggplot2::geom_segment(
      data = subjects,
      ggplot2::aes(x = 0, xend = .data$tau, y = .data$row, yend = .data$row),
      colour = "grey85", linewidth = 0.35)
  }
  if (!is.null(events)) {
    p <- p + ggplot2::geom_point(
      data = events,
      ggplot2::aes(x = .data$caltime, y = .data$row, shape = .data$repair,
                   colour = .data$group),
      size = width * 0.75, stroke = 0.7)
  }
  p <- p +
    ggplot2::geom_point(
      data = subjects,
      ggplot2::aes(x = .data$tau, y = .data$row),
      shape = 4, size = width * 0.8, stroke = 0.8, colour = "grey45") +
    ggplot2::scale_shape_manual(
      values = c(perfect = 16, imperfect = 2), name = NULL,
      guide = if (!is.null(events) && length(unique(events$repair)) > 1)
        "legend" else "none") +
    ggplot2::labs(
      x = "Calendar time",
      y = switch(sortevents,
                 followup = "Subject (ordered by follow-up)",
                 events   = "Subject (ordered by number of events)",
                 none     = "Subject"),
      caption = "x marks the end of follow-up") +
    theme_gcmrec() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())

  p <- if (continuous.group) {
    p + ggplot2::scale_colour_viridis_c(name = NULL)
  } else {
    p + ggplot2::scale_colour_manual(
      values = if (is.null(var)) gcmrec_palette(1) else
        gcmrec_palette(nlevels(subjects$group)),
      name = NULL,
      guide = if (is.null(var)) "none" else "legend")
  }
  p
}
