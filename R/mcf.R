#' Mean cumulative function of recurrent events
#'
#' Estimates the mean cumulative function (MCF), that is, the expected
#' number of events accumulated by a subject up to each time. It is the
#' Nelson-Aalen estimator adapted to recurrent events, and it answers the
#' question *"how many events does an average subject accumulate?"*
#' without assuming any model.
#'
#' The estimator adds, at each event time \eqn{t}, the number of events
#' observed at \eqn{t} divided by the number of subjects still under
#' follow-up, and the pointwise variance is accumulated in the same way.
#' Use it before modelling: comparing the MCF of two groups shows whether
#' their event burden really differs, and its shape suggests whether the
#' risk grows or falls with time.
#'
#' @param data A data frame with `id`, `time` and `event` variables, or a
#'   legacy gcmrec list (see [as_gcmrec_data()]).
#' @param group Optional grouping variable of length `nrow(data)`; one MCF
#'   is estimated per group.
#' @param level Confidence level of the pointwise interval.
#'
#' @return An object of class `mcf`: a data frame with columns `time`,
#'   `estimate`, `lower`, `upper`, `n.risk`, `n.event` and, if `group` was
#'   given, `group`. It has [plot()] and [print()] methods.
#'
#' @seealso [graph.caltimes()] for the individual view of the same data.
#'
#' @examples
#' data(readmission)
#'
#' # Overall event burden
#' m <- mcf(readmission)
#' m
#' plot(m)
#'
#' # Does it differ by tumour stage?
#' m.dukes <- mcf(readmission, group = readmission$dukes)
#' plot(m.dukes)
#' @export
mcf <- function(data, group = NULL, level = 0.95) {
  data <- as_gcmrec_data(data)
  if (is.null(data$id) || is.null(data$time) || is.null(data$event)) {
    stop("'data' must contain the variables id, time and event")
  }

  if (is.null(group)) {
    out <- mcf_one(data, level)
  } else {
    if (length(group) != nrow(data)) {
      stop("'group' must have one value per row of 'data'")
    }
    group <- factor(group)
    per_group <- lapply(levels(group), function(g) {
      res <- mcf_one(data[group == g, , drop = FALSE], level)
      res$group <- g
      res
    })
    out <- do.call(rbind, per_group)
    out$group <- factor(out$group, levels = levels(group))
  }

  class(out) <- c("mcf", "data.frame")
  attr(out, "level") <- level
  out
}

# MCF of a single group. Walks over the event times on the calendar
# scale and accumulates dN / n.risk, with the Nelson-Aalen variance.
mcf_one <- function(data, level) {
  id.unique <- unique(data$id)

  # Calendar time of each event and end of follow-up per subject
  per_subject <- lapply(id.unique, function(i) {
    sel <- data$id == i
    list(events = cumsum(data$time[sel & data$event == 1]),
         tau = sum(data$time[sel]))
  })
  events <- unlist(lapply(per_subject, `[[`, "events"), use.names = FALSE)
  tau <- vapply(per_subject, `[[`, numeric(1), "tau")

  times <- sort(unique(events))
  if (length(times) == 0) {
    stop("no events found in the data")
  }

  # Subjects still under follow-up and events observed at each time
  n.risk <- vapply(times, function(t) sum(tau >= t), numeric(1))
  n.event <- vapply(times, function(t) sum(events == t), numeric(1))

  increment <- ifelse(n.risk > 0, n.event / n.risk, 0)
  estimate <- cumsum(increment)
  var.mcf <- cumsum(ifelse(n.risk > 0, n.event / n.risk^2, 0))

  z <- abs(qnorm((1 - level) / 2))
  # Log-scale interval to keep the lower limit positive
  se.log <- sqrt(var.mcf) / pmax(estimate, 1e-12)

  data.frame(
    time = times,
    estimate = estimate,
    lower = pmax(estimate * exp(-z * se.log), 0),
    upper = estimate * exp(z * se.log),
    n.risk = n.risk,
    n.event = n.event
  )
}

#' @param x An object of class `mcf`.
#' @param ... Ignored; kept for compatibility with the generic.
#' @return For `print`, invisibly `x`.
#' @rdname mcf
#' @export
print.mcf <- function(x, ...) {
  cat("Mean cumulative function\n")
  if (!is.null(x$group)) {
    for (g in levels(x$group)) {
      sub <- x[x$group == g, ]
      cat(sprintf("  %-16s %3d event times, %.2f events per subject by t=%g\n",
                  g, nrow(sub), sub$estimate[nrow(sub)], max(sub$time)))
    }
  } else {
    cat(sprintf("  %d event times, %.2f events per subject by t=%g\n",
                nrow(x), x$estimate[nrow(x)], max(x$time)))
  }
  invisible(x)
}

#' @param conf.int Draw the confidence band?
#' @rdname mcf
#' @export
plot.mcf <- function(x, conf.int = TRUE, ...) {
  level <- attr(x, "level")
  has.group <- !is.null(x$group)

  # The MCF is a step function starting at (0, 0)
  curves <- if (has.group) split(x, x$group) else list(x)
  steps <- lapply(names(curves) %||% seq_along(curves), function(nm) {
    df <- curves[[nm]]
    st <- data.frame(time = c(0, rep(df$time, each = 2)),
                     estimate = c(0, 0, rep(df$estimate, each = 2))[
                       seq_len(2 * nrow(df) + 1)])
    st$lower <- c(0, rep(df$lower, each = 2))[seq_len(nrow(st))]
    st$upper <- c(0, rep(df$upper, each = 2))[seq_len(nrow(st))]
    if (has.group) st$group <- nm
    st
  })
  st <- do.call(rbind, steps)
  if (has.group) st$group <- factor(st$group, levels = levels(x$group))

  p <- ggplot2::ggplot(st, ggplot2::aes(x = .data$time))
  if (conf.int) {
    p <- p + if (has.group) {
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower,
                                        ymax = .data$upper,
                                        fill = .data$group),
                           alpha = 0.12, colour = NA)
    } else {
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower,
                                        ymax = .data$upper),
                           fill = gcmrec_palette(1), alpha = 0.15)
    }
  }
  p <- p + if (has.group) {
    ggplot2::geom_line(ggplot2::aes(y = .data$estimate,
                                    colour = .data$group),
                       linewidth = 0.8)
  } else {
    ggplot2::geom_line(ggplot2::aes(y = .data$estimate),
                       colour = gcmrec_palette(1), linewidth = 0.8)
  }

  if (has.group) {
    n.groups <- nlevels(st$group)
    p <- p +
      ggplot2::scale_colour_manual(values = gcmrec_palette(n.groups),
                                   name = NULL) +
      ggplot2::scale_fill_manual(values = gcmrec_palette(n.groups),
                                 name = NULL)
  }

  p +
    ggplot2::labs(x = "Calendar time",
                  y = "Expected number of events per subject",
                  caption = if (conf.int)
                    paste0(level * 100, "% pointwise confidence band")) +
    theme_gcmrec()
}

# Null-coalescing operator, to avoid depending on rlang
`%||%` <- function(a, b) if (is.null(a)) b else a
