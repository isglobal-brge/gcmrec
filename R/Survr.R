#' Create a survival recurrent object
#'
#' Creates the response object used by [gcmrec()] for recurrent event data.
#' Every subject must have exactly one censored time (its last record, with
#' `event = 0`); use [addCenTime()] when the end of follow-up coincides
#' with the last occurrence.
#'
#' @param id Subject identifier, repeated for each recurrence.
#' @param time Inter-occurrence or censoring time.
#' @param event Event indicator: 1 for an event, 0 for the censoring time.
#'
#' @return An object of class `Survr`: a matrix with columns `id`, `time`
#'   and `event`.
#'
#' @seealso [gcmrec()], [addCenTime()]
#'
#' @examples
#' data(readmission)
#' with(readmission, head(Survr(id, time, event)))
#' @export
Survr <- function(id, time, event) {
  if (length(unique(id)) != length(event[event == 0])) {
    stop("Data doesn't match. Every subject must have a censored time")
  }
  if (length(unique(event)) > 2 || max(event) != 1 || min(event) != 0) {
    stop("event must be 0-1")
  }

  ans <- cbind(id, time, event)
  oldClass(ans) <- "Survr"
  invisible(ans)
}

#' Test for a survival recurrent object
#'
#' @param x An R object.
#' @return `TRUE` if `x` inherits from class `Survr`.
#' @seealso [Survr()]
#' @examples
#' is.Survr(Survr(id = c(1, 1), time = c(5, 3), event = c(1, 0)))
#' @export
is.Survr <- function(x) {
  inherits(x, "Survr")
}
