#' Add a censored time equal to 0
#'
#' Adds a new row with a censored time equal to 0 for those subjects whose
#' end of follow-up coincides with the last occurrence (i.e. whose last
#' record is an event). [Survr()] requires every subject to have exactly
#' one censored record.
#'
#' @param datin Data frame containing id, time and event variables; other
#'   covariates are allowed and are copied to the added row.
#' @param id,time,event Column positions of the id, time and event
#'   variables in `datin`.
#'
#' @return A data frame with an added row (censored time equal to 0) for
#'   those subjects whose end of follow-up coincides with the last
#'   occurrence.
#'
#' @seealso [Survr()], [gcmrec()]
#'
#' @examples
#' # The second subject has no censored record: its follow-up ends at the
#' # last event
#' dat <- data.frame(
#'   id    = c(1, 1, 2, 2),
#'   time  = c(5, 3, 7, 4),
#'   event = c(1, 0, 1, 1)
#' )
#' addCenTime(dat)
#' @export
addCenTime <- function(datin, id = 1, time = 2, event = 3) {
  per_subject <- lapply(unique(datin[, id]), function(i) {
    subject <- datin[datin[, id] == i, , drop = FALSE]
    nrecs <- nrow(subject)
    if (subject[nrecs, event] == 1) {
      newrec <- subject[nrecs, ]
      newrec[, time] <- 0
      newrec[, event] <- 0
      subject <- rbind(subject, newrec)
    }
    subject
  })
  do.call(rbind, per_subject)
}
