#' Convert the legacy gcmrec list format to a data frame
#'
#' Converts a data set in the legacy list format (elements `n` and
#' `subject`, where each subject holds `k`, `tau`, `gaptimes` and a
#' covariate matrix; see the [hydraulic] and [GeneratedData] data sets)
#' into a data frame with columns `id`, `time`, `event` and `covar.*`.
#'
#' @param data A list with elements `n` (number of subjects) and `subject`
#'   (a list with one entry per subject).
#'
#' @return A data frame with one row per recurrence plus the censoring
#'   time of each subject.
#'
#' @seealso [as_gcmrec_data()], [gcmrec()]
#'
#' @examples
#' data(hydraulic)
#' head(List.to.Dataframe(hydraulic))
#' @export
List.to.Dataframe <- function(data) {
  per_subject <- lapply(seq_len(data$n), function(i) {
    subject <- data$subject[[i]]
    time_cen <- subject$tau - sum(subject$gaptime)
    cbind(
      id = i,
      time = c(subject$gaptime[-1], time_cen),
      event = c(rep(1, subject$k - 1), 0),
      t(subject$cov)
    )
  })

  out <- data.frame(do.call(rbind, per_subject))
  nvar <- nrow(data$subject[[1]]$cov)
  names(out) <- c("id", "time", "event", paste0("covar.", seq_len(nvar)))
  out
}
