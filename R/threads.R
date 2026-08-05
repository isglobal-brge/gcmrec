#' Number of threads used by gcmrec
#'
#' Gets or sets the maximum number of threads used by the parallel parts of
#' the package, currently the jackknife standard errors, which refit the
#' model once per subject.
#'
#' By default the package lets OpenMP decide, which means it honours the
#' `OMP_NUM_THREADS` environment variable. Set this option to cap the
#' number of threads regardless of that variable, for instance when
#' running inside a shared machine or a job scheduler.
#'
#' If the package was built without OpenMP support (which is the case with
#' the default compiler on macOS), the computation runs sequentially and
#' this setting has no effect.
#'
#' @param n Maximum number of threads, or `NULL` to restore the default
#'   (let OpenMP decide). If missing, the current setting is returned
#'   without changing it.
#'
#' @return The current setting, invisibly when it is being changed.
#'   `NULL` means "let OpenMP decide".
#'
#' @examples
#' gcmrecThreads()        # current setting
#'
#' old <- gcmrecThreads(2)
#' gcmrecThreads(old)     # restore
#' @export
gcmrecThreads <- function(n) {
  if (missing(n)) {
    return(getOption("gcmrec.threads"))
  }
  if (!is.null(n)) {
    n <- as.integer(n)
    if (length(n) != 1 || is.na(n) || n < 1) {
      stop("'n' must be a single positive integer, or NULL")
    }
  }
  old <- getOption("gcmrec.threads")
  options(gcmrec.threads = n)
  invisible(old)
}

# Value passed to the C++ core: 0 means "let OpenMP decide"
resolve_threads <- function() {
  n <- getOption("gcmrec.threads")
  if (is.null(n)) 0L else as.integer(n)
}
