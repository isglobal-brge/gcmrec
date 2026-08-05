#' Normalize input data for gcmrec
#'
#' Generic that converts the supported input formats of [gcmrec()] and
#' [graph.caltimes()] into the canonical `data.frame` layout with `id`,
#' `time` and `event` columns. New input formats can be supported by
#' adding methods to this generic.
#'
#' Currently supported inputs:
#' * `data.frame` (including tibbles and `data.table`s): returned as is.
#' * Legacy list format with elements `n` and `subject` (see
#'   [List.to.Dataframe()]): converted to a data frame. Any other list is
#'   returned unchanged so it can be used as a variable environment for
#'   the model formula, as in base R modelling functions.
#'
#' @param data The data in one of the supported formats.
#' @param ... Additional arguments passed to methods.
#'
#' @return A `data.frame` (or the original object when no conversion
#'   applies).
#'
#' @seealso [gcmrec()], [List.to.Dataframe()]
#'
#' @examples
#' data(hydraulic)
#' head(as_gcmrec_data(hydraulic))
#' @export
as_gcmrec_data <- function(data, ...) {
  UseMethod("as_gcmrec_data")
}

#' @rdname as_gcmrec_data
#' @export
as_gcmrec_data.data.frame <- function(data, ...) {
  data
}

#' @rdname as_gcmrec_data
#' @export
as_gcmrec_data.list <- function(data, ...) {
  if (!is.null(data$n) && !is.null(data$subject)) {
    List.to.Dataframe(data)
  } else {
    # Ordinary list of variables: model.frame() knows how to handle it
    data
  }
}

#' @rdname as_gcmrec_data
#' @export
as_gcmrec_data.default <- function(data, ...) {
  stop("'data' must be a data.frame or a legacy gcmrec list ",
       "(see ?List.to.Dataframe); got an object of class '",
       paste(class(data), collapse = "/"), "'")
}
