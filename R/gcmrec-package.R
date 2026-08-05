#' gcmrec: General Class of Models for Recurrent Event Data
#'
#' Parameter estimation for the general semiparametric model for recurrent
#' event data proposed by Peña and Hollander (2004). The model incorporates
#' an effective age function encoding the impact of interventions after each
#' event occurrence, the effect of accumulating event occurrences, a link
#' function for possibly time-dependent covariates, and optional gamma
#' frailties to induce dependence among inter-event times.
#'
#' The main fitting function is [gcmrec()]. Recurrent event responses are
#' built with [Survr()]. Example data sets: [readmission], [lymphoma],
#' [hydraulic] and [GeneratedData].
#'
#' @references
#' Peña, E. and Hollander, M. (2004). Models for recurrent events in
#' reliability and survival analysis. In *Mathematical Reliability: An
#' Expository Perspective* (Chapter 6, pp. 105--123). Kluwer Academic
#' Publishers.
#'
#' González, J.R., Peña, E. and Slate, E. (2005). Modelling intervention
#' effects after cancer relapses. *Statistics in Medicine*, 24(24), 3959--3975.
#'
#' @useDynLib gcmrec, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @import survival
#' @importFrom stats .getXlevels delete.response is.empty.model model.extract
#'   model.frame model.matrix na.pass pchisq predict printCoefmat qnorm rbinom
#'   terms
#' @importFrom graphics lines
#' @importFrom grDevices colorRampPalette
#' @importFrom rlang .data
#' @keywords internal
"_PACKAGE"
