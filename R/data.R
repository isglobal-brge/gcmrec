#' Rehospitalization of colorectal cancer patients
#'
#' Rehospitalization times after surgery in patients diagnosed with
#' colorectal cancer.
#'
#' @format A data frame with 861 rows and 10 columns:
#' \describe{
#'   \item{id}{identifier of each subject, repeated for each recurrence}
#'   \item{enum}{which readmission}
#'   \item{t.start}{start of interval (0 or previous recurrence time)}
#'   \item{t.stop}{recurrence or censoring time}
#'   \item{time}{inter-occurrence or censoring time}
#'   \item{event}{censoring status: 1 for each readmission, 0 for the last
#'     (censored) record of each subject}
#'   \item{chemo}{did the patient receive chemotherapy? 1: No; 2: Yes}
#'   \item{sex}{gender: 1: Males, 2: Females}
#'   \item{dukes}{Dukes' tumoral stage: 1: A-B; 2: C; 3: D}
#'   \item{charlson}{Charlson comorbidity index (time-dependent covariate):
#'     0: index 0; 1: index 1-2; 3: index >= 3}
#' }
#'
#' @source González, J.R., Fernández, E., Moreno, V. et al. (2005). Gender
#'   differences in hospital readmission among colorectal cancer patients.
#'   *Journal of Epidemiology and Community Health*, 59(6), 506--511.
#'
#' @examples
#' data(readmission)
#' head(readmission)
"readmission"

#' Indolent non-Hodgkin's lymphomas
#'
#' Cancer relapse times after first treatment in patients diagnosed with
#' low grade lymphoma.
#'
#' @format A data frame with 112 rows and 9 columns:
#' \describe{
#'   \item{id}{identifier of each subject, repeated for each recurrence}
#'   \item{time}{inter-occurrence or censoring time}
#'   \item{event}{censoring status: 1 for each relapse, 0 for the last
#'     (censored) record of each subject}
#'   \item{enum}{which relapse}
#'   \item{delay}{delay between first symptom and date of first treatment}
#'   \item{age}{age at diagnosis}
#'   \item{sex}{gender: 1: Males, 2: Females}
#'   \item{distrib}{lesions involved at diagnosis: 0: Single, 1: Localized,
#'     2: More than one nodal site, 3: Generalized}
#'   \item{effage}{response achieved after the treatment at each relapse:
#'     `"CR"` complete remission, `"PR"` partial remission, `"SD"` stable
#'     disease or null response. Used as the `cancer` argument of
#'     [gcmrec()]}
#' }
#'
#' @source González, J.R., Peña, E. and Slate, E. (2005). Modelling
#'   intervention effects after cancer relapses. *Statistics in Medicine*,
#'   24(24), 3959--3975.
#'
#'   Servitje, O., Gallardo, F., Estrach, T. et al. (2002). Primary
#'   cutaneous marginal zone B-cell lymphoma. *British Journal of
#'   Dermatology*, 147, 1147--1158.
#'
#' @examples
#' data(lymphoma)
#' head(lymphoma)
"lymphoma"

#' Hydraulic load-haul-dump (LHD) subsystems
#'
#' Hydraulic load-haul-dump (LHD) subsystems used in moving ore and rock in
#' underground mines in Sweden. The data set provides the calendar times
#' (in hours), excluding repair or down times, of the successive failures
#' of 6 such systems during the two-year development phase.
#'
#' @format A list in the legacy gcmrec format, with elements `n` (number of
#'   systems) and `subject` (a list with `k`, `tau`, `caltimes`,
#'   `gaptimes`, effective age components and covariates for each system).
#'   Use [List.to.Dataframe()] or pass it directly to [gcmrec()], which
#'   converts it via [as_gcmrec_data()].
#'
#' @source Kumar, D. and Klefsjö, B. (1992). Reliability analysis of
#'   hydraulic systems of LHD machines using the power law process model.
#'   *Reliability Engineering and System Safety*, 35, 217--224.
#'
#' @examples
#' data(hydraulic)
#' head(List.to.Dataframe(hydraulic))
"hydraulic"

#' Simulated data set generated under the minimal repair model
#'
#' Recurrence times simulated under the minimal repair model with
#' probability of perfect repair equal to 0.6, including the generated
#' effective age information. Useful as an example of the `effageData`
#' argument of [gcmrec()].
#'
#' @format A list in the legacy gcmrec format (elements `n`, `paramlist`
#'   and `subject`), where each subject includes the effective age
#'   components `intercepts`, `slopes`, `lastperrep`, `perrepind`,
#'   `effagebegin` and `effage`.
#'
#' @examples
#' data(GeneratedData)
#' temp <- List.to.Dataframe(GeneratedData)
#' head(temp)
"GeneratedData"
