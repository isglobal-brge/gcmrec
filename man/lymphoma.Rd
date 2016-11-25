\name{lymphoma}
\alias{lymphoma}
\non_function{}
\title{Indolent non-Hodgkin's lymphomas}
\usage{data(lymphoma)}
\description{
 This contains cancer relapses times after first treatment in patients diagnosed with low grade lymphoma
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{id}{identificator of each subject. Repeated for each recurrence}
    \item{time}{interocurrence or censoring time}
    \item{event}{censoring status. All event are 1 for each subject 
                 excepting last one that it is 0}
    \item{enum}{which lymphoma}
    \item{delay}{delay between first symptom and date of first treatment as a
                      continuous variable}
    \item{age}{age at diagnosis}
    \item{sex}{gender: 1:Males 2:Females}
    \item{distrib}{lesions involved at diagnosis (0=Single, 1=Localized, 2=More than
                       one nodal site, 3=Generalized}
    \item{effage}{response achieved after treatment upon relapses, coded as CR: Complete remission, 
                       PR: Partial remission or SD: stable disease or null response.}
  }
 
}
\source{

JR Gonzalez, E Peña, E Slate (2005). Modelling treatment effect after cancer relapses, with application to
recurrences in indolent non-Hodgkin's lymphomas. \emph{Stat Med}, submitted.

O. Servitje, F. Gallardo, T. Estrach, et al. (2002). Primary cutaneous marginal zone B-cell lymphoma: a clinical,
histopathological, immunophenotypic and molecular genetic study of 22 cases. \emph{Br J Dermatol},
147:1147-1158.
}

\keyword{datasets}
\keyword{survival}
