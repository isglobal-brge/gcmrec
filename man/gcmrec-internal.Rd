% --- gcmrec-internal.Rd ---

\name{gcmrec-internal}
\alias{gcmrec-internal}
\alias{formatData}
\alias{formatData.effage}
\alias{formatData.i}
\alias{generlmi}
\alias{List.to.Dataframe}

\title{Internal gcmrec functions}
\description{Internal gcmrec functions}
\usage{

formatData(id, time, event, covariates, parameffage, cancer) 
formatData.effage(id, time, status, covariates, effageData)
formatData.i(id, time, event, covariates, parameffage, cancer = NULL)
generlmi(perrep)
List.to.Dataframe(data) 
}

\details{These are not to be called by the user}
\keyword{survival}
\keyword{internal}