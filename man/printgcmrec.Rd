\name{print.gcmrec}
\alias{print.gcmrec}

\title{Print a Short Summary of parameter estimates of a general class of models
       for reccurrent event data}
\description{
  Prints a short summary of 'gcmrec' object 
}
\usage{
\method{print}{gcmrec}(x, digits = max(options()$digits - 4, 3), ...)
}
\arguments{
  \item{x}{the result of a call to the gcmrec function}
  \item{digits}{number of digits to print}
  \item{\dots}{ other unusued arguments }
}
\value{
 x, with the invisible flag set
}


\seealso{ \code{\link{summary.gcmrec}},  \code{\link{gcmrec}} }

\examples{

data(lymphoma)
mod<-gcmrec(Survr(id,time,event)~as.factor(distrib),data=lymphoma,s=1000)
print(mod)


}

\keyword{survival }
