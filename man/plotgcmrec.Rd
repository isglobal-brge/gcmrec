\name{plot.gcmrec}
\alias{plot.gcmrec}
\alias{lines.gcmrec}


\title{Plot Method for an object of class `gcmrec'.}

\description{
  Plots estimated baseline survival and hazard functions from an object of class `gcmrec'. 
}

\usage{
\method{plot}{gcmrec}(x, type.plot = "surv", ...)
}

\arguments{
  \item{x}{ Object of class \code{gcmrec} (output from calling \code{gcmrec}
    function).}
  \item{type.plot}{ a character string specifying the type of curve. Possible
                   value are "hazard", or "survival". The default
                   is "hazard". Only the first words are required, e.g "haz", "su"}
  \item{\dots}{ Other graphical parameters}
}

\value{
     Print a plot of class \code{gcmrec}
}

\seealso{\code{\link{print.gcmrec}}}


\examples{

data(lymphoma)

mod<-gcmrec(Survr(id,time,event)~as.factor(distrib),data=lymphoma,s=1000)

# baseline survivor function

plot(mod)

# baseline hazard function

plot(mod,type="haz")




}

\keyword{ survival }
