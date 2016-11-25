\name{graph.caltimes}
\alias{graph.caltimes}
\title{ Plot of recurrent events }

\description{
  Plots calendar times at succesive recurrences from a data set. Information about effective age and
  categories of covariates are allowed. 
}


\usage{
graph.caltimes(data, var = NULL, effageData = NULL, width = 2, 
        lines = TRUE, sortevents = TRUE, ...)
}

\arguments{
  \item{data}{ data frame containing id, time, event variables and some other covariates}
  \item{var}{ categorical variable  }
  \item{effageData}{ effective age function information }
  \item{width}{ point width }
  \item{lines}{ Are horizontal lines printed? The default is TRUE }
  \item{sortevents}{ Are events sorted? The default is TRUE }
  \item{\dots}{ other graphical parameters }
}

\examples{

# with data in a data frame 
library(survival)
data(bladder2)
bladder2$time<-bladder2$stop-bladder2$start

graph.caltimes(bladder2)

# or data in a list

data(hydraulic)
graph.caltimes(hydraulic)

# We can print some covariate as follows:
graph.caltimes(bladder2,bladder2$rx)
}

\keyword{survival}

