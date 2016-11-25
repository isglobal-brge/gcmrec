\name{addCenTime}
\alias{addCenTime}
\title{ Add censored time equal to 0}
\description{
  Add a new line to the dataframe with a censored time equal to 0 when
  the end of follow-up coincides to the last occurrence
}
\usage{
addCenTime(datin)
}

\arguments{
  \item{datin}{Dataframe containing id, time and event variables. 
               Another covariates are allowed}
}
\value{
  A data frame with an added line (censored time equal to 0) for those subjects where 
  the end of follow-up coincides to the last occurrence  

}


\examples{
data(bladder, package="survival")

# we compute the interocurrence time
bladder2$time<-bladder2$stop-bladder2$start


# If we execute:
#     gcmrec(Survr(id,time,event)~rx+size+number,data=bladder2,s=2060)

# We will obtain the following error message:
#   Error in Survr(id, time, event) : Data doesn't match... 


# This means that we have some patients without right-censoring time. So,
# we understand that the last event coincides with the end of study. 
# Consequently,we need to add a line with time 0 and status value equal 
# to 0, too. To do so, we can use the function "addCenTime" as follows:


# for example:
#   bladder2[bladder2$id==12,]

#     id rx number size start stop event enum time
#  45 12  1      1    1     0    3     1    1    3
#  46 12  1      1    1     3   16     1    2   13
#  47 12  1      1    1    16   23     1    3    7

# there is no censored time for 12th patient. So, if we execute

bladderOK<-addCenTime(bladder2) 


# we would get

#      id rx number size start stop event enum time
#  45  12  1      1    1     0    3     1    1    3
#  46  12  1      1    1     3   16     1    2   13
#  47  12  1      1    1    16   23     1    3    7
#  471 12  1      1    1    16   23     0    3    0

}
\keyword{survival}% at least one, from doc/KEYWORDS

