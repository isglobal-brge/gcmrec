%%%%%%%%%%%%%%
%%%  gcmrec
%%%%%%%%%%%%%%

\name{gcmrec}
\alias{gcmrec}
\alias{Survr}
\alias{is.Survr}
\title{ General Class of Models for recurrent event data }
\description{
  Fits the parameters for the general semiparametric model for recurrent events
  proposed by Peña and Hollander (2004). This class of models
  incorporates an effective age function which encodes the changes that occur after each event
  occurrence such as the impact of an intervention, it allows for the modeling of the impact of accumulating
  event occurrences on the unit, it admits a link function in which the effect of possibly
  time-dependent covariates are incorporated, and it allows the incorporation of unobservable
  frailty components which induce dependencies among the inter-event times for each unit.
}
\usage{
gcmrec(formula, data, effageData = NULL, s, Frailty = FALSE, 
    alphaSeed, betaSeed, xiSeed, tol = 10^(-6), maxit = 100, 
    rhoFunc = "alpha to k", typeEffage = "perfect", 
    maxXi = "Newton-Raphson", se = "Information matrix", 
    cancer = NULL)
}
\arguments{
  \item{formula}{A formula object. If a formula object is supplied it must
          have a Survr object as the response on the left of the '~'
          operator and a term on the right. Covariates are needed.}
  \item{data}{
                     A data frame in which to interpret the variables named 
                     in the formula. This data frame must contain the
                     variables called "id","time" and "event" for subject identification,
                     time of interocurrence, and censored status (coded 1: event, 0:censored),
                     respectively. Furthermore, we can have some covariates. 

                     Alternatively, it can also be a list containing the elements "n" and "subjects". Number of 
                     subjects must be recorded in "n". The element "subject" must have the 
                     following elements: subj, k, tau, caltimes, gaptimes, 
                     intercepts, slopes, lastperrep, perrepind, effage, effagebegin, 
                     and covariate including this information:

              \describe{
               \item{subj:}{Subject number or identificator. }
               \item{k:}{Number of recurrences (time 0 must be included). }
               \item{tau:}{Administrative time, time of study termination. }
               \item{caltimes:}{Calendar times at each recurrence (time 0 must be included). } 
               \item{gaptimes:}{Gap times at each recurrence (time 0 must be included). }
               \item{intercepts:}{Intercept value for the effect after each recurrence. }
               \item{slopes:}{Slope value for the effect after each each recurrence.}             
               \item{lastperrep:}{Element from Brown and Proschan minimal repair model. }
               \item{perrepind:}{Element from Brown and Proschan minimal repair model. }
               \item{effagebegin:}{Initial value for effective age. }
               \item{effage:}{Effective age after each recurrence. }
               \item{covariate:}{covariate value at each recurrence.}
             }
               See either GeneratedData or hydraulic data sets as an example.  
            }  

               
  \item{effageData}{List containing the information about effective age. The list must have 
                    the elements described in the option 2 of data argument. If NULL we 
                    generate these elements under perfect repair model or minimal repair one depending on the 
                    'typeEffage' argument (see below).    }
  \item{s}{A selected calendar time. }
  \item{Frailty}{Logical value. Is model with frailties fitted? If so parameters for General 
                 Class of Models with frailty component are estimated.  }
  \item{alphaSeed}{Seed value for \eqn{\alpha}{%
                              alpha}. }
  \item{betaSeed}{Seed value for \eqn{\beta}{%
                               beta}. }
  \item{xiSeed}{Seed value for \eqn{\xi}{%
                                xi}. }
  \item{tol}{Tolerance for maximization procedures. }
  \item{maxit}{Maximum number of iterations in maximization procedures.}
  \item{rhoFunc}{A character string specifying the effects attributable to the 
          accumulating event occurrences, \eqn{\rho(k;\alpha)}{% 
          rho(k;alpha)}. Possible values are "Identity" for \eqn{\rho(k;\alpha)=1}{%
          rho(k;alpha)=1} or "alpha to k" for \eqn{\rho(k;\alpha)=\alpha^k}{%
          rho(k;alpha)=alpha^k}. The default is "alpha to k". Only the first words 
          are required, e.g "Id","a". Future versions will include other functions such as Markovian model
          for tumor occurrences, \eqn{\rho(k;\alpha)=\alpha-k+1}{%
          rho(k;alpha)=alpha-k+1} proposed by Gail et al. (1980).} 
  \item{typeEffage}{Effective age function. Possible value are "perfect" or "minimal" for
                 perfect repair model or minimal repair model, respectively. The default
                 is "perfect". Only the first words are required, e.g "p","m"}       
  \item{maxXi}{Maximization method for marginal likelihood with respect to \eqn{\xi}{%
                xi}. Possible values are "Newton-Raphson" for Newton-Raphson maximization 
                procedure or "Brent" for Brent's method maximization in one dimension. 
                The default value is "Newton-Raphson". Only the first words are required,
                e.g. "N","B"}  
  \item{se}{Standard errors of parameters. Possible values are 'Information matrix' or 'Jacknife' for inverse of the partial 
                likelihood information matrix or jacnife estimates, respectively.}
  \item{cancer}{Effective age for fitting a cancer model proposed by Gonzalez et al (2005). This variable contains the information
                      of the effect of treatments administred after cancer relapses coded as "CR", "PR" or "SD" depending on if complete,
                      partial, or null response (stable disease) is achieved. See lymphoma data set as an example. }                                                   
}

\details{
  Estimation with frailties are implemented using expectation-maximization (EM) algorithm. In this
  procedure, we need to maximize the marginal likelihood with respect to \eqn{\xi}{%
  xi}. This maximization is a one-dimensional maximization without derivarives. First we bracket the 
  maximizing value, and the we obtain it using Brent's method in one dimension. When we implement 
  this algorithm, we re-parameterize \eqn{\xi}{%
  xi} using \eqn{\xi^*=log(\xi)}{%
  xi*=log(xi)} to alleviate the problem of getting negative estimates for \eqn{\xi}{%
  xi}. Iteration is terminated when successive values of \eqn{\xi/(1+\xi)}{%
  xi/(1+xi)} differ by no more than the "tol" parameter. Maybe estimation under frailty 
  model can be not too fast.   
  
}

\value{
     a gcmrec object. Methods defined for gcmrec objects are
     provided for print and plot.
}

\references{
            Pe\~na, E. and M. Hollander (2004). Mathematical Reliability: An Expository Perspective,
            Chapter 6. Models for Recurrent Events in Reliability and Survival Analysis, pp. 105-123.
            Kluwer Academic Publishers.
            
            M. Gail, T Santner, and C Brown (1980). An analysis of comparative carfinogenesis 
            experiments based on multiple times to tumor. \emph{Biometrics} \bold{36}, 255-266.          
         
            JR Gonzalez, E Peña, E Slate (2005). Modelling treatment effect after cancer relapses, with application to
            recurrences in indolent non-Hodgkin's lymphomas. \emph{Stat Med}, submitted.

            R. Brent. \emph{Algorithms for Minimization Without Derivatives}. Prentice-Hall,
            New York, 1973.
   }

 

\examples{



###################################
## Models using different data formats
###################################

#
#    Data input as a data frame
#

#   We use the well-known bladder cancer data set from survival package

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

bladderOK<-addCenTime(bladder2) 

# Now, we can fit the model using this new data set:

gcmrec(Survr(id,time,event)~rx+size+number,data=bladderOK,s=2060)



#
#    Data as a list. See either GeneratedData or hydraulic data 
#                    sets as an example.
#

#
# We can fit the model by transforming our data in a data frame 
# using "List.to.Dataframe" function:
#

data(hydraulic)
hydraulicOK<-List.to.Dataframe(hydraulic)
gcmrec(Survr(id,time,event)~covar.1+covar.2,data=hydraulicOK,s=4753)



#
# Our model allows us to incorporate effective age information
#
# To illustrate this example, we will use a simulated data set generated 
# under the minimal repair model with probability of perfect repair equal to 0.6
#
# As we have the data in a list, first we need to obtain a data frame containing
# the time, event, and covariates information:
#

data(GeneratedData)
temp<-List.to.Dataframe(GeneratedData)

# then, we can fit the model incorporating the information about the effective 
#        age in the effageData argument:

gcmrec(Survr(id,time,event)~covar.1+covar.2, data=temp, 
                    effageData=GeneratedData, s=100)




#####################################################################
##  How to fit minimal or perfect repair models, with and without frailties
#####################################################################

# Model with frailties 

mod.Fra<-gcmrec(Survr(id,time,event)~rx+size+number,data=bladderOK,s=2060,Frailty=TRUE)
print(mod.Fra)


# effective age function: perfect repair and minimal repair models 
# (models without frailties)

data(readmission)

# perfect
mod.per<-gcmrec(Survr(id,time,event)~as.factor(dukes),data=readmission,
    s=3000,typeEffage="per")
print(mod.per)

# minimal
mod.min<-gcmrec(Survr(id,time,event)~as.factor(dukes),data=readmission,
   s=3000,typeEffage="min")
print(mod.min)


#####################################################################
##  How to fit models with \rho function equal to identity 
#####################################################################

data(lymphoma)

gcmrec(Survr(id, time, event) ~ as.factor(distrib), 
    data = lymphoma, s = 1000, Frailty = TRUE, rhoFunc = "Ident")


#####################################################################
##  How to fit cancer model
#####################################################################


mod.can<-gcmrec(Survr(id,time,event)~as.factor(distrib), data=lymphoma,
     s=1000, Frailty=TRUE,  cancer=lymphoma$effage)

#  standard errors can be obtained by adding se="Jacknife".
# This procedure can be very time consuming

 }

\keyword{survival}






