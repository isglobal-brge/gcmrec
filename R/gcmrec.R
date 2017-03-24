"Survr" <-
function (id,time,event) 
{
 
    if(length(unique(id))!=length(event[event==0]))
      {
        stop("Data doesn't match. Every subject must have a censored time")
      }

    if(length(unique(event))>2 | max(event)!=1 | min(event)!=0)
      {
        stop("event must be 0-1")
      }

    ans<-cbind(id,time,event)

    oldClass(ans) <- "Survr"
    invisible(ans)

}

"is.Survr" <-
function(x)
inherits(x, "Survr")

"gcmrec"<-
function (formula, data, effageData = NULL, s, Frailty = FALSE, 
    alphaSeed, betaSeed, xiSeed, tol = 10^(-6), maxit = 100, 
    rhoFunc = "alpha to k", typeEffage = "perfect", maxXi = "Newton-Raphson", 
    se = "Information matrix", cancer = NULL) 
{
    rho.type <- charmatch(rhoFunc, c("Identity","alpha to k"), 
        nomatch = 0)
    if (rho.type == 0) {
        stop("estimator must be 'alpha to k' or 'Identity' ")
    }
    effage.type <- charmatch(typeEffage, c("perfect", "minimal"))
    if (effage.type == 0) {
        stop("typeEffage must be perfect or minimal")
    }
    if (effage.type == 2) {
        effage.type <- 0
    }
    maxXi.type <- charmatch(maxXi, c("Newton-Raphson", "Brent"))
    if (maxXi.type == 0) {
        stop("maxXi must be Newton-Raphson or Brent")
    }
    se.type <- charmatch(se, c("Information matrix", "Jacknife"))
    if (se.type == 0) {
        stop("se must be Information matrix or Jacknife")
    }
    call <- match.call()
    if ((mode(call[[2]]) == "call" && call[[2]][[1]] == as.name("Survr")) || 
        inherits(formula, "Survr")) {
        stop("formula.default(object): invalid formula")
    }
    m <- match.call(expand.dots = FALSE)

    m$s <- m$alphaSeed <- m$betaSeed <- m$xiSeed <- m$tol <- m$maxit <- m$rhoFunc <- m$Frailty <- m$effageData <- m$typeEffage <- m$maxXi <- m$se <- m$cancer <- m$... <- NULL
    Terms <- terms(formula, "strata")
    ord <- attr(Terms, "order")
    if (length(ord) & any(ord != 1)) 
        stop("Interaction terms are not valid for this function")
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.parent())
    n <- nrow(m)
    Y <- model.extract(m, "response")
    if (!is.Survr(Y)) 
        stop("Response must be a survival recurrent object")
    offset <- attr(Terms, "offset")
    tt <- length(offset)
    offset <- if (tt == 0) 
        rep(0, nrow(Y))
    else if (tt == 1) 
        m[[offset]]
    else {
        ff <- m[[offset[1]]]
        for (i in 2:tt) ff <- ff + m[[offset[i]]]
        ff
    }
    ll <- attr(Terms, "term.labels")
    mt <- attr(m, "terms")
    X <- if (!is.empty.model(mt)) 
        model.matrix(mt, m, contrasts)
    if (ncol(X) == 1) 
        stop("model need some covariates")
    cov <- X[, -1]
    if (is.null(effageData)) {
        dataOK <- formatData(Y[, 1], Y[, 2], Y[, 3], cov, effage.type, 
            cancer)
    }
    else {
        dataOK <- formatData.effage(Y[, 1], Y[, 2], Y[, 3], cov, 
            effageData)
    }
    if (max(dataOK$k) > 201) 
        stop("Procedures are implemented for number of \n recurrences less than 200")
    nvar <- ifelse(!is.null(ncol(cov)), ncol(cov), 1)
    if (missing(alphaSeed)) 
        alphaSeed <- 1
    if (missing(betaSeed)) 
        betaSeed <- rep(0, nvar)
    if (missing(xiSeed)) 
        xiSeed <- 1
    if (!Frailty) {
        if (rho.type == 2) {
            ans <- .Fortran("newtraph", as.double(s), as.integer(dataOK$n), 
                as.integer(nvar), as.integer(dataOK$k), as.integer(sum(dataOK$k)), 
                as.double(dataOK$tau), as.double(dataOK$caltimes), 
                as.double(dataOK$gaptimes), as.double(dataOK$censored), 
                as.double(dataOK$intercepts), as.double(dataOK$slopes), 
                as.double(dataOK$lastperrep), as.double(dataOK$perrepind), 
                as.double(dataOK$effagebegin), as.double(dataOK$effage), 
                as.double(dataOK$cov), as.double(alphaSeed), 
                as.double(betaSeed), as.double(rep(1, dataOK$n)), 
                as.double(offset), as.integer(rho.type), as.integer(rep(0, 
                  dataOK$n)), as.double(tol), as.integer(maxit), 
                loglik = as.double(0), estim = as.double(rep(0, 
                  nvar + 1)), info = as.double(matrix(0, nvar + 
                  1, nvar + 1)), search = as.integer(0), kiter = as.integer(0), 
                PACKAGE = "gcmrec")
        }
        if (rho.type == 1) {
            ans <- .Fortran("nrbeta", as.double(s), as.integer(dataOK$n), 
                as.integer(nvar), as.integer(dataOK$k), as.integer(sum(dataOK$k)), 
                as.double(dataOK$tau), as.double(dataOK$caltimes), 
                as.double(dataOK$gaptimes), as.double(dataOK$censored), 
                as.double(dataOK$intercepts), as.double(dataOK$slopes), 
                as.double(dataOK$lastperrep), as.double(dataOK$perrepind), 
                as.double(dataOK$effagebegin), as.double(dataOK$effage), 
                as.double(dataOK$cov), as.double(alphaSeed), 
                as.double(betaSeed), as.double(rep(1, dataOK$n)), 
                as.double(offset), as.integer(rho.type), as.integer(rep(0, 
                  dataOK$n)), as.double(tol), as.integer(maxit), 
                loglik = as.double(0), estim = as.double(rep(0, 
                  nvar)), info = as.double(matrix(0, nvar, 
                  nvar)), search = as.integer(0), kiter = as.integer(0), 
                PACKAGE = "gcmrec")
        }
        if (se.type == 2) {
            jack <- .Fortran("jacknife", as.double(s), as.integer(dataOK$n), 
                as.integer(nvar), as.integer(dataOK$k), as.integer(sum(dataOK$k)), 
                as.double(dataOK$tau), as.double(dataOK$caltimes), 
                as.double(dataOK$gaptimes), as.double(dataOK$censored), 
                as.double(dataOK$intercepts), as.double(dataOK$slopes), 
                as.double(dataOK$lastperrep), as.double(dataOK$perrepind), 
                as.double(dataOK$effagebegin), as.double(dataOK$effage), 
                as.double(dataOK$cov), as.double(ans$estim[1]), 
                as.double(ans$estim[-1]), as.double(rep(1, dataOK$n)), 
                as.double(offset), as.integer(rho.type), as.integer(rep(0, 
                  dataOK$n)), as.double(tol), as.integer(maxit), 
                estiEndJack = as.double(matrix(0, dataOK$n, nvar + 
                  1)), PACKAGE = "gcmrec")
            seJack <- matrix(jack$estiEndJack, dataOK$n, nvar + 
                1)
            JackEst <- apply(seJack, 2, mean)
            diff <- seJack - rep(1, dataOK$n) %*% matrix(JackEst, 
                1, nvar + 1)
            JackEstCov <- (dataOK$n - 1)/dataOK$n * t(diff) %*% 
                diff
        }
        if (ans$search != 1) 
            warning("Algorithm did not converge in: ", maxit, 
                " iterations")
        fit <- list(loglik = ans$loglik)
        fit$coef <- ans$estim
        if (se.type == 1 & rho.type == 2) 
            fit$var <- matrix(ans$info, nvar + 1, nvar + 1)
        if (se.type == 1 & rho.type == 1) 
            fit$var <- matrix(ans$info, nvar, nvar)
        if (se.type == 2) 
            fit$var <- JackEstCov
        fit$n <- ans[[2]]
        fit$nk <- ans[[5]]
        fit$kiter <- ans$kiter
        fit$Xi <- NULL
        fit$frailties <- rep(1, dataOK$n)
        fit$search <- ans$search
    }
    diseff <- sort(unique(dataOK$effage))
    ndiseff <- length(diseff)
    if (Frailty) {
        ans <- .Fortran("estimwithfrailty", s = as.double(s), 
            n = as.integer(dataOK$n), nvar = as.integer(nvar), 
            k = as.integer(dataOK$k), nk = as.integer(sum(dataOK$k)), 
            as.double(dataOK$tau), as.double(dataOK$caltimes), 
            as.double(dataOK$gaptimes), as.double(dataOK$censored), 
            as.double(dataOK$intercepts), as.double(dataOK$slopes), 
            as.double(dataOK$lastperrep), as.double(dataOK$perrepind), 
            as.double(dataOK$effagebegin), as.double(dataOK$effage), 
            as.integer(ndiseff), as.double(diseff), as.double(dataOK$cov), 
            as.double(alphaSeed), as.double(betaSeed), as.double(xiSeed), 
            as.double(rep(1, dataOK$n)), as.double(offset), as.integer(rho.type), 
            as.double(tol), as.integer(maxit), as.integer(maxXi.type), 
            estim = as.double(rep(0, nvar + 2 + dataOK$n)), control = as.integer(c(0, 
                0)), loglik = as.double(0), PACKAGE = "gcmrec")
        if (se.type == 2) {
            jack <- .Fortran("jacknife2", s = as.double(s), n = as.integer(dataOK$n), 
                nvar = as.integer(nvar), k = as.integer(dataOK$k), 
                nk = as.integer(sum(dataOK$k)), as.double(dataOK$tau), 
                as.double(dataOK$caltimes), as.double(dataOK$gaptimes), 
                as.double(dataOK$censored), as.double(dataOK$intercepts), 
                as.double(dataOK$slopes), as.double(dataOK$lastperrep), 
                as.double(dataOK$perrepind), as.double(dataOK$effagebegin), 
                as.double(dataOK$effage), as.integer(ndiseff), 
                as.double(diseff), as.double(dataOK$cov), as.double(ans$estim[1]), 
                as.double(ans$estim[2:(nvar + 1)]), as.double(ans$estim[nvar + 
                  2]), as.double(rep(1, dataOK$n)), as.double(offset), 
                as.integer(rho.type), as.double(tol), as.integer(maxit), 
                as.integer(maxXi.type), estiEndJack = as.double(matrix(0, 
                  dataOK$n, nvar + 2)), PACKAGE = "gcmrec")
            seJack <- matrix(jack$estiEndJack, dataOK$n, nvar + 
                2)
            JackEst <- apply(seJack, 2, mean)
            diff <- seJack - rep(1, dataOK$n) %*% matrix(JackEst, 
                1, nvar + 2)
            JackEstCov <- (dataOK$n - 1)/dataOK$n * t(diff) %*% 
                diff
        }
        fit <- list(loglik = ans$loglik)
        fit$search <- ans$control[1]
        if (fit$search != 1) 
            warning("Algorithm did not converge after: ", maxit, 
                " iterations")
    
        if (rho.type==1)  
          fit$coef <- ans$estim[2:(nvar + 1)]
        if (rho.type==2)  
          fit$coef <- ans$estim[1:(nvar + 1)]

        if (se.type == 1 & rho.type==1) 
            fit$var <- matrix(NA, nvar + 1, nvar + 1)
        if (se.type == 1 & rho.type==2) 
            fit$var <- matrix(NA, nvar + 2, nvar + 2)
        if (se.type == 2) 
            fit$var <- JackEstCov

        fit$n <- ans[[2]]
        fit$nk <- ans[[5]]
        fit$kiter <- ans$control[2]
        fit$Xi <- ans$estim[nvar + 2]
        fit$frailties <- ans$estim[(nvar + 3):(fit$n + nvar + 
            2)]
    }
    fit$method <- "Newton-Raphson"
    fit$rho.type <- rho.type
    if (fit$search == 1) {
        EstLambSurv <- .Fortran("estlambsurv", as.double(s), 
            as.integer(dataOK$n), as.integer(nvar), as.integer(dataOK$k), 
            as.integer(sum(dataOK$k)), as.double(dataOK$tau), 
            as.double(dataOK$caltimes), as.double(dataOK$gaptimes), 
            as.double(dataOK$censored), as.double(dataOK$intercepts), 
            as.double(dataOK$slopes), as.double(dataOK$lastperrep), 
            as.double(dataOK$perrepind), as.double(dataOK$effagebegin), 
            as.double(dataOK$effage), as.integer(ndiseff), as.double(diseff), 
            as.double(dataOK$cov), as.double(fit$coef[1]), as.double(fit$coef[-1]), 
            as.double(fit$frailties), as.double(offset), as.integer(rho.type), 
            Lamb = as.double(rep(0, ndiseff)), DeltaLamb = as.double(rep(0, 
                ndiseff)), Surv = as.double(rep(0, ndiseff)), 
            PACKAGE = "gcmrec")
        fit$Lambda <- EstLambSurv$Lamb
        fit$Survival <- EstLambSurv$Surv
    }
    fit$se.type <- se.type
    fit$diseff <- diseff
    fit$terms <- Terms
    fit$call <- call
    if (length(fit$coef) > 2) {
        mod.X <- X[, -1]
    }
    else {
        mod.X <- matrix(X[, -1])
    }
    if (rho.type==2)
      names(fit$coef) <- c("alpha", colnames(X)[-1])
    if (rho.type==1)
      names(fit$coef) <- c(colnames(X)[-1])

    class(fit) <- "gcmrec"
    fit
}


"plot.gcmrec" <-
function (x, type.plot="surv", ...) 
{
    dostep <- function(x, y) {
        n <- length(x)
        if (n > 2) {
            dupy <- c(TRUE, diff(y[-n]) != 0, TRUE)
            n2 <- sum(dupy)
            xrep <- rep(x[dupy], c(1, rep(2, n2 - 1)))
            yrep <- rep(y[dupy], c(rep(2, n2 - 1), 1))
            list(x = xrep, y = yrep)
        }
        else if (n == 1) 
            list(x = x, y = y)
        else list(x = x[c(1, 2, 2)], y = y[c(1, 1, 2)])
    }

  plot.type <- charmatch(type.plot, c("survival","hazard"), 
        nomatch = 0)
    if (plot.type == 0) {
        stop("estimator must be hazard or survival")
    }

  if(plot.type==1)
   {
    y.lab <- c("Baseline survivor function")
    y <- x$Surv
       plot(dostep(x$diseff, y), type = "l", ylim = c(0, max(y)), 
                xlab = "Time", ylab = y.lab)
   }        

  if(plot.type==2)
   {
    y.lab <- c("Baseline hazard function")
    y <- x$Lam
       plot(dostep(x$diseff, y), type = "l", ylim = c(0, max(y)), 
                xlab = "Time", ylab = y.lab)
   }        

    
    return(invisible())
}

"lines.gcmrec"<-
function (x, type.plot="surv", ...) 
{
    dostep <- function(x, y) {
        n <- length(x)
        if (n > 2) {
            dupy <- c(TRUE, diff(y[-n]) != 0, TRUE)
            n2 <- sum(dupy)
            xrep <- rep(x[dupy], c(1, rep(2, n2 - 1)))
            yrep <- rep(y[dupy], c(rep(2, n2 - 1), 1))
            list(x = xrep, y = yrep)
        }
        else if (n == 1) 
            list(x = x, y = y)
        else list(x = x[c(1, 2, 2)], y = y[c(1, 1, 2)])
    }

  plot.type <- charmatch(type.plot, c("survival","hazard"), 
        nomatch = 0)
    if (plot.type == 0) {
        stop("estimator must be hazard or survival")
    }

  if(plot.type==1)
   {
    y.lab <- c("Baseline survivor function")
    y <- x$Surv
    lines(dostep(x$diseff, y), ...)        
   } 

  if(plot.type==2)
   {
    y.lab <- c("Baseline hazard function")
    y <- x$Lam
    lines(dostep(x$diseff, y), ...)        
   }
    return(invisible())
}


"graph.caltimes"<-
function(data, var=NULL, effageData=NULL, width = 2, lines = TRUE, sortevents = TRUE, ...)
{

if (!is.data.frame(data))
  {
    data<-List.to.Dataframe(data)
  }

id.unique<-unique(data$id)
nsubjs<-length(id.unique)


if (!is.null(effageData))
  {
    perfect <- lapply(effageData, function(x) x[["perrepind"]])
  }
else
  {
    perfect <- rep(1,nrow(data))
  }

tauval<-NULL
for(i in 1:nsubjs)
 {
   tauval[i]<-sum(data$time[data$id==id.unique[i]])
 }

oldcex <- par("cex")
par(cex = 1.25)


plot(c(0,max(tauval)), c(1,nsubjs), type = "n", xlab = "Calendar Time", ylab = "Subject",...)
for (i in 1:nsubjs)
 {
  thiscaltimes <- cumsum(data$time[data$id==id.unique[i] & data$event==1])
  thisperfect <- perfect[[i]]
  thistau <- tauval[i]

  if ((!is.null(var)) && (length(unique(var))<=5))
    {  
      thiscovvals <- var[data$id==id.unique[i]]
      colors <- thiscovvals[1] 
    }
  else 
    {
     colors<-1
    } 

  pchs <- ifelse(thisperfect == 0, 2, 1)
  points(thiscaltimes, rep(i, length(thiscaltimes)), col = colors, pch = pchs, lwd = width) 
  points(thistau, i, col="darkgreen",pch=4, lwd = width)
  if (lines) abline(h = i, lty = 3, lwd = .5)
 }
par(cex = oldcex)
}



"addCenTime"<-
function(datin, id=1, time=2, event=3)
{
 #
 # datin should have id,time,event variables
 # 
 
 ids <- unique(datin[,id])
 datout <- NULL
 for (i in ids)
 {
  thissubj <- datin[datin[,id] == i, , drop=FALSE]
  nrecs <- nrow(thissubj)
  if (thissubj[nrecs, event] == 1)
  {
   newrec <- thissubj[nrecs,]
   newrec[, time] <- 0
   newrec[, event] <- 0
   thissubj <- rbind(thissubj, newrec)
  }
  datout <- rbind(datout, thissubj)
 }
 datout
}


"List.to.Dataframe"<-
function (data) 
{
    time.cen <- data$subject[[1]]$tau - sum(data$subject[[1]]$gaptime)
    time <- c(data$subject[[1]]$gaptime[-1], time.cen)
    event <- c(rep(1, data$subject[[1]]$k - 1), 0)
    id <- rep(1, data$subject[[1]]$k)
    dataEnd <- cbind(id = id, time = time, event = event, t(data$subject[[1]]$cov))
    i <- 2
    while (i <= data$n) {
        time.cen <- data$subject[[i]]$tau - sum(data$subject[[i]]$gaptime)
        time <- c(data$subject[[i]]$gaptime[-1], time.cen)
        event <- c(rep(1, data$subject[[i]]$k - 1), 0)
        id <- rep(i, data$subject[[i]]$k)
        temp <- cbind(id = id, time = time, event = event, 
            t(data$subject[[i]]$cov))
        dataEnd <- rbind(dataEnd, temp)
        i <- i + 1
    }
    dataEnd <- data.frame(dataEnd)
    nvar <- nrow(data$subject[[1]]$cov)
    names.Cov <- paste("covar.", 1:nvar, sep = "")
    names(dataEnd) <- c(names(dataEnd)[1:3], names.Cov)
    dataEnd
}


"formatData" <-
function (id, time, event, covariates, parameffage,cancer) 
{
    covariates <- data.frame(covariates)
    n <- length(unique(id))
    id.distinct <- unique(id)
    temp <- formatData.i(id[id == id.distinct[1]], time[id == 
        id.distinct[1]], event[id == id.distinct[1]], covariates[id == 
        id.distinct[1], ], parameffage,cancer[id == id.distinct[1]])


    k <- temp$k
    tau <- temp$tau
    caltimes <- temp$caltimes
    gaptimes <- temp$gaptimes
    censored <- temp$censored
    intercepts <- temp$intercepts
    slopes <- temp$slopes
    lastperrep <- temp$lastperrep
    perrepind <- temp$perrepind
    effagebegin <- temp$effagebegin
    effage <- temp$effage
    covariate <- temp$covariate
    for (i in 2:n) {
        temp <- formatData.i(id[id == id.distinct[i]], time[id == 
            id.distinct[i]], event[id == id.distinct[i]], covariates[id == 
            id.distinct[i], ], parameffage, cancer[id == id.distinct[i]])

        k <- c(k, temp$k)
        tau <- c(tau, temp$tau)
        caltimes <- c(caltimes, temp$caltimes)
        gaptimes <- c(gaptimes, temp$gaptimes)
        censored <- c(censored, temp$censored)
        intercepts <- c(intercepts, temp$intercepts)
        slopes <- c(slopes, temp$slopes)
        lastperrep <- c(lastperrep, temp$lastperrep)
        perrepind <- c(perrepind, temp$perrepind)
        effagebegin <- c(effagebegin, temp$effagebegin)
        effage <- c(effage, temp$effage)
        covariate <- cbind(covariate, temp$covariate)
    }
    ans <- list(n = n, k = k, tau = tau, caltimes = caltimes, 
        gaptimes = gaptimes, censored = censored, intercepts = intercepts, 
        slopes = slopes, lastperrep = lastperrep, perrepind = perrepind, 
        effagebegin = effagebegin, effage = effage, covariate = covariate)
    return(ans)
}


"formatData.i" <-
function (id, time, event, covariates, parameffage, cancer = NULL) 
{
    covariates <- data.frame(covariates)
    k <- length(id)
    tau <- sum(time)
    caltimes <- cumsum(c(0, time[event == 1]))
    gaptimes <- c(0, time[event == 1])
    censored <- time[event == 0]
    lastperrep <- c(1)
    eff <- generlmi(parameffage)
    if (is.null(cancer)) {
        intercepts <- eff$intercept
        slopes <- eff$slope
        effagebegin <- eff$intercept
        effageend <- eff$intercept
        perrepind <- c(1)
    }
    else {
        if (cancer[1] == "CR") {
            intercepts <- 0
            slopes <- 1
            effagebegin <- intercepts
            effageend <- intercepts
            perrepind <- c(1)
        }

#
# This is in the case of not CR after first treatment. We must change A_0=0
#
        if (cancer[1] != "CR") {
            intercepts <- 0
            slopes <- 1
            effagebegin <- intercepts
            effageend <- intercepts
            perrepind <- c(1)
        }
    }
    if (k >= 2) {
        for (i in 2:k) {
            eff <- generlmi(parameffage)
            intercepts <- c(intercepts, eff$intercept)
            slopes <- c(slopes, eff$slope)
            if (is.null(cancer)) {
                perrepind <- c(perrepind, eff$perrepind)
                if (eff$perrepind == 1) {
                  lastperrep <- c(lastperrep, i)
                }
                if (eff$perrepind == 0) {
                  lastperrep <- c(lastperrep, lastperrep[i - 
                    1])
                }
                effagebegin <- c(effagebegin, intercepts[i - 
                  1] + slopes[i] * (caltimes[i] - caltimes[lastperrep[i]]))
                effageend <- c(effageend, intercepts[i - 1] + 
                  slopes[i] * (caltimes[i] - caltimes[lastperrep[i - 
                    1]]))
            }
            else {
                if (cancer[i] == "CR") {
                  lastperrep <- c(lastperrep, i)
                  perrepind <- c(perrepind, 1)
                }
                if (cancer[i] != "CR") {
                  lastperrep <- c(lastperrep, lastperrep[i - 
                    1])
                  perrepind <- c(perrepind, 0)
                }
 
if(cancer[i] == "CR")
 {
  effagebegin <- c(effagebegin, 0)
  effageend <- c(effageend,gaptimes[i])
 }

if (cancer[i] == "SD")
 { 
  control<-effageend[i-1]
  effagebegin <- c(effagebegin, control+gaptimes[i])
  effageend <- c(effageend, control+gaptimes[i])
 }

if (cancer[i] == "PR")
 { 
  control<-effageend[i-1]
  effagebegin <- c(effagebegin, control+0.5*gaptimes[i])
  effageend <- c(effageend, control+0.5*gaptimes[i])
 }




            }
        }
    }
    covariate <- t(covariates)
    ans <- list(k = k, tau = tau, caltimes = caltimes, gaptimes = gaptimes, 
        censored = censored, intercepts = intercepts, slopes = slopes, 
        lastperrep = lastperrep, perrepind = perrepind, effagebegin = effagebegin, 
        effage = effageend, covariate = matrix(covariate, ncol(covariates), 
            k))
    return(ans)
}



"generlmi"<-
function(perrep)
{
ll <- 0
mm <- 1
ii <- rbinom(1, 1, perrep)
return(list(intercept = ll, slope = mm, perrepind = ii))
}



"formatData.effage"<-
function (id, time, status, covariates,effageData) 
{

# effageData: 
#    list containing effective age information:
#      effageData$n effageData$subject
#
#    effageData$subject:  
#       list containing:
#               intercepts, slopes, lastperrep, 
#               perrepind, effagebegin, effage




# some controls
if (effageData$n!=length(unique(id))) stop('data does not match')

# ...........


id.unique<-unique(id)
n<-length(id.unique)


covariates<-data.frame(covariates)

kk<-table(id)


tau<-sum(time[id==id.unique[1]])
caltimes<-cumsum(c(0,time[id==id.unique[1] & status==1]))
gaptimes<-c(0,time[id==id.unique[1] & status==1])
censored<-time[id==id.unique[1] & status==0]
covariate<-t(covariates[id==id.unique[1],])


intercepts<-effageData$subject[[1]]$intercepts
slopes<-effageData$subject[[1]]$slopes
lastperrep<-effageData$subject[[1]]$lastperrep
perrepind<-effageData$subject[[1]]$perrepind
effagebegin<-effageData$subject[[1]]$effagebegin
effage<-effageData$subject[[1]]$effage


for (i in 2:n)
{
       intercepts<-c(intercepts,effageData$subject[[i]]$intercepts)
       slopes<-c(slopes,effageData$subject[[i]]$slopes)
       lastperrep<-c(lastperrep,effageData$subject[[i]]$lastperrep)
       perrepind<-c(perrepind,effageData$subject[[i]]$perrepind)
       effagebegin<-c(effagebegin,effageData$subject[[i]]$effagebegin)
       effage<-c(effage,effageData$subject[[i]]$effage)

       tau<-c(tau,sum(time[id==id.unique[i]]))
       caltimes<-c(caltimes,cumsum(c(0,time[id==id.unique[i] & status==1])))
       gaptimes<-c(gaptimes,c(0,time[id==id.unique[i] & status==1]))
       censored<-c(censored,time[id==id.unique[i] & status==0])
       covariate<-cbind(covariate,t(covariates[id==id.unique[i],]))
}

ans<-list(n=n, k=kk, tau=tau, caltimes=caltimes, gaptimes=gaptimes, censored=censored, 
intercepts=intercepts, slopes=slopes, lastperrep=lastperrep, perrepind = perrepind, 
effagebegin=effagebegin, effage=effage, covariate=covariate)

return(ans)

}


"summary.gcmrec" <-
function(object,level=.95, len=6, d=2, lab="hr", ...)
    {

      x <- object
      if (!inherits(x, "gcmrec")) 
         stop("Invalid data")
      
      z<-abs(qnorm((1-level)/2))
      co <- x$coef[-1]
      se <- sqrt(diag(x$var))[-1]
      or <- exp(co)
      li <- exp(co-z * se)
      ls <- exp(co+z * se)
      r <- cbind(or, li, ls)

      dimnames(r) <- list(names(co), c(lab, paste(level*100,"%",sep=""), "C.I."))
      
      n<-r
   
      dd <- dim(n)
      n[n > 999.99] <- Inf
      a <- formatC(n, d, len,format="f")

      dim(a) <- dd
      if(length(dd) == 1){
          dd<-c(1,dd)
          dim(a)<-dd
          lab<-" "
          }
      else
          lab <- dimnames(n)[[1]]
      
      mx <- max(nchar(lab)) + 1
      cat(paste(rep(" ",mx),collapse=""),paste("   ",dimnames(n)[[2]]),"\n")
      for(i in (1):dd[1]) {
          lab[i] <- paste(c(rep(" ", mx - nchar(lab[i])), lab[i]),collapse = "")
      cat(lab[i], a[i, 1], "(", a[i, 2], "-", a[i, 3], ") \n")
      }
      
}



"print.gcmrec" <-
function (x, digits = max(options()$digits - 4, 3), ...) 
{
    if (!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl)
        cat("\n")
    }
    if (!is.null(x$fail)) {
        cat(" gmcrec failed.", x$fail, "\n")
        return()
    }
    savedig <- options(digits = digits)
    on.exit(options(savedig))
   
    label.Rho<-ifelse(x$rho.type==1,"Identity", "Alpha to k")

    if(x$rho.type==1)
       {
         x$coef<-c(1,x$coef)
       }
    coef <- x$coef[-1]

    se <- sqrt(diag(x$var))
    if(x$rho.type==1)
       {
         se<-c(NA,se)
       }

    if (is.null(x$Xi)) 
      se <- se[-1]
    if (!is.null(x$Xi)) 
      {
        n.Xi<-nrow(x$var)
        se <- se[-c(1,n.Xi)]
      }
    if (is.null(coef) | is.null(se)) 
        stop("Input is not valid")
    tmp <- cbind(coef, exp(coef), se, coef/se, signif(1 - pchisq((coef/se)^2, 
        1), digits - 1))

    if (x$se.type==1)
      dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)", 
        "se(coef)", "z", "p"))
    if (x$se.type==2)
      dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)", 
        "se(coef) Jacknife", "z", "p"))
       
    cat("\n")
    prmatrix(tmp)
    if (is.null(x$df)) 
        df <- sum(!is.na(coef))
    else df <- round(sum(x$df), 2)
    cat("\n")
    cat("  General class model parameter estimates", "\n")
    cat("    rho function: ", label.Rho, "\n")
    alpha <- x$coef[1]

    if (x$rho.type==2)
       se <- sqrt(x$var[1,1])
    if (x$rho.type==1)
       se <- "--"

    if (x$se.type==1)
      cat("      alpha (s.e.): ", alpha, " (", se, ")", "\n", sep = "")
    if (x$se.type==2)
      cat("      alpha (s.e. Jacknife): ", alpha, " (", se, ")", "\n", sep = "")
    
    if (!is.null(x$Xi))
      {
        se <- sqrt(x$var[n.Xi,n.Xi])  
        cat("    Frailty parameter, Xi (s.e. Jacknife): ", x$Xi," (", se, ")", "\n")
      }
    cat(" \n")

    if (is.null(x$Xi))
     {
       cat(paste("  log-likelihood=", round(x$loglik, 2)), "\n")
       cat("  n=", x$n, "\n")
       cat("  n times=", x$nk, "\n")
       cat("  number of iterations: ", x$kiter, " Newton-Raphson \n")
       cat("\n")
     }

    if (!is.null(x$Xi))
     {
       cat(paste("  Marginal log-likelihood=", round(x$loglik, 2)), "\n")
       cat("  n=", x$n, "\n")
       cat("  n times=", x$nk, "\n")
       cat("  number of iterations: ", x$kiter, " EM steps \n")
       cat("\n")
     }
    invisible()
}

