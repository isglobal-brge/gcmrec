# Internal per-subject data preparation for the C++ core.
#
# formatData() aggregates the information of all subjects into the
# "flat" vectors (concatenated by subject) consumed by the estimation
# core; formatData.i() handles one subject. With external effageData,
# formatData.effage() is used instead.

# Aggregates the data of all subjects into flat format.
#
# id, time, event: vectors of the Survr response (ordered by subject).
# covariates: matrix/data.frame of covariates, one row per record.
# parameffage: probability of perfect repair (1 = perfect, 0 = minimal)
#   used by generlmi().
# cancer: optional CR/PR/SD vector for the cancer model.
formatData <- function(id, time, event, covariates, parameffage, cancer) {
  covariates <- data.frame(covariates)
  id_distinct <- unique(id)

  per_subject <- lapply(id_distinct, function(i) {
    sel <- id == i
    formatData.i(id[sel], time[sel], event[sel],
                 covariates[sel, , drop = FALSE], parameffage, cancer[sel])
  })

  pull <- function(field) unlist(lapply(per_subject, `[[`, field),
                                 use.names = FALSE)

  list(
    n = length(id_distinct),
    k = pull("k"),
    tau = pull("tau"),
    caltimes = pull("caltimes"),
    gaptimes = pull("gaptimes"),
    censored = pull("censored"),
    intercepts = pull("intercepts"),
    slopes = pull("slopes"),
    lastperrep = pull("lastperrep"),
    perrepind = pull("perrepind"),
    effagebegin = pull("effagebegin"),
    effage = pull("effage"),
    covariate = do.call(cbind, lapply(per_subject, `[[`, "covariate"))
  )
}

# Prepares the data of a single subject: calendar times, effective age
# at the start and end of each period, and Brown-Proschan repair model
# indicators. With `cancer`, the effective age is given by the response
# to treatment (CR: reset to 0; PR: advances half the gap; SD: advances
# the full gap).
formatData.i <- function(id, time, event, covariates, parameffage,
                         cancer = NULL) {
  covariates <- data.frame(covariates)
  k <- length(id)
  tau <- sum(time)
  caltimes <- cumsum(c(0, time[event == 1]))
  gaptimes <- c(0, time[event == 1])
  censored <- time[event == 0]
  lastperrep <- 1

  eff <- generlmi(parameffage)
  intercepts <- if (is.null(cancer)) eff$intercept else 0
  slopes <- if (is.null(cancer)) eff$slope else 1
  effagebegin <- intercepts
  effageend <- intercepts
  perrepind <- 1

  if (k >= 2) {
    for (i in 2:k) {
      eff <- generlmi(parameffage)
      intercepts <- c(intercepts, eff$intercept)
      slopes <- c(slopes, eff$slope)

      if (is.null(cancer)) {
        perrepind <- c(perrepind, eff$perrepind)
        lastperrep <- c(lastperrep,
                        if (eff$perrepind == 1) i else lastperrep[i - 1])
        effagebegin <- c(effagebegin, intercepts[i - 1] +
                           slopes[i] * (caltimes[i] -
                                          caltimes[lastperrep[i]]))
        effageend <- c(effageend, intercepts[i - 1] +
                         slopes[i] * (caltimes[i] -
                                        caltimes[lastperrep[i - 1]]))
      } else {
        if (cancer[i] == "CR") {
          lastperrep <- c(lastperrep, i)
          perrepind <- c(perrepind, 1)
          effagebegin <- c(effagebegin, 0)
          effageend <- c(effageend, gaptimes[i])
        } else {
          lastperrep <- c(lastperrep, lastperrep[i - 1])
          perrepind <- c(perrepind, 0)
          # PR advances half the gap; SD (no response) the full gap
          step <- if (cancer[i] == "PR") 0.5 * gaptimes[i] else gaptimes[i]
          effagebegin <- c(effagebegin, effageend[i - 1] + step)
          effageend <- c(effageend, effageend[i - 1] + step)
        }
      }
    }
  }

  list(
    k = k,
    tau = tau,
    caltimes = caltimes,
    gaptimes = gaptimes,
    censored = censored,
    intercepts = intercepts,
    slopes = slopes,
    lastperrep = lastperrep,
    perrepind = perrepind,
    effagebegin = effagebegin,
    effage = effageend,
    covariate = matrix(t(covariates), ncol(covariates), k)
  )
}

# Generates the intercept, slope and perfect repair indicator for one
# period. perrep is the probability of perfect repair: with 1 (perfect
# model) or 0 (minimal) it is deterministic; intermediate values sample
# from the Brown-Proschan model (uses the session RNG).
generlmi <- function(perrep) {
  list(intercept = 0, slope = 1, perrepind = rbinom(1, 1, perrep))
}

# Variant of formatData() for when the effective age information is
# supplied by the user in effageData (a list with n and subject, each
# subject with intercepts, slopes, lastperrep, perrepind, effagebegin
# and effage).
formatData.effage <- function(id, time, status, covariates, effageData) {
  if (effageData$n != length(unique(id))) stop("data does not match")

  covariates <- data.frame(covariates)
  id_unique <- unique(id)
  n <- length(id_unique)

  per_subject <- lapply(seq_len(n), function(i) {
    sel <- id == id_unique[i]
    ev <- sel & status == 1
    list(
      k = sum(sel),
      tau = sum(time[sel]),
      caltimes = cumsum(c(0, time[ev])),
      gaptimes = c(0, time[ev]),
      censored = time[sel & status == 0],
      covariate = t(covariates[sel, , drop = FALSE])
    )
  })

  pull <- function(field) unlist(lapply(per_subject, `[[`, field),
                                 use.names = FALSE)
  pull_eff <- function(field) {
    unlist(lapply(effageData$subject, `[[`, field), use.names = FALSE)
  }

  list(
    n = n,
    k = pull("k"),
    tau = pull("tau"),
    caltimes = pull("caltimes"),
    gaptimes = pull("gaptimes"),
    censored = pull("censored"),
    intercepts = pull_eff("intercepts"),
    slopes = pull_eff("slopes"),
    lastperrep = pull_eff("lastperrep"),
    perrepind = pull_eff("perrepind"),
    effagebegin = pull_eff("effagebegin"),
    effage = pull_eff("effage"),
    covariate = do.call(cbind, lapply(per_subject, `[[`, "covariate"))
  )
}
