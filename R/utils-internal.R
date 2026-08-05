# Internal utilities shared by the plotting methods.

# Converts a curve evaluated at discrete points into the coordinates of
# a step function (so it is drawn without interpolating across jumps).
dostep <- function(x, y) {
  n <- length(x)
  if (n > 2) {
    # Drop duplicated interior points along constant stretches
    dupy <- c(TRUE, diff(y[-n]) != 0, TRUE)
    n2 <- sum(dupy)
    xrep <- rep(x[dupy], c(1, rep(2, n2 - 1)))
    yrep <- rep(y[dupy], c(rep(2, n2 - 1), 1))
    list(x = xrep, y = yrep)
  } else if (n == 1) {
    list(x = x, y = y)
  } else {
    list(x = x[c(1, 2, 2)], y = y[c(1, 1, 2)])
  }
}

# Extracts the requested baseline curve from a gcmrec fit as a data
# frame, with a pointwise confidence band if requested.
#
# The variance of Lambda0 is Aalen's (sum dN / S0^2). The limits are
# built on the log scale for Lambda (which guarantees positive limits)
# and on the log-log scale for the survivor function (limits inside
# [0, 1]), the default transformation of survival::survfit.
baseline_data <- function(x, type.plot, conf.int = TRUE, level = 0.95) {
  plot.type <- charmatch(type.plot, c("survival", "hazard"), nomatch = 0)
  if (plot.type == 0) {
    stop("estimator must be hazard or survival")
  }
  if (is.null(x$Survival)) {
    stop("baseline functions not available (the fit did not converge)")
  }

  z <- abs(qnorm((1 - level) / 2))
  var.lambda <- if (is.null(x$varLambda)) rep(NA_real_, length(x$diseff)) else
    x$varLambda

  if (plot.type == 1) {
    y <- x$Survival
    ylab <- "Baseline survivor function"
    # Log-log scale: S^exp(-+ z se / log(S))
    se.loglog <- sqrt(var.lambda) / pmax(-log(pmax(y, 1e-12)), 1e-12)
    lower <- y^exp(z * se.loglog)
    upper <- y^exp(-z * se.loglog)
  } else {
    y <- x$Lambda
    ylab <- "Baseline cumulative hazard"
    # Log scale on Lambda
    se.log <- sqrt(var.lambda) / pmax(y, 1e-12)
    lower <- y * exp(-z * se.log)
    upper <- y * exp(z * se.log)
  }

  out <- data.frame(time = x$diseff, estimate = y)
  if (conf.int && !all(is.na(var.lambda))) {
    out$lower <- pmax(lower, 0)
    out$upper <- if (plot.type == 1) pmin(upper, 1) else upper
    # The first point (time 0) carries no accumulated information
    out$lower[1] <- out$estimate[1]
    out$upper[1] <- out$estimate[1]
  }
  attr(out, "ylab") <- ylab
  out
}

# Expands a curve data frame to step coordinates. All value columns
# (estimate and limits) are expanded with the same index, so the band
# follows the step exactly: each value holds from its own time until
# the next one.
step_frame <- function(df) {
  n <- nrow(df)
  if (n < 2) {
    return(df)
  }
  time <- c(df$time[1], rep(df$time[-1], each = 2))
  idx <- c(rep(seq_len(n - 1), each = 2), n)

  out <- data.frame(time = time, estimate = df$estimate[idx])
  if (!is.null(df$lower)) {
    out$lower <- df$lower[idx]
    out$upper <- df$upper[idx]
  }
  out
}
