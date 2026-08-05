# gcmrec 2.0.0

Major rewrite of the package. The statistical model and the user interface
are unchanged, so existing analysis scripts keep working, but the numerical
core, the graphics and the documentation have been rebuilt, and several
tools for interpreting a fit have been added.

## Numerical core

* The Fortran 77 core has been rewritten in C++ using `Rcpp` and
  `RcppArmadillo`. Results are identical to version 1.0-5 (verified by a
  regression test suite against the original implementation), and fitting
  is roughly 6 times faster for the plain model, 7 times for the frailty
  model and up to 45 times for the jackknife.
* The limit of 200 recurrences per subject, imposed by the static arrays
  of the Fortran code, is gone.
* Jackknife standard errors are computed in parallel with OpenMP where
  available. The number of threads honours `OMP_NUM_THREADS` and can be
  capped with the new `gcmrecThreads()`.
* Third-party numerical routines with licences incompatible with the GPL
  have been removed. The one-dimensional maximiser used by
  `maxXi = "Brent"` is now an implementation of Brent's algorithm working
  on the log scale, which is more stable than the previous version.

## New features

* `mcf()` estimates the mean cumulative function of recurrent events, with
  pointwise confidence bands and optional grouping, plus `plot()` and
  `print()` methods.
* `predict()` returns the linear predictor, the relative risk or the
  survivor and cumulative hazard curves of given covariate profiles;
  `plotPredict()` draws them.
* `anova()` performs likelihood ratio tests between nested models, and in
  particular tests whether `alpha = 1` (no effect of accumulating
  occurrences). It refuses comparisons that are not valid likelihood ratio
  tests; see `?anova.gcmrec`.
* `plotForest()` draws the hazard ratios with their confidence intervals.
* `plotBaseline()` compares the baseline functions of several fits on the
  same axes.
* `coef()`, `vcov()` and `logLik()` methods, so `AIC()` and similar
  functions now work on a fit.
* `as_gcmrec_data()` is a new S3 generic that normalises the input data.
  As a result, `gcmrec()` and `graph.caltimes()` accept the legacy list
  format directly, without calling `List.to.Dataframe()` first, and
  support for further input formats only requires a new method.
* `theme_gcmrec()` gives the package graphics a consistent look and can be
  applied to a user's own plots.
* `gcmrec()` gains `keep.data`, which stores the prepared data in the
  fitted object so that `anova()` can refit the model.

## Graphics

* All plotting functions return `ggplot` objects, which can be restyled
  with the usual `ggplot2` syntax.
* `plot()` on a fit now draws a pointwise confidence band around the
  baseline function.
* `graph.caltimes()` has been redrawn as an event chart, with a follow-up
  segment per subject and a marker for the end of follow-up. Its
  `sortevents` argument, which previous versions accepted but ignored, now
  orders subjects by follow-up (default) or by number of events.
* `lines()` on a fit is deprecated, since `plot()` no longer draws on a
  base graphics device; use `plotBaseline()` with a list of models.

## Bug fixes

* The baseline functions estimated with `rhoFunc = "Identity"` were
  computed from an incomplete coefficient vector, which read past the end
  of the vector in the Fortran code. They are now correct.
* The jackknife of the frailty model passed the full offset and frailty
  vectors to each leave-one-out refit, so the values were misaligned with
  the reduced data.
* `summary()` dropped the first covariate coefficient when
  `rhoFunc = "Identity"`, mistaking it for `alpha`.
* `plot()` and `lines()` relied on partial matching of `$Surv` and `$Lam`
  to reach the `Survival` and `Lambda` components of a fit.
* The log-likelihood of a frailty fit was labelled as marginal; it is the
  partial log-likelihood conditional on the estimated frailties, and is
  labelled as such.
* `gcmrec()` no longer emits a spurious warning from `model.matrix()`
  about the `contrasts` argument.

## Documentation and infrastructure

* All exported functions are documented with `roxygen2`, with runnable
  examples, and the R sources are split one file per function.
* New vignette walking through a complete analysis: `vignette("gcmrec")`.
* Test suite based on `testthat`, including regression tests that check
  the C++ core against results produced by the original Fortran code.
* Data sets are shipped in `.rda` format with lazy loading.
* Native routines are registered, and the package declares its imports
  explicitly.


# gcmrec 1.0-5

* Last version released on CRAN, in 2012, with the numerical core written
  in Fortran 77.
