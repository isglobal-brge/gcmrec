# Tests for the S3 methods of the gcmrec object and for data normalization.

fit_readmission <- local({
  data(readmission, envir = environment())
  suppressWarnings(gcmrec(Survr(id, time, event) ~ as.factor(dukes),
                          data = readmission, s = 3000))
})

test_that("coef, vcov and logLik return what is expected", {
  expect_equal(coef(fit_readmission), fit_readmission$coef)
  expect_equal(vcov(fit_readmission), fit_readmission$var)

  # vcov is labeled with the coefficient names
  expect_equal(rownames(vcov(fit_readmission)),
               names(coef(fit_readmission)))
  expect_equal(colnames(vcov(fit_readmission)),
               names(coef(fit_readmission)))

  ll <- logLik(fit_readmission)
  expect_s3_class(ll, "logLik")
  expect_equal(as.numeric(ll), fit_readmission$loglik)
  expect_equal(attr(ll, "df"), length(fit_readmission$coef))
})

test_that("summary returns a printable summary.gcmrec object", {
  s <- summary(fit_readmission)
  expect_s3_class(s, "summary.gcmrec")
  # With rho = alpha^k the hazard ratio table excludes alpha
  expect_equal(nrow(s$conf.int), length(fit_readmission$coef) - 1)
  expect_true(all(s$conf.int[, "lower"] <= s$conf.int[, "exp(coef)"]))
  expect_true(all(s$conf.int[, "upper"] >= s$conf.int[, "exp(coef)"]))
  expect_output(print(s), "hr")
  # The print(fit) output still works
  expect_output(print(fit_readmission), "alpha")
})

test_that("plot returns a ggplot and lines is still available (deprecated)", {
  # Since version 2.0 plot() returns a ggplot object
  expect_s3_class(plot(fit_readmission), "ggplot")
  expect_s3_class(plot(fit_readmission, type.plot = "hazard"), "ggplot")
  expect_error(plot(fit_readmission, type.plot = "nope"),
               "hazard or survival")

  # lines() still works on base graphics, warning about its status
  pdf(NULL)
  on.exit(dev.off())
  graphics::plot(fit_readmission$diseff, fit_readmission$Survival,
                 type = "n")
  expect_warning(lines(fit_readmission, lty = 2), "deprecated")
})

test_that("as_gcmrec_data normalizes the input formats", {
  data(hydraulic)
  df <- as_gcmrec_data(hydraulic)
  expect_s3_class(df, "data.frame")
  expect_named(df, c("id", "time", "event", "covar.1", "covar.2"))

  # A data.frame passes through unchanged; other types fail with a clear message
  expect_identical(as_gcmrec_data(df), df)
  expect_error(as_gcmrec_data(1:10), "data.frame")
})

test_that("gcmrec accepts the legacy list format directly", {
  data(hydraulic)
  fit <- suppressWarnings(gcmrec(Survr(id, time, event) ~ covar.1 + covar.2,
                                 data = hydraulic, s = 4753))
  fit.df <- suppressWarnings(
    gcmrec(Survr(id, time, event) ~ covar.1 + covar.2,
           data = List.to.Dataframe(hydraulic), s = 4753)
  )
  expect_equal(fit$coef, fit.df$coef)
  expect_equal(fit$loglik, fit.df$loglik)
})

test_that("external effageData (GeneratedData) fits and converges", {
  data(GeneratedData)
  temp <- List.to.Dataframe(GeneratedData)
  fit <- suppressWarnings(
    gcmrec(Survr(id, time, event) ~ covar.1 + covar.2, data = temp,
           effageData = GeneratedData, s = 100)
  )
  expect_equal(fit$search, 1)
  expect_length(fit$coef, 3)  # alpha + 2 covariates
  expect_true(is.finite(fit$loglik))
})
