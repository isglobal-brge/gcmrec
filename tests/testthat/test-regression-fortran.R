# Numerical regression tests against the original Fortran implementation.
#
# Each test refits the model and compares coefficients, log-likelihood,
# variances and baseline functions with the reference values generated
# by tools/build-reference-fits.R. See helper-reference.R.
#
# Note: gcmrec() emits a known model.matrix warning ("non-list
# contrasts argument ignored") inherited from the original version; it
# is suppressed explicitly here until the refactor fixes it.

fit_quiet <- function(...) {
  suppressWarnings(gcmrec(...))
}

test_that("readmission, rho = alpha^k, information matrix SEs", {
  data(readmission)
  fit <- fit_quiet(Survr(id, time, event) ~ as.factor(dukes),
                   data = readmission, s = 3000)
  expect_matches_reference(fit, "readmission_alphak")
})

test_that("readmission, rho = identity", {
  data(readmission)
  fit <- fit_quiet(Survr(id, time, event) ~ as.factor(dukes),
                   data = readmission, s = 3000, rhoFunc = "Identity")
  # Lambda/Survival are not compared with the reference: the original
  # Fortran computed them with an incomplete beta (out-of-bounds read,
  # fixed in the port). Structural properties are checked instead.
  expect_matches_reference(fit, "readmission_identity",
                           check_baseline = FALSE)
  expect_true(all(diff(fit$Lambda) >= 0))
  expect_true(all(fit$Survival >= 0 & fit$Survival <= 1))
  expect_true(all(diff(fit$Survival) <= 0))
})

test_that("readmission, several covariates", {
  data(readmission)
  fit <- fit_quiet(Survr(id, time, event) ~ chemo + sex + dukes + charlson,
                   data = readmission, s = 3000)
  expect_matches_reference(fit, "readmission_multicov")
})

test_that("readmission, minimal repair (typeEffage = 'minimal')", {
  data(readmission)
  set.seed(SEED_MINIMAL)
  fit <- fit_quiet(Survr(id, time, event) ~ as.factor(dukes),
                   data = readmission, s = 3000, typeEffage = "minimal")
  expect_matches_reference(fit, "readmission_minimal")
})

test_that("readmission, frailty model (EM)", {
  skip_on_cran()  # ~15 s
  data(readmission)
  fit <- fit_quiet(Survr(id, time, event) ~ as.factor(dukes),
                   data = readmission, s = 3000, Frailty = TRUE)
  expect_matches_reference(fit, "readmission_frailty")
})

test_that("lymphoma, cancer model (CR/PR/SD effective age)", {
  data(lymphoma)
  fit <- fit_quiet(Survr(id, time, event) ~ as.factor(distrib),
                   data = lymphoma, s = 1000, cancer = lymphoma$effage)
  expect_matches_reference(fit, "lymphoma_cancer")
})

test_that("frailty with maxXi = 'Brent' agrees with Newton-Raphson", {
  skip_on_cran()  # ~30 s
  # The Brent branch is not covered by the fixtures (the original used
  # Numerical Recipes code, replaced by an in-house implementation):
  # it is validated by agreement with the Newton-Raphson maximization.
  data(readmission)
  fit_nr <- fit_quiet(Survr(id, time, event) ~ as.factor(dukes),
                      data = readmission, s = 3000, Frailty = TRUE)
  fit_br <- fit_quiet(Survr(id, time, event) ~ as.factor(dukes),
                      data = readmission, s = 3000, Frailty = TRUE,
                      maxXi = "Brent")
  expect_equal(fit_br$coef, fit_nr$coef, tolerance = 1e-3)
  expect_equal(fit_br$Xi / (1 + fit_br$Xi), fit_nr$Xi / (1 + fit_nr$Xi),
               tolerance = 1e-3)
})

test_that("readmission (60 subjects), jackknife SEs", {
  skip_on_cran()
  data(readmission)
  readm60 <- readmission[readmission$id %in% unique(readmission$id)[1:60], ]
  fit <- fit_quiet(Survr(id, time, event) ~ as.factor(dukes),
                   data = readm60, s = 3000, se = "Jacknife")
  expect_matches_reference(fit, "readmission60_jackknife")
})
