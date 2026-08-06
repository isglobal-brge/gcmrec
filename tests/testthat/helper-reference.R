# Shared helpers for the regression tests.
#
# The reference ("gold standard") was generated with the original Fortran
# implementation (v1.0-5) via tools/build-reference-fits.R. Any
# reimplementation of the numerical core must reproduce these values.

reference_fits <- readRDS(test_path("fixtures", "fortran-reference.rds"))

# Seed used for the typeEffage = "minimal" case (it uses rbinom internally)
SEED_MINIMAL <- attr(reference_fits, "meta")$seed_minimal

# Compares a fresh fit against the reference case `name`.
# Numerical components must agree to the given tolerance.
#
# check_baseline = FALSE skips the Lambda/Survival comparison: for
# rho = "Identity" the original Fortran version computed the baseline
# function with undefined behavior (the R layer sliced coef as if alpha
# were present and the Fortran read beta out of bounds), so those
# reference values are neither reproducible nor correct. The port fixes
# the bug; see the corresponding test.
expect_matches_reference <- function(fit, name, tolerance = 1e-6,
                                     check_baseline = TRUE) {
  ref <- reference_fits[[name]]

  expect_equal(fit$coef, ref$coef, tolerance = tolerance)
  expect_equal(fit$loglik, ref$loglik, tolerance = tolerance)
  # The Fortran reference carried no dimnames on the variance matrix;
  # the port does set them, so only the values are compared
  expect_equal(unname(fit$var), unname(ref$var), tolerance = tolerance)
  expect_equal(fit$n, ref$n)
  expect_equal(fit$nk, ref$nk)
  expect_equal(fit$search, ref$search)
  expect_equal(fit$diseff, ref$diseff, tolerance = tolerance)
  if (check_baseline) {
    expect_equal(fit$Lambda, ref$Lambda, tolerance = tolerance)
    expect_equal(fit$Survival, ref$Survival, tolerance = tolerance)
  }

  if (!is.null(ref$Xi)) {
    expect_equal(fit$Xi, ref$Xi, tolerance = tolerance)
    expect_equal(fit$frailties, ref$frailties, tolerance = tolerance)
  }

  invisible(fit)
}
