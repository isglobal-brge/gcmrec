# ------------------------------------------------------------------------
# build-reference-fits.R
#
# Generates the reference values ("gold standard") for the regression
# tests of the Fortran -> C++ port.
#
# Fits the datasets shipped with the package using the ORIGINAL Fortran
# implementation (v1.0-5) in every relevant configuration and saves the
# results to tests/testthat/fixtures/fortran-reference.rds.
#
# You only need this script to REGENERATE the reference values (for
# instance after adding a new case). Running the test suite does not
# require it, nor the Fortran code: the fixture already holds the values.
#
# Usage:
#   1. Get the original Fortran sources. They are not kept in this
#      repository to avoid duplication: they live in the history, under
#      the tag v1.0-5, which can be downloaded from
#      https://github.com/isglobal-brge/gcmrec/releases/tag/v1.0-5
#      (they are also in the CRAN archive).
#   2. Install that version into an auxiliary library. Current gfortran
#      needs the legacy flags:
#        printf 'FFLAGS = -O2 -std=legacy -fallow-argument-mismatch\n' > /tmp/Makevars-legacy
#        R_MAKEVARS_USER=/tmp/Makevars-legacy R CMD INSTALL --library=<lib> .
#   3. Run this script against that library:
#        Rscript tools/build-reference-fits.R <lib>
#
# Note: the typeEffage = "minimal" configuration uses rbinom() internally,
# so the seed is fixed; the port must reproduce it with the same seed
# and the same RNG.
# ------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
lib <- if (length(args) >= 1) args[[1]] else .libPaths()[1]

library(gcmrec, lib.loc = lib)

data(readmission)
data(lymphoma)

# Extract from the fitted object only what is numerically comparable
extract_fit <- function(fit) {
  fit[c("coef", "var", "loglik", "n", "nk", "kiter", "search",
        "Xi", "frailties", "diseff", "Lambda", "Survival",
        "rho.type", "se.type")]
}

reference <- list()
timings <- list()

run_case <- function(name, expr) {
  cat("Ajustando:", name, "... ")
  t0 <- proc.time()
  fit <- eval(expr)
  elapsed <- (proc.time() - t0)[["elapsed"]]
  cat(round(elapsed, 1), "s\n")
  reference[[name]] <<- extract_fit(fit)
  timings[[name]] <<- elapsed
  invisible(fit)
}

## 1. readmission, rho = alpha^k, information matrix SEs -------------------
run_case("readmission_alphak", quote(
  gcmrec(Survr(id, time, event) ~ as.factor(dukes),
         data = readmission, s = 3000)
))

## 2. readmission, rho = identity -----------------------------------------
run_case("readmission_identity", quote(
  gcmrec(Survr(id, time, event) ~ as.factor(dukes),
         data = readmission, s = 3000, rhoFunc = "Identity")
))

## 3. readmission, several covariates -------------------------------------
run_case("readmission_multicov", quote(
  gcmrec(Survr(id, time, event) ~ chemo + sex + dukes + charlson,
         data = readmission, s = 3000)
))

## 4. readmission, minimal repair (uses rbinom -> fixed seed) -------------
set.seed(20260804)
run_case("readmission_minimal", quote(
  gcmrec(Survr(id, time, event) ~ as.factor(dukes),
         data = readmission, s = 3000, typeEffage = "minimal")
))

## 5. readmission, frailty model (EM) --------------------------------------
run_case("readmission_frailty", quote(
  gcmrec(Survr(id, time, event) ~ as.factor(dukes),
         data = readmission, s = 3000, Frailty = TRUE)
))

## 6. lymphoma, cancer model (CR/PR/SD effective age) ----------------------
run_case("lymphoma_cancer", quote(
  gcmrec(Survr(id, time, event) ~ as.factor(distrib),
         data = lymphoma, s = 1000, cancer = lymphoma$effage)
))

## 7. jackknife on a subset (60 subjects, for computational cost) ----------
readm60 <- readmission[readmission$id %in% unique(readmission$id)[1:60], ]
run_case("readmission60_jackknife", quote(
  gcmrec(Survr(id, time, event) ~ as.factor(dukes),
         data = readm60, s = 3000, se = "Jacknife")
))

## Metadata -----------------------------------------------------------------
attr(reference, "meta") <- list(
  package_version = as.character(packageVersion("gcmrec")),
  r_version = R.version.string,
  rng_kind = RNGkind(),
  seed_minimal = 20260804,
  created = format(Sys.time(), tz = "UTC"),
  timings = timings
)

out <- file.path("tests", "testthat", "fixtures", "fortran-reference.rds")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
saveRDS(reference, out, version = 2)
cat("\nReferencia guardada en", out, "\n")
cat("Casos:", paste(names(reference), collapse = ", "), "\n")
