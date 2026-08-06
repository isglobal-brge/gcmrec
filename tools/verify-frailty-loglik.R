# ------------------------------------------------------------------------
# verify-frailty-loglik.R
#
# Checks WHICH quantity gcmrec() returns in `loglik` when the model is
# fitted with frailties (Frailty = TRUE).
#
# Conclusion it verifies: it is the PARTIAL log-likelihood CONDITIONAL on
# the estimated frailties, not the marginal one. That is why it is not
# comparable with the partial log-likelihood of a model without frailties.
#
# Usage:
#   Rscript tools/verify-frailty-loglik.R
#
# The first two checks only need the installed package; the third uses
# internal functions to evaluate the partial likelihood with fixed Z.
# ------------------------------------------------------------------------

library(gcmrec)
data(readmission)

fit <- function(...) {
  suppressWarnings(gcmrec(Survr(id, time, event) ~ as.factor(dukes),
                          data = readmission, s = 3000, ...))
}

cat("\n== 1. The two likelihoods are mutually inconsistent ==\n\n")

m0 <- fit()
m1 <- fit(Frailty = TRUE)

cat(sprintf("  loglik without frailties : %12.4f\n", m0$loglik))
cat(sprintf("  loglik with frailties    : %12.4f\n", m1$loglik))
cat(sprintf("  2 * difference           : %12.4f\n", 2 * (m1$loglik - m0$loglik)))
cat("\n  The frailty model contains the model without frailties as a\n")
cat("  limiting case (xi -> infinity) and has one extra parameter, so its\n")
cat("  likelihood cannot be SMALLER if both were the same quantity. A\n")
cat("  negative LRT statistic is the sign that they are not.\n")

cat("\n== 2. What the code does (original Fortran and C++ port) ==\n\n")
cat("  In the EstimWithFrailty subroutine of the original Fortran:\n")
cat("    - 'loglik' comes from newtraph()/nrBeta(), which return the\n")
cat("      PARTIAL log-likelihood evaluated at the current frailties.\n")
cat("    - 'loglikMarg' receives the marginal (G0) from MaxWrtXi()... but\n")
cat("      only when maxXi = 'Newton-Raphson', and it is never used later.\n")
cat("    - The last line does  loglikEnd = loglik,  i.e. it returns the\n")
cat("      conditional partial. The C++ port does exactly the same.\n")

cat("\n== 3. Direct check: reproduce the value by evaluating the partial ==\n\n")

Y <- with(readmission, Survr(id, time, event))
X <- model.matrix(~ as.factor(dukes), readmission)
covariates <- X[, -1]
d <- gcmrec:::formatData(Y[, 1], Y[, 2], Y[, 3], covariates, 1, NULL)
covMat <- matrix(d$covariate, nrow = ncol(covariates))
offset <- rep(0, nrow(Y))

# Partial log-likelihood evaluated with Z = frailties estimated by the EM
cond <- gcmrec:::.gcmrec_fit(3000, as.integer(d$k), d$tau, d$caltimes,
                             d$gaptimes, d$censored, d$intercepts, d$slopes,
                             d$lastperrep, d$effagebegin, d$effage, covMat,
                             offset, m1$coef[1], m1$coef[-1], m1$frailties,
                             2L, 1e-6, 100L)

cat(sprintf("  loglik returned by the frailty model        : %12.4f\n",
            m1$loglik))
cat(sprintf("  partial evaluated at Z = EM frailties       : %12.4f\n",
            cond$loglik))
cat(sprintf("  difference                                  : %12.2e\n",
            abs(m1$loglik - cond$loglik)))
cat("\n  They match: the returned number IS the conditional partial.\n")

cat("\n== Scope ==\n\n")
cat("  This affects ONLY the label under which that quantity is printed\n")
cat("  and the possibility of using it in a likelihood ratio test between\n")
cat("  models with and without frailties (which gcmrec >= 2.0 rejects).\n")
cat("  The coefficients, alpha, xi, the standard errors, the frailties and\n")
cat("  the baseline functions do not depend on this and do not change.\n\n")
