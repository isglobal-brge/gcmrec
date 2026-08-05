/// @file baseline.h
/// @brief Estimators of the baseline cumulative hazard function Lambda0 and
///   the baseline survival, plus computation of A_i / K_i for the EM algorithm.

#ifndef GCMREC_BASELINE_H
#define GCMREC_BASELINE_H

#include "gcmrec_data.h"

namespace gcmrec {

/// Baseline estimators evaluated at the distinct effective ages (diseff).
struct BaselineResult {
  arma::vec lambda;        ///< cumulative Lambda0
  arma::vec delta_lambda;  ///< increments of Lambda0
  arma::vec surv;          ///< baseline survival (product-limit)
  arma::vec var_lambda;    ///< variance of Lambda0 (Aalen): sum dN / S0^2
};

/// Aalen-Breslow-type estimator of Lambda0 and baseline survival
/// (original EstLambSurv subroutine). Requires estimated (alpha, beta).
///
/// @param diseff distinct effective ages, sorted in ascending order
BaselineResult est_lamb_surv(const RecurrentData& d, double s,
                             const arma::vec& diseff, double alpha,
                             const arma::vec& beta, const arma::vec& z,
                             int rho_type);

/// Per-subject quantities for the E step of the EM with gamma frailties
/// (original CompAK subroutine).
struct AKResult {
  arma::ivec kk;  ///< K_i = n_i(s-) (includes the record at time 0)
  arma::vec a;    ///< A_i = integral of the at-risk process against Lambda0
  arma::vec b;    ///< B_i = additive term of the marginal log-likelihood
};

AKResult comp_ak(const RecurrentData& d, double s, const arma::vec& diseff,
                 double alpha, const arma::vec& beta,
                 const arma::vec& delta_lambda, const arma::vec& z,
                 int rho_type);

}  // namespace gcmrec

#endif  // GCMREC_BASELINE_H
