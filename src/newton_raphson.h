/// @file newton_raphson.h
/// @brief Maximization of the partial log-likelihood by Newton-Raphson.

#ifndef GCMREC_NEWTON_RAPHSON_H
#define GCMREC_NEWTON_RAPHSON_H

#include "gcmrec_data.h"

namespace gcmrec {

/// Result of a Newton-Raphson maximization.
struct NRResult {
  arma::vec estim;    ///< final estimates
  arma::mat info_inv; ///< inverse of the information matrix (last iter.)
  double loglik = 0.0;
  int search = 0;     ///< 1 = convergence reached
  int kiter = 1;      ///< iterations used
};

/// Joint Newton-Raphson over (alpha, beta).
///
/// Replicates the strategy of the original Fortran: first alternates
/// partial maximizations over alpha and beta with a relaxed tolerance
/// (100 * tol) and, if they converge, finishes with a joint Newton-Raphson
/// at tolerance tol. `estim` has dimension nvar + 1 (alpha first).
NRResult newton_raphson(const RecurrentData& d, double s, double alpha_seed,
                        const arma::vec& beta_seed, const arma::vec& z,
                        int rho_type, double tol, int maxiter);

/// Newton-Raphson over beta only (rho = identity case).
///
/// `estim` has dimension nvar; alpha stays fixed at alpha_seed.
NRResult newton_raphson_beta(const RecurrentData& d, double s,
                             double alpha_seed, const arma::vec& beta_seed,
                             const arma::vec& z, int rho_type, double tol,
                             int maxiter);

}  // namespace gcmrec

#endif  // GCMREC_NEWTON_RAPHSON_H
