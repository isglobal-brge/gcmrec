/// @file frailty.h
/// @brief Gamma frailty model: E step, marginal likelihood of xi and the
///   full EM algorithm.

#ifndef GCMREC_FRAILTY_H
#define GCMREC_FRAILTY_H

#include "gcmrec_data.h"

namespace gcmrec {

/// Methods for maximizing the marginal likelihood of xi.
enum XiMethod { XI_NEWTON_RAPHSON = 1, XI_BRENT = 2 };

/// Sentinel value for xi -> infinity (frailty variance -> 0).
constexpr double XI_HUGE = 1.0e60;

/// Result of the EM algorithm with frailties.
struct EMResult {
  arma::vec estim;  ///< (alpha, beta, xi, Z_1..Z_n)
  int search = 0;   ///< 1 = convergence reached
  int kiter = 0;    ///< EM iterations used
  double loglik = 0.0;
};

/// Marginal log-likelihood of xi and its first two derivatives.
///
/// @param kk number of events per subject (already without the time-0 record)
/// @param a  A_i computed by comp_ak
/// @param b  B_i computed by comp_ak (additive constant in xi)
void loglik_xi(double xi, const arma::ivec& kk, const arma::vec& a,
               const arma::vec& b, double& g0, double& g1, double& g2);

/// E step: conditional expectation of the frailties,
/// Z_i = (xi + K_i) / (xi + A_i).
///
/// @param k_raw unadjusted K_i = n_i(s-) (the time-0 record is discounted here)
arma::vec frailty_values(double xi, const arma::ivec& k_raw,
                         const arma::vec& a);

/// Maximizes the marginal likelihood with respect to xi by Newton-Raphson
/// on eta = log(xi) (original MaxWrtXi subroutine).
double max_wrt_xi(double xi_old, const arma::ivec& k_raw, const arma::vec& a,
                  const arma::vec& b, double tol, int maxiter, int& search);

/// Maximizes the marginal likelihood with respect to xi with Brent's method.
///
/// The original used the mnbrak/brent routines from Numerical Recipes
/// (license incompatible with the GPL); here Brent_fmin from the R API is
/// used over the interval [1e-10, 1e8], so there may be minimal numerical
/// differences from version 1.0-5 for this option.
double max_xi_brent(const arma::ivec& k_raw, const arma::vec& a,
                    const arma::vec& b);

/// Full EM algorithm for the gamma frailty model
/// (original EstimWithFrailty subroutine).
EMResult estim_with_frailty(const RecurrentData& d, double s,
                            const arma::vec& diseff, double alpha_seed,
                            const arma::vec& beta_seed, double xi_seed,
                            const arma::vec& z_seed, int rho_type, double tol,
                            int maxiter, int xi_method);

}  // namespace gcmrec

#endif  // GCMREC_FRAILTY_H
