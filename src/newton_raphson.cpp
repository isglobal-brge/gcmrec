/// @file newton_raphson.cpp
/// @brief Implementation of the Newton-Raphson procedures.
///
/// Faithful transcription of the Fortran subroutines newtraph,
/// newtraphAlpha, newtraphBeta, newtraphBoth and nrBeta: same iteration
/// order, same stopping rules and same counters, to reproduce the results
/// of version 1.0-5 exactly.

#include "newton_raphson.h"

#include "score.h"

namespace gcmrec {

namespace {

/// Inverse of the information matrix. The original Fortran used a pivoted
/// inversion (MATINV); here arma::inv is used with the pseudo-inverse as a
/// safeguard against numerically singular matrices.
arma::mat invert_info(const arma::mat& info) {
  arma::mat inv;
  if (!arma::inv(inv, info)) {
    inv = arma::pinv(info);
  }
  return inv;
}

/// Partial 1-D maximization over alpha with beta fixed (newtraphAlpha).
void nr_alpha_step(const RecurrentData& d, double s, double& alpha,
                   const arma::vec& beta, const arma::vec& z, int rho_type,
                   double tol, int maxiter, int& search) {
  double distance = 1000.0;
  int kiter = 1;
  search = 0;

  while (distance > tol && kiter < maxiter) {
    ++kiter;
    const ScoreResult sc = score_func(d, s, alpha, beta, z, rho_type);
    const double alpha_new = alpha + sc.score_alpha / sc.info_alpha;
    distance = std::abs(alpha_new - alpha);
    alpha = alpha_new;
    if (distance <= tol) {
      search = 1;
    }
  }
}

/// Partial maximization over beta with alpha fixed (newtraphBeta / nrBeta).
NRResult nr_beta_loop(const RecurrentData& d, double s, double alpha,
                      const arma::vec& beta_seed, const arma::vec& z,
                      int rho_type, double tol, int maxiter) {
  NRResult out;
  arma::vec beta = beta_seed;
  arma::mat info_inv(d.nvar, d.nvar, arma::fill::zeros);
  double distance = 1000.0;
  double loglik = 0.0;
  int kiter = 1;

  while (distance > tol && kiter < maxiter) {
    ++kiter;
    const ScoreResult sc = score_func(d, s, alpha, beta, z, rho_type);
    loglik = sc.loglik;
    info_inv = invert_info(sc.info_beta);
    const arma::vec beta_new = beta + info_inv * sc.score_beta;
    distance = arma::norm(beta_new - beta, 2);
    beta = beta_new;
  }

  out.estim = beta;
  out.info_inv = info_inv;
  out.loglik = loglik;
  out.search = (distance <= tol) ? 1 : 0;
  out.kiter = kiter;
  return out;
}

/// Joint Newton-Raphson over (alpha, beta) (newtraphBoth).
NRResult nr_both(const RecurrentData& d, double s, double alpha_seed,
                 const arma::vec& beta_seed, const arma::vec& z, int rho_type,
                 double tol, int maxiter) {
  NRResult out;
  arma::vec esti(d.nvar + 1);
  esti[0] = alpha_seed;
  esti.subvec(1, d.nvar) = beta_seed;

  arma::mat info_inv(d.nvar + 1, d.nvar + 1, arma::fill::zeros);
  double distance = 1000.0;
  double loglik = 0.0;
  int kiter = 1;
  int search = 0;

  while (distance > tol && kiter < maxiter) {
    ++kiter;
    const ScoreResult sc = score_func(d, s, esti[0], esti.subvec(1, d.nvar),
                                      z, rho_type);
    loglik = sc.loglik;
    info_inv = invert_info(sc.info);
    const arma::vec esti_new = esti + info_inv * sc.score;
    distance = arma::norm(esti_new - esti, 2);
    esti = esti_new;
    if (distance <= tol) {
      search = 1;
    }
  }

  out.estim = esti;
  out.info_inv = info_inv;
  out.loglik = loglik;
  out.search = search;
  out.kiter = kiter;
  return out;
}

}  // namespace

NRResult newton_raphson(const RecurrentData& d, double s, double alpha_seed,
                        const arma::vec& beta_seed, const arma::vec& z,
                        int rho_type, double tol, int maxiter) {
  // Phase 1: alternate alpha / beta steps with relaxed tolerance.
  const double tol1 = 100.0 * tol;
  double alpha = alpha_seed;
  arma::vec beta = beta_seed;
  double distance = 1000.0;
  int kkiter = 1;
  int search = 0;

  while (distance > tol1 && kkiter < maxiter) {
    ++kkiter;

    arma::vec esti_old(d.nvar + 1);
    esti_old[0] = alpha;
    esti_old.subvec(1, d.nvar) = beta;

    int search_alpha = 0;
    nr_alpha_step(d, s, alpha, beta, z, rho_type, tol1, maxiter,
                  search_alpha);
    if (search_alpha != 1) {
      search = 121;  // original code: problems in the alpha step
    }

    const NRResult beta_step =
        nr_beta_loop(d, s, alpha, beta, z, rho_type, tol1, maxiter);
    beta = beta_step.estim;
    if (beta_step.search != 1) {
      search = 122;  // original code: problems in the beta step
    }

    arma::vec esti_new(d.nvar + 1);
    esti_new[0] = alpha;
    esti_new.subvec(1, d.nvar) = beta;
    distance = arma::norm(esti_new - esti_old, 2);
  }

  if (distance <= tol1) {
    search = 1;
  }

  // Phase 2: joint Newton-Raphson at full tolerance.
  NRResult out;
  if (search == 1) {
    out = nr_both(d, s, alpha, beta, z, rho_type, tol, maxiter);
  } else {
    out.estim = arma::join_cols(arma::vec{alpha}, beta);
    out.info_inv.zeros(d.nvar + 1, d.nvar + 1);
    out.search = 0;
  }
  out.kiter = kkiter;  // the counter exposed to R is the phase-1 one
  return out;
}

NRResult newton_raphson_beta(const RecurrentData& d, double s,
                             double alpha_seed, const arma::vec& beta_seed,
                             const arma::vec& z, int rho_type, double tol,
                             int maxiter) {
  return nr_beta_loop(d, s, alpha_seed, beta_seed, z, rho_type, tol, maxiter);
}

}  // namespace gcmrec
