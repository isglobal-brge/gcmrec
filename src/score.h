/// @file score.h
/// @brief Partial log-likelihood, score vector and information matrix.

#ifndef GCMREC_SCORE_H
#define GCMREC_SCORE_H

#include "gcmrec_data.h"

namespace gcmrec {

/// Result of a score evaluation at (alpha, beta).
///
/// `score` and `info` are the joint (alpha, beta) versions of dimension
/// nvar + 1; the separate blocks are used in the partial maximizations.
struct ScoreResult {
  double loglik = 0.0;
  double score_alpha = 0.0;
  double info_alpha = 0.0;
  arma::vec score_beta;
  arma::vec info_alpha_beta;
  arma::mat info_beta;
  arma::vec score;  ///< (score_alpha, score_beta)
  arma::mat info;   ///< joint information matrix

  explicit ScoreResult(int nvar)
      : score_beta(nvar, arma::fill::zeros),
        info_alpha_beta(nvar, arma::fill::zeros),
        info_beta(nvar, nvar, arma::fill::zeros),
        score(nvar + 1, arma::fill::zeros),
        info(nvar + 1, nvar + 1, arma::fill::zeros) {}
};

/// Evaluates log-likelihood, score and information matrix at (alpha, beta).
///
/// @param d        recurrent event data
/// @param s        selected calendar time
/// @param alpha    alpha parameter
/// @param beta     covariate coefficients
/// @param z        frailties (ones if there are no frailties)
/// @param rho_type type of rho function
ScoreResult score_func(const RecurrentData& d, double s, double alpha,
                       const arma::vec& beta, const arma::vec& z,
                       int rho_type);

}  // namespace gcmrec

#endif  // GCMREC_SCORE_H
