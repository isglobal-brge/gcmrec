/// @file at_risk.h
/// @brief Generalized at-risk process S0(s, w) and its derivatives with
///   respect to alpha and beta (equations 23-27 of Pena and Hollander, 2004).

#ifndef GCMREC_AT_RISK_H
#define GCMREC_AT_RISK_H

#include "gcmrec_data.h"

namespace gcmrec {

/// Aggregate S0(s, w) and its first and second derivatives.
///
/// Notation: `gr_*` first derivative, `gr2_*` second derivative;
/// suffix `alpha` / `beta` indicates with respect to which parameter.
struct S0Result {
  double s0 = 0.0;
  double gr_alpha = 0.0;
  double gr2_alpha = 0.0;
  arma::vec gr_beta;
  arma::vec gr2_alpha_beta;
  arma::mat gr2_beta;
  arma::vec y_subj;  ///< individual at-risk process Y_i(s, w) of each subject

  explicit S0Result(int nvar, int n)
      : gr_beta(nvar, arma::fill::zeros),
        gr2_alpha_beta(nvar, arma::fill::zeros),
        gr2_beta(nvar, nvar, arma::fill::zeros),
        y_subj(n, arma::fill::zeros) {}
};

/// Contribution of subject i to the at-risk process, Y_i(s, w), and derivatives.
struct YSubjResult {
  double y = 0.0;
  double gr_alpha = 0.0;
  double gr2_alpha = 0.0;
  arma::vec gr_beta;
  arma::vec gr2_alpha_beta;
  arma::mat gr2_beta;

  explicit YSubjResult(int nvar)
      : gr_beta(nvar, arma::fill::zeros),
        gr2_alpha_beta(nvar, arma::fill::zeros),
        gr2_beta(nvar, nvar, arma::fill::zeros) {}
};

/// Y_i(s, w) for subject i (Q_ij and R_i terms, eqs. 23-24) with derivatives.
///
/// @param d        recurrent event data
/// @param i        0-based subject index
/// @param nsi      n_i(s-), subject's records prior to s
/// @param s        selected calendar time
/// @param w        effective age at which the process is evaluated
/// @param alpha    alpha parameter
/// @param beta     vector of covariate coefficients
/// @param rho_type type of rho function (see RhoType)
YSubjResult at_risk_subj(const RecurrentData& d, int i, int nsi, double s,
                         double w, double alpha, const arma::vec& beta,
                         int rho_type);

/// S0(s, w) = sum_i Z_i * Y_i(s, w) and derivatives (eqs. 25-27).
///
/// @param ns vector n_i(s-) precomputed with nsm()
/// @param z  frailties (or ones if the model has no frailties)
S0Result at_risk(const RecurrentData& d, const arma::ivec& ns, double s,
                 double w, double alpha, const arma::vec& beta,
                 const arma::vec& z, int rho_type);

}  // namespace gcmrec

#endif  // GCMREC_AT_RISK_H
