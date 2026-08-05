/// @file score.cpp
/// @brief Score implementation (faithful transcription of Fortran scorefunc).

#include "score.h"

#include "at_risk.h"
#include "utils.h"

namespace gcmrec {

ScoreResult score_func(const RecurrentData& d, double s, double alpha,
                       const arma::vec& beta, const arma::vec& z,
                       int rho_type) {
  ScoreResult out(d.nvar);
  const arma::ivec ns = nsm(d, s);

  for (int i = 0; i < d.n; ++i) {
    const int f = d.first[i];
    // The original Fortran uses the covariates of the subject's last
    // record for the psi term of the score (time-fixed covariates).
    const arma::vec covariate = d.cov.col(f + d.k[i] - 1);

    if (ns[i] <= 1) continue;

    // Contribution of each event observed before s (j = 2..n_i(s-)).
    for (int j = 2; j <= ns[i]; ++j) {
      const double w = d.effage[f + j - 1];
      const S0Result s0r = at_risk(d, ns, s, w, alpha, beta, z, rho_type);

      out.loglik += std::log(rho(j - 2, alpha, rho_type)) +
                    std::log(psi(covariate, beta, d.offset[i])) -
                    std::log(s0r.s0);

      const double e_alpha = s0r.gr_alpha / s0r.s0;
      out.score_alpha += (j - 2) / alpha - e_alpha;
      out.info_alpha += (j - 2) / (alpha * alpha) + s0r.gr2_alpha / s0r.s0 -
                        e_alpha * e_alpha;

      const arma::vec e_beta = s0r.gr_beta / s0r.s0;
      out.score_beta += covariate - e_beta;
      out.info_alpha_beta +=
          s0r.gr2_alpha_beta / s0r.s0 - e_alpha * e_beta;
      out.info_beta += s0r.gr2_beta / s0r.s0 - e_beta * e_beta.t();
    }
  }

  // Assembly of the joint (alpha, beta) versions.
  out.score[0] = out.score_alpha;
  out.info(0, 0) = out.info_alpha;
  for (int t = 0; t < d.nvar; ++t) {
    out.score[t + 1] = out.score_beta[t];
    out.info(0, t + 1) = out.info_alpha_beta[t];
    out.info(t + 1, 0) = out.info_alpha_beta[t];
  }
  out.info.submat(1, 1, d.nvar, d.nvar) = out.info_beta;

  return out;
}

}  // namespace gcmrec
