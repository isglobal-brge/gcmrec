/// @file at_risk.cpp
/// @brief Implementation of the generalized at-risk process (eqs. 23-27).
///
/// Faithful transcription of the original Fortran AtRiskSubj/AtRisk
/// subroutines: the accumulation order and conditions are preserved to
/// reproduce the results of version 1.0-5 exactly.

#include "at_risk.h"

#include "utils.h"

namespace gcmrec {

YSubjResult at_risk_subj(const RecurrentData& d, int i, int nsi, double s,
                         double w, double alpha, const arma::vec& beta,
                         int rho_type) {
  YSubjResult out(d.nvar);
  const int f = d.first[i];  // subject's first record in the flat vectors
  const int kk = nsi;        // subject's records prior to s

  // Q_ij terms (eq. 23): already completed inter-event periods.
  // Fortran j = 2..kk; here jj = j - 1 (0-based over the records).
  for (int jj = 1; jj < kk; ++jj) {
    if (w > d.effagebegin[f + jj - 1] && w <= d.effage[f + jj]) {
      const arma::vec covariate = d.cov.col(f + jj - 1);
      const double q = rho(jj - 1, alpha, rho_type) *
                       psi(covariate, beta, d.offset[i]) /
                       d.slopes[f + jj - 1];

      // Derivatives with respect to alpha; for rho = 1 the convention of
      // the original Fortran (gradient = q) is kept, which is irrelevant
      // because in that case alpha is not estimated.
      const double gq_a =
          (rho_type == RHO_IDENTITY) ? q : (jj - 1) / alpha * q;
      const double g2q_a = (rho_type == RHO_IDENTITY)
                               ? q
                               : static_cast<double>((jj - 1) * (jj - 2)) /
                                     (alpha * alpha) * q;

      out.y += q;
      out.gr_alpha += gq_a;
      out.gr2_alpha += g2q_a;
      out.gr_beta += covariate * q;
      out.gr2_alpha_beta += covariate * gq_a;
      out.gr2_beta += q * covariate * covariate.t();
    }
  }

  // R_i term (eq. 24): period in progress at time s.
  // Note: the literal Fortran formula is transcribed (with slopes = 1 it
  // is equivalent to the expression in the paper).
  const int last = f + kk - 1;  // subject's record kk (1-based)
  const double effage_at_s =
      d.intercepts[last] + d.slopes[last] * std::min(s, d.tau[i]) -
      d.caltimes[f + static_cast<int>(d.lastperrep[last]) - 1];

  if (w > d.effagebegin[last] && w <= effage_at_s) {
    const arma::vec covariate = d.cov.col(last);
    const double r = rho(kk - 1, alpha, rho_type) *
                     psi(covariate, beta, d.offset[i]) / d.slopes[last];

    const double gr_a = (rho_type == RHO_IDENTITY) ? r : (kk - 1) / alpha * r;
    const double g2r_a = (rho_type == RHO_IDENTITY)
                             ? r
                             : static_cast<double>((kk - 1) * (kk - 2)) /
                                   (alpha * alpha) * r;

    out.y += r;
    out.gr_alpha += gr_a;
    out.gr2_alpha += g2r_a;
    out.gr_beta += covariate * r;
    out.gr2_alpha_beta += covariate * gr_a;
    out.gr2_beta += r * covariate * covariate.t();
  }

  return out;
}

S0Result at_risk(const RecurrentData& d, const arma::ivec& ns, double s,
                 double w, double alpha, const arma::vec& beta,
                 const arma::vec& z, int rho_type) {
  S0Result out(d.nvar, d.n);

  for (int i = 0; i < d.n; ++i) {
    const YSubjResult yi = at_risk_subj(d, i, ns[i], s, w, alpha, beta,
                                        rho_type);
    out.y_subj[i] = yi.y;
    out.s0 += z[i] * yi.y;
    out.gr_alpha += z[i] * yi.gr_alpha;
    out.gr2_alpha += z[i] * yi.gr2_alpha;
    out.gr_beta += z[i] * yi.gr_beta;
    out.gr2_alpha_beta += z[i] * yi.gr2_alpha_beta;
    out.gr2_beta += z[i] * yi.gr2_beta;
  }

  return out;
}

}  // namespace gcmrec
