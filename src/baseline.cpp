/// @file baseline.cpp
/// @brief Implementation of the baseline estimators and of CompAK.

#include "baseline.h"

#include "at_risk.h"
#include "utils.h"

namespace gcmrec {

BaselineResult est_lamb_surv(const RecurrentData& d, double s,
                             const arma::vec& diseff, double alpha,
                             const arma::vec& beta, const arma::vec& z,
                             int rho_type) {
  const int ndiseff = static_cast<int>(diseff.n_elem);
  BaselineResult out;
  out.lambda.zeros(ndiseff);
  out.delta_lambda.zeros(ndiseff);
  out.surv.ones(ndiseff);
  out.var_lambda.zeros(ndiseff);

  const arma::ivec ns = nsm(d, s);

  // diseff[0] is effective age 0 (initial record): Lambda(0) = 0.
  for (int i = 1; i < ndiseff; ++i) {
    // Delta N(s, t): events observed before s whose effective age is
    // exactly diseff[i] (the values come from the same doubles, so the
    // exact comparison is intentional).
    double delta_n = 0.0;
    for (int j = 0; j < d.nk(); ++j) {
      if (d.caltimes[j] <= s && d.effage[j] == diseff[i]) {
        delta_n += 1.0;
      }
    }

    const S0Result s0r =
        at_risk(d, ns, s, diseff[i], alpha, beta, z, rho_type);

    if (s0r.s0 > 0.0) {
      const double dellamb = std::min(delta_n / s0r.s0, 1.0);
      out.lambda[i] = out.lambda[i - 1] + dellamb;
      out.surv[i] = out.surv[i - 1] * (1.0 - dellamb);
      // Aalen variance of the cumulative hazard: sum dN / S0^2
      out.var_lambda[i] =
          out.var_lambda[i - 1] + delta_n / (s0r.s0 * s0r.s0);
    } else {
      out.lambda[i] = out.lambda[i - 1];
      out.surv[i] = out.surv[i - 1];
      out.var_lambda[i] = out.var_lambda[i - 1];
    }
  }

  for (int i = 1; i < ndiseff; ++i) {
    out.delta_lambda[i] = out.lambda[i] - out.lambda[i - 1];
  }

  return out;
}

AKResult comp_ak(const RecurrentData& d, double s, const arma::vec& diseff,
                 double alpha, const arma::vec& beta,
                 const arma::vec& delta_lambda, const arma::vec& z,
                 int rho_type) {
  const int ndiseff = static_cast<int>(diseff.n_elem);
  AKResult out;
  out.kk = nsm(d, s);
  out.a.zeros(d.n);
  out.b.zeros(d.n);

  // A_i = sum_t Y_i(s, t) * dLambda0(t) over the effective ages.
  for (int t = 0; t < ndiseff; ++t) {
    const S0Result s0r =
        at_risk(d, out.kk, s, diseff[t], alpha, beta, z, rho_type);
    out.a += s0r.y_subj * delta_lambda[t];
  }

  // B_i: additive term of the marginal log-likelihood of xi. It does not
  // depend on xi, so it does not affect its maximization; the behavior of
  // the current Fortran is transcribed (which stops after the first event
  // whose effective age is found in diseff).
  for (int i = 0; i < d.n; ++i) {
    const int f = d.first[i];
    const arma::vec covariate = d.cov.col(f + d.k[i] - 1);

    if (out.kk[i] <= 1) continue;

    bool done = false;
    for (int j = 2; j <= out.kk[i] && !done; ++j) {
      out.b[i] += std::log(rho(j - 2, alpha, rho_type)) +
                  std::log(psi(covariate, beta, d.offset[i]));
      for (int u = 0; u < ndiseff && !done; ++u) {
        if (d.effage[f + j - 1] == diseff[u]) {
          out.b[i] += std::log(delta_lambda[u]);
          done = true;
        }
      }
    }
  }

  return out;
}

}  // namespace gcmrec
