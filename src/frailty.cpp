/// @file frailty.cpp
/// @brief Implementation of the gamma frailty model (EM).

#include "frailty.h"

#include <cmath>
#include <limits>

#include "baseline.h"
#include "newton_raphson.h"
#include "utils.h"

namespace gcmrec {

void loglik_xi(double xi, const arma::ivec& kk, const arma::vec& a,
               const arma::vec& b, double& g0, double& g1, double& g2) {
  const int n = static_cast<int>(kk.n_elem);
  g0 = 0.0;
  g1 = 0.0;
  g2 = 0.0;

  for (int i = 0; i < n; ++i) {
    double q1 = 0.0, q2 = 0.0, q3 = 0.0;
    for (int j = 1; j <= kk[i]; ++j) {
      const double denom = xi + (j - 1);
      q1 += std::log(denom);
      q2 += 1.0 / denom;
      q3 += 1.0 / (denom * denom);
    }
    const double xa = xi + a[i];
    g0 += q1 + xi * std::log(xi) - (xi + kk[i]) * std::log(xa) + b[i];
    g1 += q2 + std::log(xi) + 1.0 - std::log(xa) - (xi + kk[i]) / xa;
    g2 += -q3 + 1.0 / xi - 1.0 / xa - (a[i] - kk[i]) / (xa * xa);
  }
}

arma::vec frailty_values(double xi, const arma::ivec& k_raw,
                         const arma::vec& a) {
  const int n = static_cast<int>(k_raw.n_elem);
  arma::vec z(n);

  if (xi == XI_HUGE) {
    z.ones();
    return z;
  }
  for (int i = 0; i < n; ++i) {
    z[i] = (xi + (k_raw[i] - 1)) / (xi + a[i]);
  }
  return z;
}

double max_wrt_xi(double xi_old, const arma::ivec& k_raw, const arma::vec& a,
                  const arma::vec& b, double tol, int maxiter, int& search) {
  const arma::ivec kk = k_raw - 1;  // discounts the record at time 0

  double xi = xi_old;
  double eta = std::log(xi);
  double xirat_old = xi / (1.0 + xi);
  double dist = 100.0;
  double xi_new = xi;
  int niter = 0;
  search = 0;

  while (dist > tol && niter < maxiter) {
    ++niter;
    double g0, g1, g2;
    loglik_xi(xi, kk, a, b, g0, g1, g2);

    // Newton-Raphson on eta = log(xi) to guarantee xi > 0.
    const double eta_new = eta - g1 / (g1 + xi * g2);
    xi_new = std::exp(eta_new);

    if (xi_new > XI_HUGE) {
      search = 1;
      return XI_HUGE;
    }
    const double xirat_new = xi_new / (1.0 + xi_new);
    dist = std::abs(xirat_old - xirat_new);
    xi = xi_new;
    eta = eta_new;
    xirat_old = xirat_new;
  }

  if (dist <= tol) {
    search = 1;
  }
  return xi_new;
}

namespace {

/// Brent's (1973) one-dimensional minimizer: combines golden-section
/// search with parabolic interpolation. Own implementation of the classic
/// algorithm (the original used Numerical Recipes code, whose license is
/// incompatible with the GPL).
///
/// @param a,b search interval
/// @param f   objective function to minimize
/// @param t   absolute tolerance
template <typename F>
double brent_fmin(double a, double b, F f, double t) {
  const double c = 0.5 * (3.0 - std::sqrt(5.0));  // complementary golden ratio
  const double eps = std::sqrt(std::numeric_limits<double>::epsilon());

  double x = a + c * (b - a);
  double w = x, v = x;
  double d = 0.0, e = 0.0;
  double fx = f(x), fw = fx, fv = fx;

  for (;;) {
    const double m = 0.5 * (a + b);
    const double tol = eps * std::abs(x) + t;
    const double t2 = 2.0 * tol;

    if (std::abs(x - m) <= t2 - 0.5 * (b - a)) break;

    // Parabolic fit through (v, w, x) if acceptable; otherwise, golden step.
    double p = 0.0, q = 0.0, r = 0.0;
    if (std::abs(e) > tol) {
      r = (x - w) * (fx - fv);
      q = (x - v) * (fx - fw);
      p = (x - v) * q - (x - w) * r;
      q = 2.0 * (q - r);
      if (q > 0.0) p = -p; else q = -q;
      r = e;
      e = d;
    }

    if (std::abs(p) < std::abs(0.5 * q * r) && p > q * (a - x) &&
        p < q * (b - x)) {
      d = p / q;
      const double u_try = x + d;
      if (u_try - a < t2 || b - u_try < t2) d = (x < m) ? tol : -tol;
    } else {
      e = (x < m) ? b - x : a - x;
      d = c * e;
    }

    const double u =
        (std::abs(d) >= tol) ? x + d : x + ((d > 0.0) ? tol : -tol);
    const double fu = f(u);

    if (fu <= fx) {
      if (u < x) b = x; else a = x;
      v = w; fv = fw;
      w = x; fw = fx;
      x = u; fx = fu;
    } else {
      if (u < x) a = u; else b = u;
      if (fu <= fw || w == x) {
        v = w; fv = fw;
        w = u; fw = fu;
      } else if (fu <= fv || v == x || v == w) {
        v = u; fv = fu;
      }
    }
  }

  return x;
}

}  // namespace

double max_xi_brent(const arma::ivec& k_raw, const arma::vec& a,
                    const arma::vec& b) {
  const arma::ivec kk = k_raw - 1;

  // The search is done on eta = log(xi) (like the Newton-Raphson in
  // max_wrt_xi): on the linear scale the marginal likelihood has a nearly
  // flat plateau towards xi -> infinity that confuses the golden-section
  // search; on the log scale the interior optimum is well resolved.
  const auto neg_marginal_loglik_eta = [&](double eta) {
    double g0, g1, g2;
    loglik_xi(std::exp(eta), kk, a, b, g0, g1, g2);
    return -g0;
  };
  const double eta_min = brent_fmin(std::log(1.0e-8), std::log(1.0e8),
                                    neg_marginal_loglik_eta, 1.0e-8);
  return std::exp(eta_min);
}

EMResult estim_with_frailty(const RecurrentData& d, double s,
                            const arma::vec& diseff, double alpha_seed,
                            const arma::vec& beta_seed, double xi_seed,
                            const arma::vec& z_seed, int rho_type, double tol,
                            int maxiter, int xi_method) {
  EMResult out;

  // Step 0: initialization with the seed values.
  arma::vec z_old = z_seed;
  double xi_old = xi_seed;
  double xirat_old = xi_old / (1.0 + xi_old);
  double alpha_old = alpha_seed;
  arma::vec beta_old = beta_seed;
  double loglik = 0.0;

  BaselineResult base_old =
      est_lamb_surv(d, s, diseff, alpha_old, beta_old, z_old, rho_type);
  AKResult ak = comp_ak(d, s, diseff, alpha_old, beta_old,
                        base_old.delta_lambda, z_old, rho_type);

  double dist_all = 10.0;
  int kiter = 0;

  while (dist_all > tol && kiter < maxiter) {
    ++kiter;

    // E step: new frailties given (alpha, beta, xi, Lambda).
    const arma::vec z_new = frailty_values(xi_old, ak.kk, ak.a);

    // M step 1: new Lambda0 with the updated frailties.
    const BaselineResult base_new =
        est_lamb_surv(d, s, diseff, alpha_old, beta_old, z_new, rho_type);

    // M step 2: new (alpha, beta) by Newton-Raphson.
    double alpha_new;
    arma::vec beta_new;
    if (rho_type == RHO_ALPHA_TO_K) {
      const NRResult nr = newton_raphson(d, s, alpha_old, beta_old, z_new,
                                         rho_type, tol, maxiter);
      alpha_new = nr.estim[0];
      beta_new = nr.estim.subvec(1, d.nvar);
      loglik = nr.loglik;
    } else {
      const NRResult nr = newton_raphson_beta(d, s, alpha_old, beta_old,
                                              z_new, rho_type, tol, maxiter);
      alpha_new = 1.0;
      beta_new = nr.estim;
      loglik = nr.loglik;
    }

    // M step 3: new xi by maximizing the marginal likelihood.
    ak = comp_ak(d, s, diseff, alpha_new, beta_new, base_new.delta_lambda,
                 z_new, rho_type);

    double xi_new;
    if (xi_method == XI_NEWTON_RAPHSON) {
      int search_xi = 0;
      xi_new = max_wrt_xi(xi_old, ak.kk, ak.a, ak.b, tol, maxiter, search_xi);
    } else {
      xi_new = max_xi_brent(ak.kk, ak.a, ak.b);
    }
    const double xirat_new =
        (xi_new == XI_HUGE) ? 1.0 : xi_new / (1.0 + xi_new);

    // Convergence: maximum distance between successive iterations
    // (xi is compared on the xi / (1 + xi) scale).
    const double dist_alpha = std::abs(alpha_old - alpha_new);
    const double dist_beta = arma::norm(beta_old - beta_new, 2);
    const double dist_xi = std::abs(xirat_old - xirat_new);
    const double dist_lamb = arma::norm(base_old.lambda - base_new.lambda, 2);
    dist_all =
        std::max(std::max(std::max(dist_lamb, dist_xi), dist_beta),
                 dist_alpha);

    // Update for the next iteration.
    alpha_old = alpha_new;
    beta_old = beta_new;
    xi_old = xi_new;
    xirat_old = xirat_new;
    z_old = z_new;
    base_old = base_new;
  }

  out.search = (dist_all <= tol) ? 1 : 0;
  out.kiter = kiter;
  out.loglik = loglik;
  out.estim.set_size(d.nvar + 2 + d.n);
  out.estim[0] = alpha_old;
  out.estim.subvec(1, d.nvar) = beta_old;
  out.estim[d.nvar + 1] = xi_old;
  out.estim.subvec(d.nvar + 2, d.nvar + 1 + d.n) = z_old;

  return out;
}

}  // namespace gcmrec
