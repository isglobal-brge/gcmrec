/// @file jackknife.cpp
/// @brief Implementation of the leave-one-out jackknife.
///
/// Each refit is independent, so the loop is parallelized with OpenMP.
/// The body makes no calls to the R API (only Armadillo and arithmetic),
/// a requirement for thread safety.

#include "jackknife.h"

#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "frailty.h"
#include "newton_raphson.h"

namespace gcmrec {

namespace {

/// Number of threads to use: `nthreads` if positive; otherwise whatever
/// OpenMP decides (which honors OMP_NUM_THREADS).
int resolve_threads(int nthreads) {
#ifdef _OPENMP
  const int available = omp_get_max_threads();
  return (nthreads > 0) ? std::min(nthreads, available) : available;
#else
  (void)nthreads;
  return 1;
#endif
}

}  // namespace

arma::mat jackknife(const RecurrentData& d, double s, double alpha_seed,
                    const arma::vec& beta_seed, int rho_type, double tol,
                    int maxiter, int nthreads) {
  arma::mat out(d.n, d.nvar + 1, arma::fill::zeros);
  const int threads = resolve_threads(nthreads);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(threads)
#endif
  for (int i = 0; i < d.n; ++i) {
    const RecurrentData di = d.without_subject(i);
    const arma::vec z(di.n, arma::fill::ones);
    const NRResult nr = newton_raphson(di, s, alpha_seed, beta_seed, z,
                                       rho_type, tol, maxiter);
    out.row(i) = nr.estim.t();
  }

  return out;
}

arma::mat jackknife_frailty(const RecurrentData& d, double s,
                            const arma::vec& diseff, double alpha_seed,
                            const arma::vec& beta_seed, double xi_seed,
                            int rho_type, double tol, int maxiter,
                            int xi_method, int nthreads) {
  arma::mat out(d.n, d.nvar + 2, arma::fill::zeros);
  const int threads = resolve_threads(nthreads);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(threads)
#endif
  for (int i = 0; i < d.n; ++i) {
    const RecurrentData di = d.without_subject(i);
    const arma::vec z(di.n, arma::fill::ones);
    const EMResult em =
        estim_with_frailty(di, s, diseff, alpha_seed, beta_seed, xi_seed, z,
                           rho_type, tol, maxiter, xi_method);
    out.row(i) = em.estim.subvec(0, d.nvar + 1).t();
  }

  return out;
}

}  // namespace gcmrec
