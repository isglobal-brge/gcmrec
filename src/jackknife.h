/// @file jackknife.h
/// @brief Jackknife (leave-one-out) standard errors by subject.

#ifndef GCMREC_JACKKNIFE_H
#define GCMREC_JACKKNIFE_H

#include "gcmrec_data.h"

namespace gcmrec {

/// Jackknife for the model without frailties: refits the model excluding
/// each subject in turn. Returns an n x (nvar + 1) matrix with
/// (alpha, beta) per row. Parallelized with OpenMP when available.
///
/// @param nthreads maximum number of threads; <= 0 lets OpenMP decide
///   (which honors the OMP_NUM_THREADS environment variable)
arma::mat jackknife(const RecurrentData& d, double s, double alpha_seed,
                    const arma::vec& beta_seed, int rho_type, double tol,
                    int maxiter, int nthreads);

/// Jackknife for the model with frailties. Returns n x (nvar + 2) with
/// (alpha, beta, xi) per row.
///
/// Note: the original Fortran passed the full offset vector to each
/// leave-one-out refit (misaligned); here the reduced vector is passed.
arma::mat jackknife_frailty(const RecurrentData& d, double s,
                            const arma::vec& diseff, double alpha_seed,
                            const arma::vec& beta_seed, double xi_seed,
                            int rho_type, double tol, int maxiter,
                            int xi_method, int nthreads);

}  // namespace gcmrec

#endif  // GCMREC_JACKKNIFE_H
