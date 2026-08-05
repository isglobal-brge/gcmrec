/// @file exports.cpp
/// @brief Rcpp interface: conversion of the data coming from R into the
///   internal structure and functions exported to the package.
///
/// All exported functions are internal to the package ("." prefix);
/// they are not part of the public API.

#include <RcppArmadillo.h>

#include "baseline.h"
#include "frailty.h"
#include "jackknife.h"
#include "newton_raphson.h"
#include "utils.h"

using namespace Rcpp;

namespace {

/// Builds the internal structure from the flat vectors produced by
/// formatData() in R.
gcmrec::RecurrentData make_data(const IntegerVector& k,
                                const NumericVector& tau,
                                const NumericVector& caltimes,
                                const NumericVector& gaptimes,
                                const NumericVector& censored,
                                const NumericVector& intercepts,
                                const NumericVector& slopes,
                                const NumericVector& lastperrep,
                                const NumericVector& effagebegin,
                                const NumericVector& effage,
                                const NumericMatrix& cov,
                                const NumericVector& offset) {
  gcmrec::RecurrentData d;
  d.n = k.size();
  d.nvar = cov.nrow();

  d.k.set_size(d.n);
  for (int i = 0; i < d.n; ++i) {
    d.k[i] = k[i];
  }

  d.tau = arma::vec(tau.begin(), tau.size());
  d.censored = arma::vec(censored.begin(), censored.size());
  d.caltimes = arma::vec(caltimes.begin(), caltimes.size());
  d.gaptimes = arma::vec(gaptimes.begin(), gaptimes.size());
  d.intercepts = arma::vec(intercepts.begin(), intercepts.size());
  d.slopes = arma::vec(slopes.begin(), slopes.size());
  d.lastperrep = arma::vec(lastperrep.begin(), lastperrep.size());
  d.effagebegin = arma::vec(effagebegin.begin(), effagebegin.size());
  d.effage = arma::vec(effage.begin(), effage.size());
  d.cov = arma::mat(cov.begin(), cov.nrow(), cov.ncol());

  // The historical offset arrives with length >= n; only the first n
  // elements are used (one per subject), as in the Fortran version.
  if (offset.size() < d.n) {
    stop("offset must have at least n elements");
  }
  d.offset = arma::vec(offset.begin(), offset.size()).head(d.n);

  d.index();
  return d;
}

}  // namespace

/// Fit of the model without frailties. With rho_type = 2 it maximizes over
/// (alpha, beta); with rho_type = 1 (identity) over beta only.
// [[Rcpp::export(.gcmrec_fit)]]
List gcmrec_fit(double s, IntegerVector k, NumericVector tau,
                NumericVector caltimes, NumericVector gaptimes,
                NumericVector censored, NumericVector intercepts,
                NumericVector slopes, NumericVector lastperrep,
                NumericVector effagebegin, NumericVector effage,
                NumericMatrix cov, NumericVector offset, double alpha_seed,
                NumericVector beta_seed, NumericVector z, int rho_type,
                double tol, int maxit) {
  const gcmrec::RecurrentData d =
      make_data(k, tau, caltimes, gaptimes, censored, intercepts, slopes,
                lastperrep, effagebegin, effage, cov, offset);
  const arma::vec beta0(beta_seed.begin(), beta_seed.size());
  const arma::vec zz(z.begin(), z.size());

  gcmrec::NRResult nr;
  if (rho_type == gcmrec::RHO_ALPHA_TO_K) {
    nr = gcmrec::newton_raphson(d, s, alpha_seed, beta0, zz, rho_type, tol,
                                maxit);
  } else {
    nr = gcmrec::newton_raphson_beta(d, s, alpha_seed, beta0, zz, rho_type,
                                     tol, maxit);
  }

  return List::create(Named("loglik") = nr.loglik,
                      Named("estim") = nr.estim,
                      Named("info") = nr.info_inv,
                      Named("search") = nr.search,
                      Named("kiter") = nr.kiter);
}

/// Baseline cumulative hazard function Lambda0 and baseline survival.
// [[Rcpp::export(.gcmrec_baseline)]]
List gcmrec_baseline(double s, IntegerVector k, NumericVector tau,
                     NumericVector caltimes, NumericVector gaptimes,
                     NumericVector censored, NumericVector intercepts,
                     NumericVector slopes, NumericVector lastperrep,
                     NumericVector effagebegin, NumericVector effage,
                     NumericMatrix cov, NumericVector offset,
                     NumericVector diseff, double alpha,
                     NumericVector beta, NumericVector z, int rho_type) {
  const gcmrec::RecurrentData d =
      make_data(k, tau, caltimes, gaptimes, censored, intercepts, slopes,
                lastperrep, effagebegin, effage, cov, offset);
  const arma::vec beta_(beta.begin(), beta.size());
  const arma::vec zz(z.begin(), z.size());
  const arma::vec diseff_(diseff.begin(), diseff.size());

  const gcmrec::BaselineResult base =
      gcmrec::est_lamb_surv(d, s, diseff_, alpha, beta_, zz, rho_type);

  return List::create(Named("Lambda") = base.lambda,
                      Named("DeltaLambda") = base.delta_lambda,
                      Named("Survival") = base.surv,
                      Named("varLambda") = base.var_lambda);
}

/// Fit of the model with gamma frailties via EM.
// [[Rcpp::export(.gcmrec_frailty_fit)]]
List gcmrec_frailty_fit(double s, IntegerVector k, NumericVector tau,
                        NumericVector caltimes, NumericVector gaptimes,
                        NumericVector censored, NumericVector intercepts,
                        NumericVector slopes, NumericVector lastperrep,
                        NumericVector effagebegin, NumericVector effage,
                        NumericMatrix cov, NumericVector offset,
                        NumericVector diseff, double alpha_seed,
                        NumericVector beta_seed, double xi_seed,
                        NumericVector z, int rho_type, double tol, int maxit,
                        int xi_method) {
  const gcmrec::RecurrentData d =
      make_data(k, tau, caltimes, gaptimes, censored, intercepts, slopes,
                lastperrep, effagebegin, effage, cov, offset);
  const arma::vec beta0(beta_seed.begin(), beta_seed.size());
  const arma::vec zz(z.begin(), z.size());
  const arma::vec diseff_(diseff.begin(), diseff.size());

  const gcmrec::EMResult em =
      gcmrec::estim_with_frailty(d, s, diseff_, alpha_seed, beta0, xi_seed,
                                 zz, rho_type, tol, maxit, xi_method);

  return List::create(Named("estim") = em.estim,
                      Named("search") = em.search,
                      Named("kiter") = em.kiter,
                      Named("loglik") = em.loglik);
}

/// Leave-one-out jackknife without frailties (parallel with OpenMP).
// [[Rcpp::export(.gcmrec_jackknife)]]
NumericMatrix gcmrec_jackknife(double s, IntegerVector k, NumericVector tau,
                               NumericVector caltimes, NumericVector gaptimes,
                               NumericVector censored,
                               NumericVector intercepts, NumericVector slopes,
                               NumericVector lastperrep,
                               NumericVector effagebegin, NumericVector effage,
                               NumericMatrix cov, NumericVector offset,
                               double alpha_seed, NumericVector beta_seed,
                               int rho_type, double tol, int maxit,
                               int nthreads) {
  const gcmrec::RecurrentData d =
      make_data(k, tau, caltimes, gaptimes, censored, intercepts, slopes,
                lastperrep, effagebegin, effage, cov, offset);
  const arma::vec beta0(beta_seed.begin(), beta_seed.size());

  const arma::mat est =
      gcmrec::jackknife(d, s, alpha_seed, beta0, rho_type, tol, maxit,
                        nthreads);
  return wrap(est);
}

/// Leave-one-out jackknife with frailties (parallel with OpenMP).
// [[Rcpp::export(.gcmrec_jackknife_frailty)]]
NumericMatrix gcmrec_jackknife_frailty(
    double s, IntegerVector k, NumericVector tau, NumericVector caltimes,
    NumericVector gaptimes, NumericVector censored, NumericVector intercepts,
    NumericVector slopes, NumericVector lastperrep,
    NumericVector effagebegin, NumericVector effage, NumericMatrix cov,
    NumericVector offset, NumericVector diseff, double alpha_seed,
    NumericVector beta_seed, double xi_seed, int rho_type, double tol,
    int maxit, int xi_method, int nthreads) {
  const gcmrec::RecurrentData d =
      make_data(k, tau, caltimes, gaptimes, censored, intercepts, slopes,
                lastperrep, effagebegin, effage, cov, offset);
  const arma::vec beta0(beta_seed.begin(), beta_seed.size());
  const arma::vec diseff_(diseff.begin(), diseff.size());

  const arma::mat est =
      gcmrec::jackknife_frailty(d, s, diseff_, alpha_seed, beta0, xi_seed,
                                rho_type, tol, maxit, xi_method, nthreads);
  return wrap(est);
}
