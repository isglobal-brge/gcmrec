/// @file utils.h
/// @brief Basic model functions: rho, psi and the event count n_i(s-).

#ifndef GCMREC_UTILS_H
#define GCMREC_UTILS_H

#include "gcmrec_data.h"

namespace gcmrec {

/// Supported types of the rho(k; alpha) function (rhoFunc argument in R).
enum RhoType { RHO_IDENTITY = 1, RHO_ALPHA_TO_K = 2 };

/// Effect of the accumulated occurrences, rho(k; alpha).
///
/// @param j        number of previous occurrences
/// @param alpha    alpha parameter of the model
/// @param rho_type RHO_IDENTITY -> 1; RHO_ALPHA_TO_K -> alpha^j
inline double rho(int j, double alpha, int rho_type) {
  return (rho_type == RHO_IDENTITY) ? 1.0 : std::pow(alpha, j);
}

/// Covariate link function, psi(x) = exp(x'beta + offset).
inline double psi(const arma::vec& covariate, const arma::vec& beta,
                  double offset) {
  return std::exp(arma::dot(covariate, beta) + offset);
}

/// Number of the subject's records with calendar time prior to s,
/// n_i(s-). Includes the initial record at time 0.
///
/// @param caltimes subject's calendar times (first element 0)
/// @param s        selected calendar time
inline int nism(const arma::vec& caltimes, int k, double s) {
  if (s > caltimes[k - 1]) {
    return k;
  }
  int count = 0;
  for (int i = 0; i < k && caltimes[i] < s; ++i) {
    ++count;
  }
  return count;
}

/// n_i(s-) for all subjects.
inline arma::ivec nsm(const RecurrentData& d, double s) {
  arma::ivec ns(d.n);
  for (int i = 0; i < d.n; ++i) {
    ns[i] = nism(d.caltimes.subvec(d.first[i], d.first[i] + d.k[i] - 1),
                 d.k[i], s);
  }
  return ns;
}

}  // namespace gcmrec

#endif  // GCMREC_UTILS_H
