/// @file gcmrec_data.h
/// @brief Recurrent event data structure shared by the entire numerical
///   core of the gcmrec package.
///
/// C++ port (2026) of the original Fortran core by Juan R. Gonzalez (2003),
/// based on the general model of Pena and Hollander (2004).

#ifndef GCMREC_DATA_H
#define GCMREC_DATA_H

#include <RcppArmadillo.h>

namespace gcmrec {

/// Recurrent event data in "flat" format indexed by subject.
///
/// The flat vectors concatenate the records of all subjects, ordered by
/// subject (the record at time 0 included). The j-th (0-based) record of
/// subject i lives at position `first[i] + j`.
struct RecurrentData {
  int n;                 ///< number of subjects
  int nvar;              ///< number of covariates
  arma::ivec k;          ///< records per subject (includes time 0)
  arma::vec tau;         ///< administrative time per subject
  arma::vec censored;    ///< censoring time per subject
  arma::vec offset;      ///< offset per subject
  arma::vec caltimes;    ///< calendar times S_ij (flat)
  arma::vec gaptimes;    ///< inter-event gap times T_ij (flat)
  arma::vec intercepts;  ///< intercepts of the effective age function
  arma::vec slopes;      ///< slopes of the effective age function
  arma::vec lastperrep;  ///< index (1-based) of the last perfect repair
  arma::vec effagebegin; ///< effective age at the start of each period
  arma::vec effage;      ///< effective age at the end of each period
  arma::mat cov;         ///< covariates, dimension nvar x nk
  arma::ivec first;      ///< starting position (0-based) of each subject

  /// Total number of records (sum of k).
  int nk() const { return static_cast<int>(caltimes.n_elem); }

  /// Recomputes `first` from `k`. Call after modifying `k`.
  void index() {
    first.set_size(n);
    int pos = 0;
    for (int i = 0; i < n; ++i) {
      first[i] = pos;
      pos += k[i];
    }
  }

  /// Copy of the data excluding subject `drop` (for jackknife).
  RecurrentData without_subject(int drop) const {
    RecurrentData out;
    out.n = n - 1;
    out.nvar = nvar;

    const int lo = first[drop];              // subject's first record
    const int hi = lo + k[drop] - 1;         // subject's last record

    out.k = drop_elem(k, drop);
    out.tau = drop_elem(tau, drop);
    out.censored = drop_elem(censored, drop);
    out.offset = drop_elem(offset, drop);

    out.caltimes = drop_range(caltimes, lo, hi);
    out.gaptimes = drop_range(gaptimes, lo, hi);
    out.intercepts = drop_range(intercepts, lo, hi);
    out.slopes = drop_range(slopes, lo, hi);
    out.lastperrep = drop_range(lastperrep, lo, hi);
    out.effagebegin = drop_range(effagebegin, lo, hi);
    out.effage = drop_range(effage, lo, hi);

    out.cov = cov;
    out.cov.shed_cols(lo, hi);

    out.index();
    return out;
  }

 private:
  template <typename V>
  static V drop_elem(const V& x, int pos) {
    V out = x;
    out.shed_row(pos);
    return out;
  }

  static arma::vec drop_range(const arma::vec& x, int lo, int hi) {
    arma::vec out = x;
    out.shed_rows(lo, hi);
    return out;
  }
};

}  // namespace gcmrec

#endif  // GCMREC_DATA_H
