// todo (roll >= 1.2.1): use namespace roll for check_rows_equal

#ifndef ROLLOPTIM_CHECK_H
#define ROLLOPTIM_CHECK_H

#include <RcppArmadillo.h>
using namespace Rcpp;

namespace rolloptim {

// scalar checks: finite double
inline void check_finite(const double& value, const char* name) {

  if (std::isnan(value) || std::isinf(value)) {
    stop("'%s' must be finite value", name);
  }

}

// scalar checks: non-negative finite double
inline void check_nonneg_finite(const double& value, const char* name) {

  if (std::isnan(value) || std::isinf(value) || (value < 0)) {
    stop("'%s' must be non-negative finite value", name);
  }

}

// dimension checks: equality
inline void check_rows_equal(const int& n_rows_a, const int& n_rows_b,
                             const char* name_a, const char* name_b) {

  if (n_rows_a != n_rows_b) {
    stop("number of rows in '%s' must equal the number of rows in '%s'", name_a, name_b);
  }

}

inline void check_cols_equal(const int& n_cols_a, const int& n_cols_b,
                             const char* name_a, const char* name_b) {

  if (n_cols_a != n_cols_b) {
    stop("number of columns in '%s' must equal the number of columns in '%s'", name_a, name_b);
  }

}

inline void check_square(const int& n_rows, const int& n_cols,
                         const char* name) {

  if (n_rows != n_cols) {
    stop("dimensions of '%s' must be square", name);
  }

}

// dimension checks: bounds
inline void check_bounds(const double& lower, const double& upper) {

  check_finite(lower, "lower");
  check_finite(upper, "upper");

  if (lower > upper) {
    stop("value of 'lower' must be less than or equal to value of 'upper'");
  }

}

inline void check_total(const int& n_cols, const double& total,
                        const double& lower, const double& upper) {

  double lower_total = (double)n_cols * lower;
  double upper_total = (double)n_cols * upper;

  check_finite(total, "total");

  if ((total < lower_total) || (total > upper_total)) {
    stop("'total' must be between %f and %f", lower_total, upper_total);
  }

}

// target check for optimization
inline void check_target(const SEXP& x, const SEXP& target,
                         const char* name) {

  if (Rf_isNull(x) || Rf_isNull(target)) {
    stop("requires both '%s' and 'target'", name);
  }
  
}

}

#endif