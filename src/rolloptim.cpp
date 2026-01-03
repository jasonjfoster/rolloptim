#include "rolloptim.h"

void check_rows(const int& n_rows_mu, const int& n_slices_sigma) {
  
  if (n_rows_mu != n_slices_sigma) {
    stop("number of rows in 'mu' must equal the number of slices in 'sigma'");
  }
  
}

void check_cols(const int& n_cols_mu, const int& n_cols_sigma) {
  
  if (n_cols_mu != n_cols_sigma) {
    stop("number of columns in 'mu' must equal the number of columns in 'sigma'");
  }
  
}

void check_sigma(const int& n_rows_sigma, const int& n_cols_sigma) {
  
  if (n_rows_sigma != n_cols_sigma) {
    stop("dimensions of 'sigma' must be square");
  }
  
}

void check_target(const SEXP& mu, const SEXP& target,
                  const char* name) {
  
  if (Rf_isNull(mu) || Rf_isNull(target)) {
    stop("requires both '%s' and 'target'", name);
  }
  
}

// [[Rcpp::export(.roll_min_var)]]
NumericMatrix roll_min_var(const NumericVector& sigma, const SEXP& mu,
                           const SEXP& target, const double& total,
                           const double& lower, const double& upper) {

  if (Rf_isNull(mu) && Rf_isNull(target)) {

    SEXP sexp_dim_sigma = sigma.attr("dim");

    // if (Rf_isNull(sexp_dim_sigma)) {
    //   stop("'sigma' must be cube with slices of covariance matrices");
    // }

    IntegerVector dim_sigma = sigma.attr("dim");
    int n_rows = dim_sigma[2];
    int n_cols = dim_sigma[1];

    // check 'sigma' argument for errors
    check_sigma(dim_sigma[0], dim_sigma[1]);

    int n_size = n_cols + 1 + n_cols + n_cols;
    arma::cube arma_sigma(sigma.begin(), n_cols, n_cols, n_rows);
    arma::vec arma_lower(n_cols);
    arma::vec arma_upper(n_cols);
    arma::vec arma_ones(n_cols);
    arma::mat arma_diag(n_cols, n_cols);
    arma::mat arma_A(n_size, n_size);
    arma::vec arma_b(n_size);
    arma::mat arma_weights(n_rows, n_cols);
    
    arma_ones.ones();
    arma_diag.eye();
    arma_lower.fill(lower);
    arma_upper.fill(upper);
    
    // total constraint
    arma_A.submat(n_cols, 0, n_cols, n_cols - 1) = trans(arma_ones);
    arma_A.submat(0, n_cols, n_cols - 1, n_cols) = arma_ones;
    arma_b(n_cols) = total;
    
    // lower constraints
    arma_A.submat(n_cols + 1, 0,
                  n_cols + 1 + n_cols - 1, n_cols - 1) = arma_diag;
    arma_A.submat(0, n_cols + 1,
                  n_cols - 1, n_cols + 1 + n_cols - 1) = arma_diag;
    arma_b.subvec(n_cols + 1, n_cols + 1 + n_cols - 1) = arma_lower;
    
    // upper constraints
    arma_A.submat(n_cols + 1 + n_cols, 0,
                  n_cols + 1 + n_cols + n_cols - 1, n_cols - 1) = arma_diag;
    arma_A.submat(0, n_cols + 1 + n_cols,
                  n_cols - 1, n_cols + 1 + n_cols + n_cols - 1) = arma_diag;
    arma_b.subvec(n_cols + 1 + n_cols,
                  n_cols + 1 + n_cols + n_cols - 1) = arma_upper;
    
    // compute rolling optimizations
    rolloptim::RollMinVarMuFALSE roll_min_var(arma_sigma, n_rows, n_cols,
                                              total, arma_lower, arma_upper, arma_ones, arma_diag,
                                              arma_A, arma_b,
                                              arma_weights);
    parallelFor(0, n_rows, roll_min_var);
    
    NumericMatrix result(wrap(arma_weights));
    List dimnames_sigma = sigma.attr("dimnames");
    if (dimnames_sigma.size() > 1) {
      result.attr("dimnames") = List::create(R_NilValue, dimnames_sigma[1]);
    }
    
    return result;

  } else {

    check_target(mu, target, "mu");

    NumericMatrix rcpp_mu(mu);
    NumericVector rcpp_target(target);

    int n_rows = rcpp_mu.nrow();
    int n_cols = rcpp_mu.ncol();
    int n_size = n_cols + 2 + n_cols + n_cols;
    arma::cube arma_sigma;
    SEXP sexp_dim_sigma = sigma.attr("dim");

    if (!Rf_isNull(sexp_dim_sigma)) {

      IntegerVector dim_sigma(sexp_dim_sigma);

      check_rows(n_rows, dim_sigma[2]);
      check_cols(n_cols, dim_sigma[1]);
      check_sigma(dim_sigma[0], dim_sigma[1]);

      arma_sigma = arma::cube(sigma.begin(), n_cols, n_cols, n_rows);

    } else {
      
      // if (n_cols != 1) {
      //   stop("number of columns in 'mu' must be one for univariate sigma");
      // }

      check_rows(n_rows, sigma.size());

      arma_sigma = arma::cube(1, 1, n_rows);
      for (int i = 0; i < n_rows; i++) {
        arma_sigma(0, 0, i) = sigma[i];
      }

    }

    arma::mat arma_mu(rcpp_mu.begin(), n_rows, n_cols);
    arma::vec arma_target(rcpp_target.begin(), rcpp_target.size());
    arma::vec arma_lower(n_cols);
    arma::vec arma_upper(n_cols);
    arma::vec arma_ones(n_cols);
    arma::mat arma_diag(n_cols, n_cols);
    arma::mat arma_A(n_size, n_size);
    arma::vec arma_b(n_size);
    arma::mat arma_weights(n_rows, n_cols);

    arma_ones.ones();
    arma_diag.eye();
    arma_lower.fill(lower);
    arma_upper.fill(upper);

    // total constraint
    arma_A.submat(n_cols, 0, n_cols, n_cols - 1) = trans(arma_ones);
    arma_A.submat(0, n_cols, n_cols - 1, n_cols) = arma_ones;
    arma_b(n_cols) = total;

    // mu constraint
    arma_A.submat(n_cols + 1, 0, n_cols + 1, n_cols - 1).zeros();
    arma_A.submat(0, n_cols + 1, n_cols - 1, n_cols + 1).zeros();
    arma_b(n_cols + 1) = 0;
    
    // lower constraints
    arma_A.submat(n_cols + 2, 0,
                  n_cols + 2 + n_cols - 1, n_cols - 1) = arma_diag;
    arma_A.submat(0, n_cols + 2,
                  n_cols - 1, n_cols + 2 + n_cols - 1) = arma_diag;
    arma_b.subvec(n_cols + 2, n_cols + 2 + n_cols - 1) = arma_lower;
    
    // upper constraints
    arma_A.submat(n_cols + 2 + n_cols, 0,
                  n_cols + 2 + n_cols + n_cols - 1, n_cols - 1) = arma_diag;
    arma_A.submat(0, n_cols + 2 + n_cols,
                  n_cols - 1, n_cols + 2 + n_cols + n_cols - 1) = arma_diag;
    arma_b.subvec(n_cols + 2 + n_cols,
                  n_cols + 2 + n_cols + n_cols - 1) = arma_upper;
    
    // compute rolling optimizations
    rolloptim::RollMinVarMuTRUE roll_min_var(arma_sigma, arma_mu, n_rows, n_cols,
                                             total, arma_target,
                                             arma_lower, arma_upper, arma_ones, arma_diag,
                                             arma_A, arma_b,
                                             arma_weights);
    parallelFor(0, n_rows, roll_min_var);

    NumericMatrix result(wrap(arma_weights));
    List dimnames = rcpp_mu.attr("dimnames");
    result.attr("dimnames") = dimnames;
    result.attr("index") = rcpp_mu.attr("index");
    result.attr(".indexCLASS") = rcpp_mu.attr(".indexCLASS");
    result.attr(".indexTZ") = rcpp_mu.attr(".indexTZ");
    result.attr("tclass") = rcpp_mu.attr("tclass");
    result.attr("tzone") = rcpp_mu.attr("tzone");
    result.attr("class") = rcpp_mu.attr("class");
    
    return result;

  }
  
}

// [[Rcpp::export(.roll_max_mean)]]
NumericMatrix roll_max_mean(const NumericMatrix& mu, const SEXP& sigma,
                            const SEXP& target, const double& total,
                            const double& lower, const double& upper) {

  if (!Rf_isNull(sigma) || !Rf_isNull(target)) {
    warning("'target' is only supported for linear constraints");
  }
  
  int n_rows = mu.nrow();
  int n_cols = mu.ncol();
  int n_size = n_cols + 1 + n_cols + n_cols;
  arma::mat arma_mu(mu.begin(), n_rows, n_cols);
  arma::vec arma_lower(n_cols);
  arma::vec arma_upper(n_cols);
  arma::vec arma_ones(n_cols);
  arma::mat arma_diag(n_cols, n_cols);
  arma::mat arma_A(n_size, n_size);
  arma::vec arma_b(n_size);
  arma::mat arma_weights(n_rows, n_cols);
  
  arma_ones.ones();
  arma_diag.eye();
  arma_lower.fill(lower);
  arma_upper.fill(upper);
  
  // total constraint
  arma_A.submat(n_cols, 0, n_cols, n_cols - 1) = trans(arma_ones);
  arma_A.submat(0, n_cols, n_cols - 1, n_cols) = arma_ones;
  arma_b(n_cols) = total;
  
  // lower constraints
  arma_A.submat(n_cols + 1, 0,
                n_cols + 1 + n_cols - 1, n_cols - 1) = arma_diag;
  arma_A.submat(0, n_cols + 1,
                n_cols - 1, n_cols + 1 + n_cols - 1) = arma_diag;
  arma_b.subvec(n_cols + 1, n_cols + 1 + n_cols - 1) = arma_lower;
  
  // upper constraints
  arma_A.submat(n_cols + 1 + n_cols, 0,
                n_cols + 1 + n_cols + n_cols - 1, n_cols - 1) = arma_diag;
  arma_A.submat(0, n_cols + 1 + n_cols,
                n_cols - 1, n_cols + 1 + n_cols + n_cols - 1) = arma_diag;
  arma_b.subvec(n_cols + 1 + n_cols,
                n_cols + 1 + n_cols + n_cols - 1) = arma_upper;
  
  // compute rolling optimizations
  rolloptim::RollMaxMeanSigmaFALSE roll_max_mean(arma_mu, n_rows, n_cols,
                                                 total, arma_lower, arma_upper, arma_ones, arma_diag,
                                                 arma_A, arma_b, arma_weights);
  parallelFor(0, n_rows, roll_max_mean);
  
  NumericMatrix result(wrap(arma_weights));
  List dimnames = mu.attr("dimnames");
  result.attr("dimnames") = dimnames;
  result.attr("index") = mu.attr("index");
  result.attr(".indexCLASS") = mu.attr(".indexCLASS");
  result.attr(".indexTZ") = mu.attr(".indexTZ");
  result.attr("tclass") = mu.attr("tclass");
  result.attr("tzone") = mu.attr("tzone");
  result.attr("class") = mu.attr("class");
  
  return result;
  
}

// [[Rcpp::export(.roll_max_utility)]]
NumericMatrix roll_max_utility(const NumericMatrix& mu, const NumericVector& sigma,
                               const double& lambda, const double& total,
                               const double& lower, const double& upper) {
  
  int n_rows = mu.nrow();
  int n_cols = mu.ncol();
  int n_size = n_cols + 1 + n_cols + n_cols;
  arma::cube arma_sigma;
  SEXP sexp_dim_sigma = sigma.attr("dim");

  if (!Rf_isNull(sexp_dim_sigma)) {

    IntegerVector dim_sigma(sexp_dim_sigma);

    check_rows(n_rows, dim_sigma[2]);
    check_cols(n_cols, dim_sigma[1]);
    check_sigma(dim_sigma[0], dim_sigma[1]);

    arma_sigma = arma::cube(sigma.begin(), n_cols, n_cols, n_rows);

  } else {

      // if (n_cols != 1) {
      //   stop("number of columns in 'mu' must be one for univariate sigma");
      // }

    check_rows(n_rows, sigma.size());

    arma_sigma = arma::cube(1, 1, n_rows);
    for (int i = 0; i < n_rows; i++) {
      arma_sigma(0, 0, i) = sigma[i];
    }

  }

  arma::mat arma_mu(mu.begin(), n_rows, n_cols);
  arma::vec arma_lower(n_cols);
  arma::vec arma_upper(n_cols);
  arma::vec arma_ones(n_cols);
  arma::mat arma_diag(n_cols, n_cols);
  arma::mat arma_A(n_size, n_size);
  arma::vec arma_b(n_size);
  arma::mat arma_weights(n_rows, n_cols);
  
  arma_ones.ones();
  arma_diag.eye();
  arma_lower.fill(lower);
  arma_upper.fill(upper);
  
  // total constraint
  arma_A.submat(n_cols, 0, n_cols, n_cols - 1) = trans(arma_ones);
  arma_A.submat(0, n_cols, n_cols - 1, n_cols) = arma_ones;
  arma_b(n_cols) = total;
  
  // lower constraints
  arma_A.submat(n_cols + 1, 0,
                n_cols + 1 + n_cols - 1, n_cols - 1) = arma_diag;
  arma_A.submat(0, n_cols + 1,
                n_cols - 1, n_cols + 1 + n_cols - 1) = arma_diag;
  arma_b.subvec(n_cols + 1, n_cols + 1 + n_cols - 1) = arma_lower;
  
  // upper constraints
  arma_A.submat(n_cols + 1 + n_cols, 0,
                n_cols + 1 + n_cols + n_cols - 1, n_cols - 1) = arma_diag;
  arma_A.submat(0, n_cols + 1 + n_cols,
                n_cols - 1, n_cols + 1 + n_cols + n_cols - 1) = arma_diag;
  arma_b.subvec(n_cols + 1 + n_cols,
                n_cols + 1 + n_cols + n_cols - 1) = arma_upper;
  
  // compute rolling optimizations
  rolloptim::RollMaxUtility roll_max_utility(arma_mu, arma_sigma, n_rows, n_cols,
                                             lambda, total, arma_lower, arma_upper, arma_ones, arma_diag,
                                             arma_A, arma_b, arma_weights);
  parallelFor(0, n_rows, roll_max_utility);
  
  NumericMatrix result(wrap(arma_weights));
  List dimnames = mu.attr("dimnames");
  result.attr("dimnames") = dimnames;
  result.attr("index") = mu.attr("index");
  result.attr(".indexCLASS") = mu.attr(".indexCLASS");
  result.attr(".indexTZ") = mu.attr(".indexTZ");
  result.attr("tclass") = mu.attr("tclass");
  result.attr("tzone") = mu.attr("tzone");
  result.attr("class") = mu.attr("class");
  
  return result;
  
}

// [[Rcpp::export(.roll_min_rss)]]
NumericMatrix roll_min_rss(const NumericVector& xx, const NumericVector& xy,
                           const double& total, const double& lower,
                           const double& upper) {
  
  int n_rows;
  int n_cols;
  SEXP sexp_dim_xx = xx.attr("dim");
  SEXP sexp_dim_xy = xy.attr("dim");

  if (!Rf_isNull(sexp_dim_xx) && !Rf_isNull(sexp_dim_xy)) {

    IntegerVector dim_xx(sexp_dim_xx);
    IntegerVector dim_xy(sexp_dim_xy);

    // if (dim_xx.size() == 3) {
    
      n_rows = dim_xx[2];
      n_cols = dim_xx[1];

    // } else if (dim_xx.size() == 2) {

    //     n_rows = dim_xx[0];
    //     n_cols = dim_xx[1];

    // }

    check_sigma(dim_xx[0], dim_xx[1]);

    // if (dim_xy.size() == 3) {

      check_rows(n_rows, dim_xy[2]);
      check_cols(n_cols, dim_xy[0]);

    // } else if (dim_xy.size() == 2) {

    //   check_rows(n_rows, dim_xy[1]);
    //   check_cols(n_cols, dim_xy[0]);

    // }

  } else {

    n_rows = xx.size();
    n_cols = 1;   

  }
  
  int n_size = n_cols + 1 + n_cols + n_cols;
  arma::cube arma_xx(xx.begin(), n_cols, n_cols, n_rows);
  arma::mat arma_xy(xy.begin(), n_cols, n_rows);
  arma::vec arma_lower(n_cols);
  arma::vec arma_upper(n_cols);
  arma::vec arma_ones(n_cols);
  arma::mat arma_diag(n_cols, n_cols);
  arma::mat arma_A(n_size, n_size);
  arma::vec arma_b(n_size);
  arma::mat arma_weights(n_rows, n_cols);

  arma_ones.ones();
  arma_diag.eye();
  arma_lower.fill(lower);
  arma_upper.fill(upper);
  
  // total constraint
  arma_A.submat(n_cols, 0, n_cols, n_cols - 1) = trans(arma_ones);
  arma_A.submat(0, n_cols, n_cols - 1, n_cols) = arma_ones;
  arma_b(n_cols) = total;
  
  // lower constraints
  arma_A.submat(n_cols + 1, 0,
                n_cols + 1 + n_cols - 1, n_cols - 1) = arma_diag;
  arma_A.submat(0, n_cols + 1,
                n_cols - 1, n_cols + 1 + n_cols - 1) = arma_diag;
  arma_b.subvec(n_cols + 1, n_cols + 1 + n_cols - 1) = arma_lower;
  
  // upper constraints
  arma_A.submat(n_cols + 1 + n_cols, 0,
                n_cols + 1 + n_cols + n_cols - 1, n_cols - 1) = arma_diag;
  arma_A.submat(0, n_cols + 1 + n_cols,
                n_cols - 1, n_cols + 1 + n_cols + n_cols - 1) = arma_diag;
  arma_b.subvec(n_cols + 1 + n_cols,
                n_cols + 1 + n_cols + n_cols - 1) = arma_upper;
  
  // compute rolling optimizations
  rolloptim::RollMaxUtility roll_max_utility(trans(arma_xy), arma_xx, n_rows, n_cols,
                                             1, total, arma_lower, arma_upper, arma_ones, arma_diag,
                                             arma_A, arma_b, arma_weights);
  parallelFor(0, n_rows, roll_max_utility);
  
  NumericMatrix result(wrap(arma_weights));
  List dimnames_xx = xx.attr("dimnames");
  if (dimnames_xx.size() > 1) {
    result.attr("dimnames") = List::create(R_NilValue, dimnames_xx[1]);
  }
  
  return result;
  
}