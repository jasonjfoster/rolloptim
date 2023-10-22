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

// [[Rcpp::export(.roll_min_var)]]
NumericMatrix roll_min_var(const NumericVector& sigma, const double& total,
                           const double& lower, const double& upper) {
  
  IntegerVector dim_sigma = sigma.attr("dim");
  int n_rows_sigma = dim_sigma[2];
  int n_cols_sigma = dim_sigma[1];
  int n_size = n_cols_sigma + 1 + n_cols_sigma + n_cols_sigma;
  arma::cube arma_sigma(sigma.begin(), n_cols_sigma, n_cols_sigma, n_rows_sigma);
  arma::vec arma_lower(n_cols_sigma);
  arma::vec arma_upper(n_cols_sigma);
  arma::vec arma_ones(n_cols_sigma);
  arma::mat arma_diag(n_cols_sigma, n_cols_sigma);
  arma::mat arma_A(n_size, n_size);
  arma::vec arma_b(n_size);
  arma::mat arma_weights(n_rows_sigma, n_cols_sigma);
  
  // check 'sigma' argument for errors
  check_sigma(dim_sigma[0], dim_sigma[1]);
  
  arma_ones.ones();
  arma_diag.eye();
  arma_lower.fill(lower);
  arma_upper.fill(upper);
  
  // total constraint
  arma_A.submat(n_cols_sigma, 0, n_cols_sigma, n_cols_sigma - 1) = trans(arma_ones);
  arma_A.submat(0, n_cols_sigma, n_cols_sigma - 1, n_cols_sigma) = arma_ones;
  arma_b(n_cols_sigma) = total;
  
  // lower constraints
  arma_A.submat(n_cols_sigma + 1, 0,
                n_cols_sigma + 1 + n_cols_sigma - 1, n_cols_sigma - 1) = arma_diag;
  arma_A.submat(0, n_cols_sigma + 1,
                n_cols_sigma - 1, n_cols_sigma + 1 + n_cols_sigma - 1) = arma_diag;
  arma_b.subvec(n_cols_sigma + 1, n_cols_sigma + 1 + n_cols_sigma - 1) = arma_lower;
  
  // upper constraints
  arma_A.submat(n_cols_sigma + 1 + n_cols_sigma, 0,
                n_cols_sigma + 1 + n_cols_sigma + n_cols_sigma - 1, n_cols_sigma - 1) = arma_diag;
  arma_A.submat(0, n_cols_sigma + 1 + n_cols_sigma,
                n_cols_sigma - 1, n_cols_sigma + 1 + n_cols_sigma + n_cols_sigma - 1) = arma_diag;
  arma_b.subvec(n_cols_sigma + 1 + n_cols_sigma,
                n_cols_sigma + 1 + n_cols_sigma + n_cols_sigma - 1) = arma_upper;
  
  // compute rolling optimizations
  rolloptim::RollMinVar roll_min_var(arma_sigma, n_rows_sigma, n_cols_sigma,
                                     total, arma_lower, arma_upper, arma_ones, arma_diag,
                                     arma_A, arma_b, arma_weights);
  parallelFor(0, n_rows_sigma, roll_min_var);
  
  NumericMatrix result(wrap(arma_weights));
  List dimnames_sigma = sigma.attr("dimnames");
  if (dimnames_sigma.size() > 1) {
    result.attr("dimnames") = List::create(R_NilValue, dimnames_sigma[1]);
  }
  
  return result;
  
}

// [[Rcpp::export(.roll_max_mean)]]
NumericMatrix roll_max_mean(const NumericMatrix& mu, const double& total,
                            const double& lower, const double& upper) {
  
  int n_rows_mu = mu.nrow();
  int n_cols_mu = mu.ncol();
  int n_size = n_cols_mu + 1 + n_cols_mu + n_cols_mu;
  arma::mat arma_mu(mu.begin(), n_rows_mu, n_cols_mu);
  arma::vec arma_lower(n_cols_mu);
  arma::vec arma_upper(n_cols_mu);
  arma::vec arma_ones(n_cols_mu);
  arma::mat arma_diag(n_cols_mu, n_cols_mu);
  arma::mat arma_A(n_size, n_size);
  arma::vec arma_b(n_size);
  arma::mat arma_weights(n_rows_mu, n_cols_mu);
  
  arma_ones.ones();
  arma_diag.eye();
  arma_lower.fill(lower);
  arma_upper.fill(upper);
  
  // total constraint
  arma_A.submat(n_cols_mu, 0, n_cols_mu, n_cols_mu - 1) = trans(arma_ones);
  arma_A.submat(0, n_cols_mu, n_cols_mu - 1, n_cols_mu) = arma_ones;
  arma_b(n_cols_mu) = total;
  
  // lower constraints
  arma_A.submat(n_cols_mu + 1, 0,
                n_cols_mu + 1 + n_cols_mu - 1, n_cols_mu - 1) = arma_diag;
  arma_A.submat(0, n_cols_mu + 1,
                n_cols_mu - 1, n_cols_mu + 1 + n_cols_mu - 1) = arma_diag;
  arma_b.subvec(n_cols_mu + 1, n_cols_mu + 1 + n_cols_mu - 1) = arma_lower;
  
  // upper constraints
  arma_A.submat(n_cols_mu + 1 + n_cols_mu, 0,
                n_cols_mu + 1 + n_cols_mu + n_cols_mu - 1, n_cols_mu - 1) = arma_diag;
  arma_A.submat(0, n_cols_mu + 1 + n_cols_mu,
                n_cols_mu - 1, n_cols_mu + 1 + n_cols_mu + n_cols_mu - 1) = arma_diag;
  arma_b.subvec(n_cols_mu + 1 + n_cols_mu,
                n_cols_mu + 1 + n_cols_mu + n_cols_mu - 1) = arma_upper;
  
  // compute rolling optimizations
  rolloptim::RollMaxMean roll_max_mean(arma_mu, n_rows_mu, n_cols_mu,
                                       total, arma_lower, arma_upper, arma_ones, arma_diag,
                                       arma_A, arma_b, arma_weights);
  parallelFor(0, n_rows_mu, roll_max_mean);
  
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
  
  int n_rows_mu = mu.nrow();
  int n_cols_mu = mu.ncol();
  int n_size = n_cols_mu + 1 + n_cols_mu + n_cols_mu;
  IntegerVector dim_sigma = sigma.attr("dim");
  arma::mat arma_mu(mu.begin(), n_rows_mu, n_cols_mu);
  arma::cube arma_sigma(sigma.begin(), n_cols_mu, n_cols_mu, n_rows_mu);
  arma::vec arma_lower(n_cols_mu);
  arma::vec arma_upper(n_cols_mu);
  arma::vec arma_ones(n_cols_mu);
  arma::mat arma_diag(n_cols_mu, n_cols_mu);
  arma::mat arma_A(n_size, n_size);
  arma::vec arma_b(n_size);
  arma::mat arma_weights(n_rows_mu, n_cols_mu);
  
  // check 'mu' and 'sigma' arguments for errors
  if (dim_sigma.size() == 3) {
    check_rows(n_rows_mu, dim_sigma[2]);
  } else {
    check_rows(n_rows_mu, 1);
  }
  
  check_cols(n_cols_mu, dim_sigma[1]);
  check_sigma(dim_sigma[0], dim_sigma[1]);
  
  arma_ones.ones();
  arma_diag.eye();
  arma_lower.fill(lower);
  arma_upper.fill(upper);
  
  // total constraint
  arma_A.submat(n_cols_mu, 0, n_cols_mu, n_cols_mu - 1) = trans(arma_ones);
  arma_A.submat(0, n_cols_mu, n_cols_mu - 1, n_cols_mu) = arma_ones;
  arma_b(n_cols_mu) = total;
  
  // lower constraints
  arma_A.submat(n_cols_mu + 1, 0,
                n_cols_mu + 1 + n_cols_mu - 1, n_cols_mu - 1) = arma_diag;
  arma_A.submat(0, n_cols_mu + 1,
                n_cols_mu - 1, n_cols_mu + 1 + n_cols_mu - 1) = arma_diag;
  arma_b.subvec(n_cols_mu + 1, n_cols_mu + 1 + n_cols_mu - 1) = arma_lower;
  
  // upper constraints
  arma_A.submat(n_cols_mu + 1 + n_cols_mu, 0,
                n_cols_mu + 1 + n_cols_mu + n_cols_mu - 1, n_cols_mu - 1) = arma_diag;
  arma_A.submat(0, n_cols_mu + 1 + n_cols_mu,
                n_cols_mu - 1, n_cols_mu + 1 + n_cols_mu + n_cols_mu - 1) = arma_diag;
  arma_b.subvec(n_cols_mu + 1 + n_cols_mu,
                n_cols_mu + 1 + n_cols_mu + n_cols_mu - 1) = arma_upper;
  
  // compute rolling optimizations
  rolloptim::RollMaxUtility roll_max_utility(arma_mu, arma_sigma, n_rows_mu, n_cols_mu,
                                             lambda, total, arma_lower, arma_upper, arma_ones, arma_diag,
                                             arma_A, arma_b, arma_weights);
  parallelFor(0, n_rows_mu, roll_max_utility);
  
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