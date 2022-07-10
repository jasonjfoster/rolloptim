#include "rollport.h"

// [[Rcpp::export(.roll_min_risk)]]
NumericMatrix roll_min_risk(const NumericMatrix& mu, const NumericVector& sigma,
                            const double& total, const double& lower,
                            const double& upper) {
  
  int n_rows_mu = mu.nrow();
  int n_cols_mu = mu.ncol();
  arma::mat arma_mu(mu.begin(), n_rows_mu, n_cols_mu);
  arma::cube arma_sigma(sigma.begin(), n_cols_mu, n_cols_mu, n_rows_mu);
  arma::vec arma_lower(n_cols_mu);
  arma::vec arma_upper(n_cols_mu);
  arma::vec arma_ones(n_cols_mu);
  arma::mat arma_diag(n_cols_mu, n_cols_mu);
  arma::mat arma_weights(n_rows_mu, n_cols_mu);
  
  arma_ones.ones();
  arma_diag.eye();
  arma_lower.fill(lower);
  arma_upper.fill(upper);
  
  // compute rolling portfolio optimization
  RollMinRisk roll_min_risk(arma_mu, arma_sigma, n_rows_mu, n_cols_mu,
                            total, arma_lower, arma_upper,
                            arma_ones, arma_diag, arma_weights);
  parallelFor(0, n_rows_mu, roll_min_risk);
  
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

// [[Rcpp::export(.roll_max_return)]]
NumericMatrix roll_max_return(const NumericMatrix& mu, const NumericVector& sigma,
                              const double& total, const double& lower,
                              const double& upper) {
  
  int n_rows_mu = mu.nrow();
  int n_cols_mu = mu.ncol();
  arma::mat arma_mu(mu.begin(), n_rows_mu, n_cols_mu);
  arma::cube arma_sigma(sigma.begin(), n_cols_mu, n_cols_mu, n_rows_mu);
  arma::vec arma_lower(n_cols_mu);
  arma::vec arma_upper(n_cols_mu);
  arma::vec arma_ones(n_cols_mu);
  arma::mat arma_diag(n_cols_mu, n_cols_mu);
  arma::mat arma_weights(n_rows_mu, n_cols_mu);
  
  // // check 'x' and 'y' arguments for errors
  // check_lm(n_rows_xy, yy.nrow());
  // check both rows and cols (remove "_mu" from names?)
  
  arma_ones.ones();
  arma_diag.eye();
  arma_lower.fill(lower);
  arma_upper.fill(upper);
  
  // compute rolling portfolio optimization
  RollMaxReturn roll_max_return(arma_mu, arma_sigma, n_rows_mu, n_cols_mu,
                                total, arma_lower, arma_upper,
                                arma_ones, arma_diag, arma_weights);
  parallelFor(0, n_rows_mu, roll_max_return);
  
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

// [[Rcpp::export(.roll_max_ratio)]]
NumericMatrix roll_max_ratio(const NumericMatrix& mu, const NumericVector& sigma,
                             const double& gamma, const double& total,
                             const double& lower, const double& upper) {
  
  int n_rows_mu = mu.nrow();
  int n_cols_mu = mu.ncol();
  arma::mat arma_mu(mu.begin(), n_rows_mu, n_cols_mu);
  arma::cube arma_sigma(sigma.begin(), n_cols_mu, n_cols_mu, n_rows_mu);
  arma::vec arma_lower(n_cols_mu);
  arma::vec arma_upper(n_cols_mu);
  arma::vec arma_ones(n_cols_mu);
  arma::mat arma_diag(n_cols_mu, n_cols_mu);
  arma::mat arma_weights(n_rows_mu, n_cols_mu);
  
  arma_ones.ones();
  arma_diag.eye();
  arma_lower.fill(lower);
  arma_upper.fill(upper);
  
  // compute rolling portfolio optimization
  RollMaxRatio roll_max_ratio(arma_mu, arma_sigma, n_rows_mu, n_cols_mu,
                              gamma, total, arma_lower, arma_upper,
                              arma_ones, arma_diag, arma_weights);
  parallelFor(0, n_rows_mu, roll_max_ratio);
  
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