/*** R
library(quantmod)
library(roll)
library(RiskPortfolios)
library(testthat)
library(microbenchmark)

set.seed(5640)
n_vars <- 3
n_obs <- 15
n_size <- n_obs * n_vars
dates <- rev(seq(Sys.Date(), length.out = n_obs, by = "-1 day"))

test_width <- 5
test_total <- 1
test_lower <- 0
test_upper <- 1

test_zoo <- zoo(matrix(rnorm(n_size) / 1000, nrow = n_obs, ncol = n_vars), dates)
colnames(test_zoo) <- paste0("x", 1:n_vars)

mu <- roll_prod(1 + test_zoo, test_width, min_obs = 1) - 1
sigma <- roll_cov(test_zoo, width = test_width, min_obs = 1) * 1000
*/

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;

// 'Worker' function for computing the rolling optimization
struct RollMinRiskEq : public Worker {
  
  const arma::mat arma_mu;      // source
  const arma::cube arma_sigma;
  const int n_rows_mu;
  const int n_cols_mu;
  const double total;
  const arma::vec arma_lower;
  const arma::vec arma_upper;
  const arma::vec arma_ones;
  const arma::mat arma_diag;
  arma::mat& arma_weights;      // destination (pass by reference)
  
  // initialize with source and destination
  RollMinRiskEq(const arma::mat arma_mu, const arma::cube arma_sigma,
                const int n_rows_mu, const int n_cols_mu,
                const double total, const arma::vec arma_lower,
                const arma::vec arma_upper, const arma::vec arma_ones,
                const arma::mat arma_diag, arma::mat& arma_weights)
    : arma_mu(arma_mu), arma_sigma(arma_sigma),
      n_rows_mu(n_rows_mu), n_cols_mu(n_cols_mu),
      total(total), arma_lower(arma_lower),
      arma_upper(arma_upper), arma_ones(arma_ones),
      arma_diag(arma_diag), arma_weights(arma_weights) { }
  
  // function call operator that iterates by slice
  void operator()(std::size_t begin_slice, std::size_t end_slice) {
    for (std::size_t i = begin_slice; i < end_slice; i++) {
      
      int n_solve = 0;
      int n_size = n_cols_mu;
      long double target = max(arma_mu.row(i));
      long double obj = arma::datum::inf;
      long double obj_prev = arma::datum::inf;
      arma::mat mu = arma_mu.row(i);
      arma::mat sigma = arma_sigma.slice(i);
      arma::mat A = sigma;
      arma::vec b(n_size);
      arma::vec weights(n_size);
      
      // total constraint
      n_size += 1;
      
      A = join_rows(A, arma_ones);
      
      b.resize(n_size);
      b(n_size - 1) = total;
      
      // target constraint
      n_size += 1;
      
      A = join_rows(A, trans(mu));
      
      b.resize(n_size);
      b(n_size - 1) = target;
      
      // lower constraints
      n_size += n_cols_mu;
      
      A = join_rows(A, arma_diag);
      
      b.resize(n_size);
      b.subvec(n_size - n_cols_mu, n_size - 1) = arma_lower;
      
      // upper constraints
      n_size += n_cols_mu;
      
      A = join_rows(A, arma_diag);
      
      b.resize(n_size);
      b.subvec(n_size - n_cols_mu, n_size - 1) = arma_upper;
      
      A.resize(n_size, n_size);
      
      // coefficients matrix is symmetric
      arma::mat A_trans = trans(A);
      arma::uvec lower_tri = trimatl_ind(size(A));
      A(lower_tri) = A_trans(lower_tri);
      
      // number of index combinations
      for (int j = 0; j < pow((long double)2.0, (long double)2.0 * n_cols_mu); j++) {
        
        n_size = n_cols_mu + 2;
        arma::uvec arma_ix = arma::linspace<arma::uvec>(0, n_size - 1, n_size);
        
        // find the index combination
        for (int k = 0; k < 2 * n_cols_mu; k++) {
          
          if (!(j & (int)pow((long double)2.0, k))) {
            
            n_size += 1;
            
            arma_ix.resize(n_size);
            arma_ix(n_size - 1) = n_cols_mu + 2 + k;
            
          }
          
        }
        
        arma::mat A_subset = A(arma_ix, arma_ix);
        arma::vec b_subset = b(arma_ix);
        
        // check if solution is found 
        bool status_solve = arma::solve(weights, A_subset, b_subset,
                                        arma::solve_opts::no_approx);
        
        // don't find approximate solution for rank deficient system
        if (status_solve) {
          
          arma_ix = arma::linspace<arma::uvec>(0, n_cols_mu - 1, n_cols_mu);
          
          // weights
          arma::vec weights_subset = weights(arma_ix);
          arma::mat trans_weights = trans(weights_subset);
          
          // check if constraints are satisfied
          if (all(weights_subset >= arma_lower[0]) && all(weights_subset <= arma_upper[0])) {
            
            n_solve += 1;
            
            // objective value
            obj = as_scalar(trans_weights * sigma * weights_subset);
            
            if (obj <= obj_prev) {
              
              obj_prev = obj;
              
              arma_weights.row(i) = trans_weights;
              
            }
            
          }
          
        }
        
      }
      
      if (n_solve == 0) {
        
        arma::rowvec no_solution(n_cols_mu);
        no_solution.fill(NA_REAL);
        
        arma_weights.row(i) = no_solution;
        
      }
      
    }
  }
  
};

// 'Worker' function for computing the rolling optimization
struct RollMinRiskGeq : public Worker {

  const arma::mat arma_mu;      // source
  const arma::cube arma_sigma;
  const int n_rows_mu;
  const int n_cols_mu;
  const double total;
  const arma::vec arma_lower;
  const arma::vec arma_upper;
  const arma::vec arma_ones;
  const arma::mat arma_diag;
  arma::mat& arma_weights;      // destination (pass by reference)

  // initialize with source and destination
  RollMinRiskGeq(const arma::mat arma_mu, const arma::cube arma_sigma,
                 const int n_rows_mu, const int n_cols_mu,
                 const double total, const arma::vec arma_lower,
                 const arma::vec arma_upper, const arma::vec arma_ones,
                 const arma::mat arma_diag, arma::mat& arma_weights)
    : arma_mu(arma_mu), arma_sigma(arma_sigma),
      n_rows_mu(n_rows_mu), n_cols_mu(n_cols_mu),
      total(total), arma_lower(arma_lower),
      arma_upper(arma_upper), arma_ones(arma_ones),
      arma_diag(arma_diag), arma_weights(arma_weights) { }

  // function call operator that iterates by slice
  void operator()(std::size_t begin_slice, std::size_t end_slice) {
    for (std::size_t i = begin_slice; i < end_slice; i++) {
      
      int n_solve = 0;
      int n_size = n_cols_mu;
      long double target = -arma::datum::inf;
      long double obj = arma::datum::inf;
      long double obj_prev = arma::datum::inf;
      arma::mat mu = arma_mu.row(i);
      arma::mat sigma = arma_sigma.slice(i);
      arma::mat A = sigma;
      arma::vec b(n_size);
      arma::vec weights(n_size);
      
      // total constraint
      n_size += 1;
      
      A = join_rows(A, arma_ones);
      
      b.resize(n_size);
      b(n_size - 1) = total;
      
      // target constraint
      n_size += 1;
      
      A = join_rows(A, trans(mu));
      
      b.resize(n_size);
      b(n_size - 1) = target;

      // lower constraints
      n_size += n_cols_mu;
      
      A = join_rows(A, arma_diag);
      
      b.resize(n_size);
      b.subvec(n_size - n_cols_mu, n_size - 1) = arma_lower;
      
      // upper constraints
      n_size += n_cols_mu;
      
      A = join_rows(A, arma_diag);
      
      b.resize(n_size);
      b.subvec(n_size - n_cols_mu, n_size - 1) = arma_upper;
      
      A.resize(n_size, n_size);
      
      // coefficients matrix is symmetric
      arma::mat A_trans = trans(A);
      arma::uvec lower_tri = trimatl_ind(size(A));
      A(lower_tri) = A_trans(lower_tri);
      
      // number of index combinations
      for (int j = 0; j < pow((long double)2.0, (long double)2.0 * n_cols_mu + 1); j++) {
        
        n_size = n_cols_mu + 1;
        arma::uvec arma_ix = arma::linspace<arma::uvec>(0, n_size - 1, n_size);
        
        // find the index combination
        for (int k = 0; k < 2 * n_cols_mu + 1; k++) {
          
          if (!(j & (int)pow((long double)2.0, k))) {
            
            n_size += 1;
            
            arma_ix.resize(n_size);
            arma_ix(n_size - 1) = n_cols_mu + 1 + k;
            
          }
          
        }
        
        arma::mat A_subset = A(arma_ix, arma_ix);
        arma::vec b_subset = b(arma_ix);
        
        // check if solution is found 
        bool status_solve = arma::solve(weights, A_subset, b_subset,
                                        arma::solve_opts::no_approx);
        
        // don't find approximate solution for rank deficient system
        if (status_solve) {
          
          arma_ix = arma::linspace<arma::uvec>(0, n_cols_mu - 1, n_cols_mu);
          
          // weights
          arma::vec weights_subset = weights(arma_ix);
          arma::mat trans_weights = trans(weights_subset);
          
          // check if constraints are satisfied
          if (all(weights_subset >= arma_lower[0]) && all(weights_subset <= arma_upper[0]) &&
              (sum(mu * weights_subset) >= target)) {
            
            n_solve += 1;
            
            // objective value
            obj = as_scalar(trans_weights * sigma * weights_subset);
            
            if (obj <= obj_prev) {
              
              obj_prev = obj;
              
              arma_weights.row(i) = trans_weights;
              
            }
            
          }

        }
        
      }
      
      if (n_solve == 0) {
        
        arma::rowvec no_solution(n_cols_mu);
        no_solution.fill(NA_REAL);
        
        arma_weights.row(i) = no_solution;
        
      }

    }
  }

};

// https://stackoverflow.com/a/31725653
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
  
  // compute rolling optimizations
  RollMinRiskEq roll_max_return(arma_mu, arma_sigma, n_rows_mu, n_cols_mu,
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
  
  // compute rolling optimizations
  RollMinRiskGeq roll_min_risk(arma_mu, arma_sigma, n_rows_mu, n_cols_mu,
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

/*** R
##' Rolling Portfolio Optimizations to Minimize Risk
##'
##' A function for computing the rolling portfolio optimizations to minimize risk.
##' 
##' @param mu matrix. Rows are returns and columns are variables.
##' @param sigma cube. Slices are covariance matrices.
##' @examples

roll_result <- .roll_max_return(mu, sigma, test_total, test_lower, test_upper)
roll_result <- .roll_min_risk(mu, sigma, test_total, test_lower, test_upper)
*/

/*** R
rollapplyr_cube <- function(f, type, mu, sigma, total = 1, lower = 0, upper = 1) {
  
  mu_attr <- attributes(mu)
  
  n_rows_mu <- nrow(mu)
  n_cols_mu <- ncol(mu)
  result <- matrix(as.numeric(NA), n_rows_mu, n_cols_mu)
  
  for (i in 2:n_rows_mu) {
    
    status_solve <- tryCatch(f(mu[i, ], Sigma = sigma[ , , i],
                               control = list(type = type, constraint = "user",
                                              gamma = sqrt(.Machine$double.eps),
                                              LB = rep(lower, n_cols_mu),
                                              UB = rep(upper, n_cols_mu))),
                             
                             error = function(x) {
                               rep(NA, n_cols_mu)
                             })
    
    result[i, ] <- status_solve
    
  }
  
  attributes(result) <- mu_attr
  
  return(result)
  
}
*/

/*** R
rp_result <- rollapplyr_cube(optimalPortfolio, "mv", mu, sigma,
                             test_total, test_lower, test_upper)

rp_result <- rollapplyr_cube(optimalPortfolio, "minvol", mu, sigma,
                             test_total, test_lower, test_upper)
*/

/*** R
expect_equal(roll_result[4:15], rp_result[4:15])
*/

/*** R
microbenchmark("roll" = .roll_min_risk(mu, sigma,  test_total, test_lower, test_upper),
               "rp" = rollapplyr_cube(optimalPortfolio, mu, sigma,
                                      test_total, test_lower, test_upper))
*/