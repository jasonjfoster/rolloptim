/*** R
library(quantmod)
library(roll)
library(CVXR)
library(testthat)
library(microbenchmark)

set.seed(5640)
n_vars <- 3
n_obs <- 15
width <- 5

target <- -1
total <- 1
lower <- 0
upper <- 1

n_size <- n_obs * n_vars
dates <- rev(seq(Sys.Date(), length.out = n_obs, by = "-1 day"))
test_zoo <- zoo(matrix(rnorm(n_size), nrow = n_obs, ncol = n_vars), dates)
colnames(test_zoo) <- paste0("x", 1:n_vars)

mu <- roll_prod(1 + test_zoo, width, min_obs = 1) - 1
sigma <- roll_cov(test_zoo, width = width, min_obs = 1)

# tickers <- c("ITOT", "TLT", "IEF", "IYR")
# cols <- c("Weight", "Return", "Risk", "Ratio", "Stress")
# width_st <- 252
# width_lt <- 252 * 5
# width_z <- 252 * 30
# scale <- list("periods" = 252, "overlap" = 5)
# start <- rep(1, length(tickers))
# 
# getSymbols(tickers, src = "tiingo", from = "1900-01-01", adjust = TRUE, api.key = Sys.getenv("TIINGO_API_KEY"))
# prices_xts <- do.call(merge, c(lapply(tickers, function(i) Cl(get(i))), all = TRUE))
# colnames(prices_xts) <- tickers
# 
# returns_xts <- na.omit(diff(log(prices_xts)))
# overlap_xts <- roll_mean(returns_xts, scale[["overlap"]], min_obs = 1)
# 
# mu_xts <- roll_prod(1 + returns_xts, width_st, min_obs = 1) - 1
# sigma_xts <- roll_cov(overlap_xts, width = width_lt, min_obs = 1) * scale[["periods"]] * scale[["overlap"]]
*/

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;

// 'Worker' function for computing the rolling statistic using an online algorithm
struct RollMinRisk : public Worker {

  const arma::mat arma_mu;      // source
  const arma::cube arma_sigma;
  const int n_rows_mu;
  const int n_cols_mu;
  const double target;
  const double total;
  const arma::vec arma_lower;
  const arma::vec arma_upper;
  const arma::vec arma_ones;
  const arma::mat arma_diag;
  arma::mat& arma_weights;      // destination (pass by reference)

  // initialize with source and destination
  RollMinRisk(const arma::mat arma_mu, const arma::cube arma_sigma,
              const int n_rows_mu, const int n_cols_mu,
              const double target, const double total,
              const arma::vec arma_lower, const arma::vec arma_upper,
              const arma::vec arma_ones, const arma::mat arma_diag,
              arma::mat& arma_weights)
    : arma_mu(arma_mu), arma_sigma(arma_sigma),
      n_rows_mu(n_rows_mu), n_cols_mu(n_cols_mu),
      target(target), total(total),
      arma_lower(arma_lower), arma_upper(arma_upper),
      arma_ones(arma_ones), arma_diag(arma_diag),
      arma_weights(arma_weights) { }

  // function call operator that iterates by slice
  void operator()(std::size_t begin_slice, std::size_t end_slice) {
    for (std::size_t i = begin_slice; i < end_slice; i++) {
      
      int n_solve = 0;
      int n_size = n_cols_mu;
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
      
      for (int j = 0; j < pow((long double)2.0, (long double)2.0 * n_cols_mu + 1); j++) {
        
        n_size = n_cols_mu + 1;
        arma::uvec arma_ix = arma::linspace<arma::uvec>(0, n_size - 1, n_size);
        
        for (int k = 0; k < 2 * n_cols_mu + 1; k++) {
          
          if ((j & (int)pow((long double)2.0, k))) {
            
            n_size += 1;
            
            arma_ix.resize(n_size);
            arma_ix(n_size - 1) = n_cols_mu + 1 + k;
            
          }
          
        }
        
        arma::mat A_subset = A(arma_ix, arma_ix);
        arma::vec b_subset = b(arma_ix);
        
        bool status_solve = arma::solve(weights, A_subset, b_subset,
                                        arma::solve_opts::no_approx);
        
        if (status_solve) {
          
          arma_ix = arma::linspace<arma::uvec>(0, n_cols_mu - 1, n_cols_mu);
          
          // weights
          arma::vec weights_subset = weights(arma_ix);
          arma::mat trans_weights = trans(weights_subset);
          
          if (all(weights_subset >= arma_lower[0]) && all(weights_subset <= arma_upper[0]) &&
              (sum(mu * weights_subset) >= target)) { // CHANGE BACK TO >=!
            
            n_solve += 1;
            
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
// [[Rcpp::export(.roll_min_risk)]]
NumericMatrix roll_min_risk(const NumericMatrix& mu, const NumericVector& sigma,
                            const double& target, const double& total,
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
  
  RollMinRisk roll_min_risk(arma_mu, arma_sigma, n_rows_mu, n_cols_mu,
                            target, total, arma_lower, arma_upper,
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

roll_result <- .roll_min_risk(mu, sigma, target, total, lower, upper)

i <- 5
roll_temp <- as.numeric(.roll_min_risk(mu[i, ], sigma[ , , i], target, total, lower, upper))
sum(mu[i, ] * roll_temp)
sqrt(t(roll_temp) %*% sigma[ , , i] %*% roll_temp)
*/

/*** R
rollapplyr_cube <- function(f, mu, sigma, target, total = 1, lower = 0, upper = 1) {
  
  mu_attr <- attributes(mu)
  
  n_rows_mu <- nrow(mu)
  n_cols_mu <- ncol(mu)
  result <- matrix(as.numeric(NA), n_rows_mu, n_cols_mu)
  
  for (i in 2:n_rows_mu) { # CHANGE TO 1?
    
    result[i, ] <- f(mu[i, ], sigma[ , , i], target, total, lower, upper)
    
  }
  
  attributes(result) <- mu_attr
  
  return(result)
  
}

min_risk_optim_cvxr <- function(mu, sigma, target, total = 1, lower = 0, upper = 1) {
  
  mu <- as.numeric(mu)
  
  params <- Variable(length(mu))
  
  cons <- list(params >= lower, params <= upper,
               sum(params) == total, sum(mu * params) >= target)
  
  obj <- Minimize(quad_form(params, sigma))
  
  result <- solve(Problem(obj, cons))$getValue(params)
  
  if (anyNA(result)) {
    result <- rep(NA, length(mu))
  }
  
  return(as.numeric(result))
  
}

cvxr_result <- rollapplyr_cube(min_risk_optim_cvxr, mu, sigma, target, total, lower, upper)

i <- 5
cvxr_temp <- min_risk_optim_cvxr(mu[i, ], sigma[ , , i], target, total, lower, upper)
sum(mu[i, ] * cvxr_temp)
sqrt(t(cvxr_temp) %*% sigma[ , , i] %*% cvxr_temp)
*/

/*** R
rollapplyr_min_risk <- function(object, expected, mu, sigma) {
  
  n_rows <- nrow(object)
  result_ls <- list()
  
  for (i in 1:n_rows) {
    
    weights_object <- as.numeric(object[i, ])
    weights_expected <- as.numeric(expected[i, ])
    
    temp <- NULL
    temp <- tryCatch(expect_equal(weights_object, weights_expected, ignore_attr = TRUE),
                     
                     error = function(x) {
                       
                       mu_object <- sum(mu[i, ] * weights_object)
                       mu_expected <- sum(mu[i, ] * weights_expected)
                       
                       sigma_object <- t(weights_object) %*% sigma[ , , i] %*% weights_object
                       sigma_expected <- t(weights_expected) %*% sigma[ , , i] %*% weights_expected 
                       
                       if (all(weights_expected >= 0) &&
                           ((mu_object < mu_expected) || (sigma_object > sigma_expected))) {
                         paste0("Error: not optimal (see row ", i, ")")
                       }
                       
                     })
    
    if (is.character(temp)) {
      result_ls <- append(result_ls, list(temp))
    }
    
  }
  
  result <- do.call(rbind, result_ls)
  
  return(result)
  
}

rollapplyr_min_risk(roll_result, cvxr_result, mu, sigma)

# microbenchmark("roll" = .roll_min_risk(mu, sigma, target, total, lower, upper),
#                "cvxr" = rollapplyr_cube(min_risk_optim_cvxr, mu, sigma, target, total, lower, upper))

*/