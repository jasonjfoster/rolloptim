#ifndef ROLLPORT_H
#define ROLLPORT_H

#define ARMA_WARN_LEVEL 0

#include <RcppArmadillo.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;

namespace rollport {

// 'Worker' function for computing the rolling optimization
struct RollMinVar : public Worker {
  
  const arma::cube arma_sigma;  // source
  const int n_rows_sigma;
  const int n_cols_sigma;
  const double total;
  const arma::vec arma_lower;
  const arma::vec arma_upper;
  const arma::vec arma_ones;
  const arma::mat arma_diag;
  const arma::mat arma_A;
  const arma::vec arma_b;
  arma::mat& arma_weights;      // destination (pass by reference)
  
  // initialize with source and destination
  RollMinVar(const arma::cube arma_sigma, const int n_rows_sigma,
             const int n_cols_sigma, const double total,
             const arma::vec arma_lower, const arma::vec arma_upper,
             const arma::vec arma_ones, const arma::mat arma_diag,
             const arma::mat arma_A, const arma::vec arma_b,
             arma::mat& arma_weights)
    : arma_sigma(arma_sigma), n_rows_sigma(n_rows_sigma),
      n_cols_sigma(n_cols_sigma), total(total),
      arma_lower(arma_lower), arma_upper(arma_upper),
      arma_ones(arma_ones), arma_diag(arma_diag),
      arma_A(arma_A), arma_b(arma_b),
      arma_weights(arma_weights) { }
  
  // function call operator that iterates by slice
  void operator()(std::size_t begin_slice, std::size_t end_slice) {
    for (std::size_t i = begin_slice; i < end_slice; i++) {
      
      int n_size = 0;
      int n_solve = 0;
      int n_combn = pow((long double)2.0, (long double)2.0 * n_cols_sigma);
      long double obj = arma::datum::inf;
      long double obj_prev = arma::datum::inf;
      arma::mat sigma = arma_sigma.slice(i);
      arma::mat A(arma_A.begin(), n_cols_sigma + 1 + n_cols_sigma + n_cols_sigma,
                  n_cols_sigma + 1 + n_cols_sigma + n_cols_sigma);
      arma::vec b(arma_b.begin(), n_cols_sigma + 1 + n_cols_sigma + n_cols_sigma);
      arma::vec weights(n_cols_sigma);
      
      A.submat(0, 0, n_cols_sigma - 1, n_cols_sigma - 1) = sigma;
      
      // number of binary combinations
      for (int j = 0; j < n_combn; j++) {
        
        n_size = j;
        arma::uvec arma_ix(n_cols_sigma + 2 * n_cols_sigma + 1);
        arma_ix.subvec(0, n_cols_sigma) = arma::linspace<arma::uvec>(1, n_cols_sigma + 1, n_cols_sigma + 1);
        
        // find the binary combination
        for (int k = 0; k < 2 * n_cols_sigma; k++) {
          
          if (n_size % 2 == 0) {
            arma_ix[n_cols_sigma + 1 + k] = n_cols_sigma + 1 + k;
          }
          
          n_size /= 2;
          
        }
        
        arma::uvec arma_ix_subset = find(arma_ix);
        arma::mat A_subset = A(arma_ix_subset, arma_ix_subset);
        arma::vec b_subset = b(arma_ix_subset);
        
        // check if solution is found 
        bool status_solve = arma::solve(weights, A_subset, b_subset,
                                        arma::solve_opts::no_approx);
        
        // don't find approximate solution for rank deficient system
        if (status_solve) {
          
          arma_ix = arma::linspace<arma::uvec>(0, n_cols_sigma - 1, n_cols_sigma);
          
          // weights
          arma::vec weights_subset = weights(arma_ix);
          arma::mat trans_weights = trans(weights_subset);
          
          // check if constraints are satisfied
          if ((weights_subset.min() - arma_lower[0] >= -sqrt(arma::datum::eps)) &&
              (weights_subset.max() - arma_upper[0] <= sqrt(arma::datum::eps))) {
            
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
        
        arma::rowvec no_solution(n_cols_sigma);
        no_solution.fill(NA_REAL);
        
        arma_weights.row(i) = no_solution;
        
      }
      
    }
  }
  
};

// 'Worker' function for computing the rolling optimization
struct RollMaxMean : public Worker {
  
  const arma::mat arma_mu;      // source
  const int n_rows_mu;
  const int n_cols_mu;
  const double total;
  const arma::vec arma_lower;
  const arma::vec arma_upper;
  const arma::vec arma_ones;
  const arma::mat arma_diag;
  const arma::mat arma_A;
  const arma::vec arma_b;
  arma::mat& arma_weights;      // destination (pass by reference)
  
  // initialize with source and destination
  RollMaxMean(const arma::mat arma_mu, const int n_rows_mu,
              const int n_cols_mu, const double total,
              const arma::vec arma_lower, const arma::vec arma_upper,
              const arma::vec arma_ones, const arma::mat arma_diag,
              const arma::mat arma_A, const arma::vec arma_b,
              arma::mat& arma_weights)
    : arma_mu(arma_mu), n_rows_mu(n_rows_mu),
      n_cols_mu(n_cols_mu), total(total),
      arma_lower(arma_lower), arma_upper(arma_upper),
      arma_ones(arma_ones), arma_diag(arma_diag),
      arma_A(arma_A), arma_b(arma_b),
      arma_weights(arma_weights) { }
  
  // function call operator that iterates by slice
  void operator()(std::size_t begin_slice, std::size_t end_slice) {
    for (std::size_t i = begin_slice; i < end_slice; i++) {
      
      int n_size = 0;
      int n_solve = 0;
      int n_combn = pow((long double)2.0, (long double)2.0 * n_cols_mu);
      long double obj = arma::datum::inf;
      long double obj_prev = arma::datum::inf;
      arma::mat mu = arma_mu.row(i);
      arma::mat A(arma_A.begin(), n_cols_mu + 1 + n_cols_mu + n_cols_mu,
                  n_cols_mu + 1 + n_cols_mu + n_cols_mu);
      arma::vec b(arma_b.begin(), n_cols_mu + 1 + n_cols_mu + n_cols_mu);
      arma::vec weights(n_cols_mu);
      
      A.submat(0, 0, 0, n_cols_mu - 1) = mu;
      A.submat(0, 0, n_cols_mu - 1, 0) = trans(mu);
      
      // number of binary combinations
      for (int j = 0; j < n_combn; j++) {
        
        n_size = j;
        arma::uvec arma_ix(n_cols_mu + 2 * n_cols_mu + 1);
        arma_ix.subvec(0, n_cols_mu) = arma::linspace<arma::uvec>(1, n_cols_mu + 1, n_cols_mu + 1);
        
        // find the binary combination
        for (int k = 0; k < 2 * n_cols_mu; k++) {
          
          if (n_size % 2 == 0) {
            arma_ix[n_cols_mu + 1 + k] = n_cols_mu + 1 + k;
          }
          
          n_size /= 2;
          
        }
        
        arma::uvec arma_ix_subset = find(arma_ix);
        arma::mat A_subset = A(arma_ix_subset, arma_ix_subset);
        arma::vec b_subset = b(arma_ix_subset);
        
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
          if ((weights_subset.min() - arma_lower[0] >= -sqrt(arma::datum::eps)) &&
              (weights_subset.max() - arma_upper[0] <= sqrt(arma::datum::eps))) {
            
            n_solve += 1;
            
            // objective value
            obj = as_scalar(-mu * weights_subset);
            
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
struct RollMaxUtility : public Worker {
  
  const arma::mat arma_mu;      // source
  const arma::cube arma_sigma;
  const int n_rows_mu;
  const int n_cols_mu;
  const double lambda;
  const double total;
  const arma::vec arma_lower;
  const arma::vec arma_upper;
  const arma::vec arma_ones;
  const arma::mat arma_diag;
  const arma::mat arma_A;
  const arma::vec arma_b;
  arma::mat& arma_weights;      // destination (pass by reference)
  
  // initialize with source and destination
  RollMaxUtility(const arma::mat arma_mu, const arma::cube arma_sigma,
                 const int n_rows_mu, const int n_cols_mu,
                 const double lambda, const double total,
                 const arma::vec arma_lower, const arma::vec arma_upper,
                 const arma::vec arma_ones, const arma::mat arma_diag,
                 const arma::mat arma_A, const arma::vec arma_b,
                 arma::mat& arma_weights)
    : arma_mu(arma_mu), arma_sigma(arma_sigma),
      n_rows_mu(n_rows_mu), n_cols_mu(n_cols_mu),
      lambda(lambda), total(total),
      arma_lower(arma_lower), arma_upper(arma_upper),
      arma_ones(arma_ones), arma_diag(arma_diag),
      arma_A(arma_A), arma_b(arma_b),
      arma_weights(arma_weights) { }
  
  // function call operator that iterates by slice
  void operator()(std::size_t begin_slice, std::size_t end_slice) {
    for (std::size_t i = begin_slice; i < end_slice; i++) {
      
      int n_size = 0;
      int n_solve = 0;
      int n_combn = pow((long double)2.0, (long double)2.0 * n_cols_mu);
      long double obj = arma::datum::inf;
      long double obj_prev = arma::datum::inf;
      arma::mat mu = arma_mu.row(i);
      arma::mat sigma = arma_sigma.slice(i);
      arma::mat A(arma_A.begin(), n_cols_mu + 1 + n_cols_mu + n_cols_mu,
                  n_cols_mu + 1 + n_cols_mu + n_cols_mu);
      arma::vec b(arma_b.begin(), n_cols_mu + 1 + n_cols_mu + n_cols_mu);
      arma::vec weights(n_cols_mu);
      
      A.submat(0, 0, n_cols_mu - 1, n_cols_mu - 1) = lambda * sigma;
      b.subvec(0, n_cols_mu - 1) = trans(mu);
      
      // number of binary combinations
      for (int j = 0; j < n_combn; j++) {
        
        n_size = j;
        arma::uvec arma_ix(n_cols_mu + 1 + 2 * n_cols_mu);
        arma_ix.subvec(0, n_cols_mu) = arma::linspace<arma::uvec>(1, n_cols_mu + 1, n_cols_mu + 1);
        
        // find the binary combination
        for (int k = 0; k < 2 * n_cols_mu; k++) {
          
          if (n_size % 2 == 0) {
            arma_ix[n_cols_mu + 1 + k] = n_cols_mu + 1 + k;
          }
          
          n_size /= 2;
          
        }
        
        arma::uvec arma_ix_subset = find(arma_ix);
        arma::mat A_subset = A(arma_ix_subset, arma_ix_subset);
        arma::vec b_subset = b(arma_ix_subset);
        
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
          if ((weights_subset.min() - arma_lower[0] >= -sqrt(arma::datum::eps)) &&
              (weights_subset.max() - arma_upper[0] <= sqrt(arma::datum::eps))) {
            
            n_solve += 1;
            
            // objective value
            obj = 0.5 * lambda * as_scalar(trans_weights * sigma * weights_subset) -
              sum(mu * weights_subset);
            
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

}

#endif