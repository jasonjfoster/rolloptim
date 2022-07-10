#ifndef ROLLPORT_H
#define ROLLPORT_H

#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;

// 'Worker' function for computing rolling portfolio optimization
struct RollMinRisk : public Worker {
  
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
  RollMinRisk(const arma::mat arma_mu, const arma::cube arma_sigma,
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

// 'Worker' function for computing rolling portfolio optimization
struct RollMaxReturn : public Worker {
  
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
  RollMaxReturn(const arma::mat arma_mu, const arma::cube arma_sigma,
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

// 'Worker' function for computing rolling portfolio optimization
struct RollMaxRatio : public Worker {
  
  const arma::mat arma_mu;      // source
  const arma::cube arma_sigma;
  const int n_rows_mu;
  const int n_cols_mu;
  const double gamma;
  const double total;
  const arma::vec arma_lower;
  const arma::vec arma_upper;
  const arma::vec arma_ones;
  const arma::mat arma_diag;
  arma::mat& arma_weights;      // destination (pass by reference)
  
  // initialize with source and destination
  RollMaxRatio(const arma::mat arma_mu, const arma::cube arma_sigma,
               const int n_rows_mu, const int n_cols_mu,
               const double gamma, const double total,
               const arma::vec arma_lower, const arma::vec arma_upper,
               const arma::vec arma_ones, const arma::mat arma_diag,
               arma::mat& arma_weights)
    : arma_mu(arma_mu), arma_sigma(arma_sigma),
      n_rows_mu(n_rows_mu), n_cols_mu(n_cols_mu),
      gamma(gamma), total(total),
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
      arma::mat A = gamma * sigma;
      arma::vec b = trans(mu);
      arma::vec weights(n_size);
      
      // total constraint
      n_size += 1;
      
      A = join_rows(A, arma_ones);
      
      b.resize(n_size);
      b(n_size - 1) = total;
      
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
        
        n_size = n_cols_mu + 1;
        arma::uvec arma_ix = arma::linspace<arma::uvec>(0, n_size - 1, n_size);
        
        // find the index combination
        for (int k = 0; k < 2 * n_cols_mu; k++) {
          
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
          if (all(weights_subset >= arma_lower[0]) && all(weights_subset <= arma_upper[0])) {
            
            n_solve += 1;
            
            // objective value
            obj = gamma * as_scalar(trans_weights * sigma * weights_subset) -
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

#endif