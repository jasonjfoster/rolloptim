#ifndef ROLLOPTIM_H
#define ROLLOPTIM_H

#define ARMA_WARN_LEVEL 0

#include <RcppArmadillo.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;

namespace rolloptim {

// 'Worker' function for computing the rolling optimization
struct RollMinVarMuFALSE : public Worker {
  
  const arma::cube arma_sigma;  // source
  const int n_rows;
  const int n_cols;
  const double total;
  const arma::vec arma_lower;
  const arma::vec arma_upper;
  const arma::vec arma_ones;
  const arma::mat arma_diag;
  const arma::mat arma_A;
  const arma::vec arma_b;
  arma::mat& arma_weights;      // destination (pass by reference)
  
  // initialize with source and destination
  RollMinVarMuFALSE(const arma::cube arma_sigma, const int n_rows,
                    const int n_cols, const double total,
                    const arma::vec arma_lower, const arma::vec arma_upper,
                    const arma::vec arma_ones, const arma::mat arma_diag,
                    const arma::mat arma_A, const arma::vec arma_b,
                    arma::mat& arma_weights)
    : arma_sigma(arma_sigma), n_rows(n_rows),
      n_cols(n_cols), total(total),
      arma_lower(arma_lower), arma_upper(arma_upper),
      arma_ones(arma_ones), arma_diag(arma_diag),
      arma_A(arma_A), arma_b(arma_b),
      arma_weights(arma_weights) { }
  
  // function call operator that iterates by slice
  void operator()(std::size_t begin_slice, std::size_t end_slice) {
    for (std::size_t i = begin_slice; i < end_slice; i++) {
      
      int n_solve = 0;
      int n_combn = 1 << (2 * n_cols); // 2 ^ (2 * n_cols)
      long double obj_prev = arma::datum::inf;
      arma::mat sigma = arma_sigma.slice(i);
      arma::mat A(arma_A.begin(), n_cols + 1 + n_cols + n_cols,
                  n_cols + 1 + n_cols + n_cols);
      arma::vec b(arma_b.begin(), n_cols + 1 + n_cols + n_cols);
      
      A.submat(0, 0, n_cols - 1, n_cols - 1) = sigma;
      
      // number of binary combinations
      for (int j = 0; j < n_combn; j++) {
        
        arma::uvec arma_ix(n_cols + 2 * n_cols + 1);
        arma_ix.subvec(0, n_cols) = arma::linspace<arma::uvec>(1, n_cols + 1, n_cols + 1);
        
        // find the binary combination for lower and upper bounds
        for (int k = 0; k < 2 * n_cols; k++) {
          if ((j & (1 << k)) == 0) { 
            arma_ix[n_cols + 1 + k] = n_cols + 1 + k;
          }          
        }
        
        arma::uvec arma_ix_subset = find(arma_ix);
        arma::mat A_subset = A(arma_ix_subset, arma_ix_subset);
        arma::vec b_subset = b(arma_ix_subset);
        arma::vec weights(n_cols);

        // check if solution is found 
        bool status_solve = arma::solve(weights, A_subset, b_subset,
                                        arma::solve_opts::no_approx);
        
        // don't find approximate solution for rank deficient system
        if (status_solve) {
                    
          // weights
          arma::vec weights_subset = weights.subvec(0, n_cols - 1);
          
          // check if constraints are satisfied
          if ((weights_subset.min() - arma_lower[0] >= -arma::datum::eps) &&
              (weights_subset.max() - arma_upper[0] <= arma::datum::eps)) {
            
            n_solve += 1;
            
            // objective value
            long double obj = arma::as_scalar(trans(weights_subset) * sigma * weights_subset);
            
            if (obj <= obj_prev) {
              
              obj_prev = obj;
              
              arma_weights.row(i) = trans(weights_subset);
              
            }
            
          }
          
        }
        
      }
      
      if (n_solve == 0) {
        arma_weights.row(i).fill(NA_REAL);
      }
      
    }
  }
  
};

// 'Worker' function for computing the rolling optimization
struct RollMinVarMuTRUE : public Worker {
  
  const arma::cube arma_sigma;  // source
  const arma::mat arma_mu;
  const int n_rows;
  const int n_cols;
  const double total;
  const arma::vec arma_target;
  const arma::vec arma_lower;
  const arma::vec arma_upper;
  const arma::vec arma_ones;
  const arma::mat arma_diag;
  const arma::mat arma_A;
  const arma::vec arma_b;
  arma::mat& arma_weights;      // destination (pass by reference)
  
  // initialize with source and destination
  RollMinVarMuTRUE(const arma::cube arma_sigma, const arma::mat arma_mu, 
                   const int n_rows, const int n_cols,
                   const double total, const arma::vec arma_target,
                   const arma::vec arma_lower, const arma::vec arma_upper,
                   const arma::vec arma_ones, const arma::mat arma_diag,
                   const arma::mat arma_A, const arma::vec arma_b, 
                   arma::mat& arma_weights)
    : arma_sigma(arma_sigma), arma_mu(arma_mu),
      n_rows(n_rows), n_cols(n_cols),
      total(total), arma_target(arma_target),
      arma_lower(arma_lower), arma_upper(arma_upper),
      arma_ones(arma_ones), arma_diag(arma_diag),
      arma_A(arma_A), arma_b(arma_b),
      arma_weights(arma_weights) { }
  
  // function call operator that iterates by slice
  void operator()(std::size_t begin_slice, std::size_t end_slice) {
    for (std::size_t i = begin_slice; i < end_slice; i++) {
      
      int n_solve = 0;
      int n_combn = 1 << (2 * n_cols + 1); // 2 ^ (2 * n_cols + 1)
      long double obj_prev = arma::datum::inf;
      arma::mat sigma = arma_sigma.slice(i);
      arma::mat mu = arma_mu.row(i);
      arma::mat A(arma_A.begin(), n_cols + 1 + 1 + n_cols + n_cols,
                  n_cols + 1 + 1 + n_cols + n_cols);
      arma::vec b(arma_b.begin(), n_cols + 1 + 1 + n_cols + n_cols);
      
      A.submat(0, 0, n_cols - 1, n_cols - 1) = sigma;
      A.submat(n_cols + 1, 0, n_cols + 1, n_cols - 1) = mu;
      A.submat(0, n_cols + 1, n_cols - 1, n_cols + 1) = trans(mu);
      b(n_cols + 1) = arma_target(i);

      // number of binary combinations
      for (int j = 0; j < n_combn; j++) {
        
        arma::uvec arma_ix(n_cols + 2 * n_cols + 2);
        arma_ix.subvec(0, n_cols) = arma::linspace<arma::uvec>(1, n_cols + 1, n_cols + 1);

        // mu inequality
        if ((j & 1) == 0) {
          arma_ix[n_cols + 1] = n_cols + 1;
        }

        // find the binary combination for lower and upper bounds
        for (int k = 0; k < 2 * n_cols; k++) {
          if ((j & (1 << (k + 1))) == 0) { 
            arma_ix[n_cols + 2 + k] = n_cols + 2 + k;
          }          
        }

        arma::uvec arma_ix_subset = find(arma_ix);
        arma::mat A_subset = A(arma_ix_subset, arma_ix_subset);
        arma::vec b_subset = b(arma_ix_subset);
        arma::vec weights(n_cols);

        // check if solution is found 
        bool status_solve = arma::solve(weights, A_subset, b_subset,
                                        arma::solve_opts::no_approx);
        
        // don't find approximate solution for rank deficient system
        if (status_solve) {
          
          arma::vec weights_subset = weights.subvec(0, n_cols - 1);
          double total_mu = arma::as_scalar(mu * weights_subset);

          // check if constraints are satisfied
          if ((weights_subset.min() - arma_lower[0] >= -arma::datum::eps) &&
              (weights_subset.max() - arma_upper[0] <= arma::datum::eps) &&
              (arma_target(i) - total_mu <= arma::datum::eps)) {
            
            n_solve += 1;
            
            // objective value
            long double obj = arma::as_scalar(trans(weights_subset) * sigma * weights_subset);
            
            if (obj <= obj_prev) {
              
              obj_prev = obj;
              
              arma_weights.row(i) = trans(weights_subset);
              
            }
            
          }
          
        }
        
      }
      
      if (n_solve == 0) {
        arma_weights.row(i).fill(NA_REAL);
      }
      
    }
  }
  
};

// 'Worker' function for computing the rolling optimization
struct RollMaxMeanSigmaFALSE : public Worker {
  
  const arma::mat arma_mu;      // source
  const int n_rows;
  const int n_cols;
  const double total;
  const arma::vec arma_lower;
  const arma::vec arma_upper;
  const arma::vec arma_ones;
  const arma::mat arma_diag;
  const arma::mat arma_A;
  const arma::vec arma_b;
  arma::mat& arma_weights;      // destination (pass by reference)
  
  // initialize with source and destination
  RollMaxMeanSigmaFALSE(const arma::mat arma_mu, const int n_rows,
                        const int n_cols, const double total,
                        const arma::vec arma_lower, const arma::vec arma_upper,
                        const arma::vec arma_ones, const arma::mat arma_diag,
                        const arma::mat arma_A, const arma::vec arma_b,
                        arma::mat& arma_weights)
    : arma_mu(arma_mu), n_rows(n_rows),
      n_cols(n_cols), total(total),
      arma_lower(arma_lower), arma_upper(arma_upper),
      arma_ones(arma_ones), arma_diag(arma_diag),
      arma_A(arma_A), arma_b(arma_b),
      arma_weights(arma_weights) { }
  
  // function call operator that iterates by slice
  void operator()(std::size_t begin_slice, std::size_t end_slice) {
    for (std::size_t i = begin_slice; i < end_slice; i++) {
      
      int n_solve = 0;
      int n_combn = 1 << (2 * n_cols); // 2 ^ (2 * n_cols)
      long double obj_prev = arma::datum::inf;
      arma::mat mu = arma_mu.row(i);
      arma::mat A(arma_A.begin(), n_cols + 1 + n_cols + n_cols,
                  n_cols + 1 + n_cols + n_cols);
      arma::vec b(arma_b.begin(), n_cols + 1 + n_cols + n_cols);
      
      A.submat(0, 0, 0, n_cols - 1) = mu;
      A.submat(0, 0, n_cols - 1, 0) = trans(mu);
      
      // number of binary combinations
      for (int j = 0; j < n_combn; j++) {
        
        arma::uvec arma_ix(n_cols + 2 * n_cols + 1);
        arma_ix.subvec(0, n_cols) = arma::linspace<arma::uvec>(1, n_cols + 1, n_cols + 1);
        
        // find the binary combination for lower and upper bounds
        for (int k = 0; k < 2 * n_cols; k++) {
          if ((j & (1 << k)) == 0) {
            arma_ix[n_cols + 1 + k] = n_cols + 1 + k;
          }          
        }
        
        arma::uvec arma_ix_subset = find(arma_ix);
        arma::mat A_subset = A(arma_ix_subset, arma_ix_subset);
        arma::vec b_subset = b(arma_ix_subset);
        arma::vec weights(n_cols);

        // check if solution is found 
        bool status_solve = arma::solve(weights, A_subset, b_subset,
                                        arma::solve_opts::no_approx);
        
        // don't find approximate solution for rank deficient system
        if (status_solve) {
                    
          // weights
          arma::vec weights_subset = weights.subvec(0, n_cols - 1);
          
          // check if constraints are satisfied
          if ((weights_subset.min() - arma_lower[0] >= -arma::datum::eps) &&
              (weights_subset.max() - arma_upper[0] <= arma::datum::eps)) {
            
            n_solve += 1;
            
            // objective value
            long double obj = arma::as_scalar(-mu * weights_subset);
            
            if (obj <= obj_prev) {
              
              obj_prev = obj;
              
              arma_weights.row(i) = trans(weights_subset);
              
            }
            
          }
          
        }
        
      }
      
      if (n_solve == 0) {
        arma_weights.row(i).fill(NA_REAL);
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
      
      int n_solve = 0;
      int n_combn = 1 << (2 * n_cols_mu); // 2 ^ (2 * n_cols_mu)
      long double obj_prev = arma::datum::inf;
      arma::mat mu = arma_mu.row(i);
      arma::mat sigma = arma_sigma.slice(i);
      arma::mat A(arma_A.begin(), n_cols_mu + 1 + n_cols_mu + n_cols_mu,
                  n_cols_mu + 1 + n_cols_mu + n_cols_mu);
      arma::vec b(arma_b.begin(), n_cols_mu + 1 + n_cols_mu + n_cols_mu);
      
      A.submat(0, 0, n_cols_mu - 1, n_cols_mu - 1) = lambda * sigma;
      b.subvec(0, n_cols_mu - 1) = trans(mu);
      
      // number of binary combinations
      for (int j = 0; j < n_combn; j++) {
        
        arma::uvec arma_ix(n_cols_mu + 1 + 2 * n_cols_mu);
        arma_ix.subvec(0, n_cols_mu) = arma::linspace<arma::uvec>(1, n_cols_mu + 1, n_cols_mu + 1);
        
        // find the binary combination
        for (int k = 0; k < 2 * n_cols_mu; k++) {
          if ((j & (1 << k)) == 0) { 
            arma_ix[n_cols_mu + 1 + k] = n_cols_mu + 1 + k;
          }          
        }
        
        arma::uvec arma_ix_subset = find(arma_ix);
        arma::mat A_subset = A(arma_ix_subset, arma_ix_subset);
        arma::vec b_subset = b(arma_ix_subset);
        arma::vec weights(n_cols_mu);

        // check if solution is found 
        bool status_solve = arma::solve(weights, A_subset, b_subset,
                                        arma::solve_opts::no_approx);
        
        // don't find approximate solution for rank deficient system
        if (status_solve) {
                    
          // weights
          arma::vec weights_subset = weights.subvec(0, n_cols_mu - 1);
          
          // check if constraints are satisfied
          if ((weights_subset.min() - arma_lower[0] >= -arma::datum::eps) &&
              (weights_subset.max() - arma_upper[0] <= arma::datum::eps)) {
            
            n_solve += 1;
            
            // objective value
            long double obj = 0.5 * lambda * arma::as_scalar(trans(weights_subset) * sigma * weights_subset) -
              arma::dot(mu, weights_subset);
            
            if (obj <= obj_prev) {
              
              obj_prev = obj;
              
              arma_weights.row(i) = trans(weights_subset);
              
            }
            
          }
          
        }
        
      }
      
      if (n_solve == 0) {        
        arma_weights.row(i).fill(NA_REAL);
      }
      
    }
  }
  
};

}

#endif