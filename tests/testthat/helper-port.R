rollapplyr_port <- function(f, mu, sigma, ...) {
  
  mu_attr <- attributes(mu)
  
  n_rows_mu <- nrow(mu)
  n_cols_mu <- ncol(mu)
  result <- matrix(as.numeric(NA), n_rows_mu, n_cols_mu)
  
  for (i in 1:n_rows_mu) {
    
    status_solve <- tryCatch(f(mu[i, ], sigma[ , , i], ...),
                             
                             error = function(x) {
                               rep(NA, n_cols_mu)
                             })
    
    result[i, ] <- status_solve
    
  }
  
  attributes(result) <- mu_attr
  
  return(result)
  
}