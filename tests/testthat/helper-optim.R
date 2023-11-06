rollapplyr_optim <- function(f, mu = NULL, sigma = NULL) {
  
  if (!is.null(mu) && !is.null(sigma)) {

    mu_attr <- attributes(mu)

    n_rows_mu <- nrow(mu)
    n_cols_mu <- ncol(mu)
    result <- matrix(as.numeric(NA), n_rows_mu, n_cols_mu)

    for (i in 1:n_rows_mu) {

      status_solve <- tryCatch(f(mu[i, ], sigma[ , , i]),

                               error = function(x) {
                                 rep(NA, n_cols_mu)
                               })

      result[i, ] <- status_solve

    }

    attributes(result) <- mu_attr

  } else if (!is.null(mu)) {

    mu_attr <- attributes(mu)

    n_rows_mu <- nrow(mu)
    n_cols_mu <- ncol(mu)
    result <- matrix(as.numeric(NA), n_rows_mu, n_cols_mu)

    for (i in 1:n_rows_mu) {

      status_solve <- tryCatch(f(mu[i, ]),

                               error = function(x) {
                                 rep(NA, n_cols_mu)
                               })

      result[i, ] <- status_solve

    }

    attributes(result) <- mu_attr

  } else if (!is.null(sigma)) {
    
    sigma_dimnames <- dimnames(sigma)
    
    dim_sigma <- dim(sigma)
    n_rows_sigma <- dim_sigma[3]
    n_cols_sigma <- dim_sigma[2]
    result <- matrix(as.numeric(NA), n_rows_sigma, n_cols_sigma)
    
    for (i in 1:n_rows_sigma) {
      
      status_solve <- tryCatch(f(sigma[ , , i]),

                               error = function(x) {
                                 rep(NA, n_cols_sigma)
                               })
      
      result[i, ] <- status_solve
      
    }
    
    if (length(sigma_dimnames) > 1) {
      attr(result, "dimnames") <- list(NULL, sigma_dimnames[[1]])
    } 
    
  }
  
  return(result)
  
}

rollapplyr_xy <- function(f, x, y, width) {
  
  n_rows_xy <- nrow(x)
  n_cols_x <- ncol(x)
  result <- matrix(as.numeric(NA), n_rows_xy, n_cols_x)
  
  for (i in 1:n_rows_xy) {
    
    x_subset <- x[max(1, i - width + 1):i, , drop = FALSE]
    y_subset <- y[max(1, i - width + 1):i, , drop = FALSE]
    
    result[i, ] <- f(x_subset, y_subset)
    
  }
  
  return(result)
  
}