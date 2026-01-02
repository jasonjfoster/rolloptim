rollapplyr_optim <- function(f, mu = NULL, sigma = NULL, target = NULL) {
  
  if (!is.null(mu) && !is.null(sigma) && !is.null(target)) {

    mu_attr <- attributes(mu)

    n_rows <- nrow(mu)
    n_cols <- ncol(mu)
    result <- matrix(as.numeric(NA), n_rows, n_cols)

    for (i in 1:n_rows) {

      status_solve <- tryCatch(f(sigma = sigma[ , , i], mu = mu[i, ],
                                 target = target[i]),

                               error = function(x) {
                                 rep(NA, n_cols)
                               })

      result[i, ] <- status_solve

    }

    attributes(result) <- mu_attr

  } else if (!is.null(mu) && !is.null(sigma)) {

    mu_attr <- attributes(mu)

    n_rows <- nrow(mu)
    n_cols <- ncol(mu)
    result <- matrix(as.numeric(NA), n_rows, n_cols)

    for (i in 1:n_rows) {

      status_solve <- tryCatch(f(sigma = sigma[ , , i], mu = mu[i, ]),

                               error = function(x) {
                                 rep(NA, n_cols)
                               })

      result[i, ] <- status_solve

    }

    attributes(result) <- mu_attr

  } else if (!is.null(mu)) {

    mu_attr <- attributes(mu)

    n_rows <- nrow(mu)
    n_cols <- ncol(mu)
    result <- matrix(as.numeric(NA), n_rows, n_cols)

    for (i in 1:n_rows) {

      status_solve <- tryCatch(f(mu = mu[i, ]),

                               error = function(x) {
                                 rep(NA, n_cols)
                               })

      result[i, ] <- status_solve

    }

    attributes(result) <- mu_attr

  } else if (!is.null(sigma)) {
    
    sigma_dimnames <- dimnames(sigma)
    
    dim_sigma <- dim(sigma)
    n_rows <- dim_sigma[3]
    n_cols <- dim_sigma[2]
    result <- matrix(as.numeric(NA), n_rows, n_cols)
    
    for (i in 1:n_rows) {
      
      status_solve <- tryCatch(f(sigma = sigma[ , , i]),

                               error = function(x) {
                                 rep(NA, n_cols)
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
  
  n_rows <- if (is.matrix(x)) nrow(x) else length(x)
  n_cols <- if (is.matrix(x)) ncol(x) else 1
  result <- matrix(as.numeric(NA), n_rows, n_cols)
  
  for (i in 1:n_rows) {
    
    x_subset <- x[max(1, i - width + 1):i, , drop = FALSE]
    y_subset <- y[max(1, i - width + 1):i, , drop = FALSE]
    
    result[i, ] <- f(x_subset, y_subset)
    
  }
  
  return(result)
  
}