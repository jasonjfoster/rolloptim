rollapplyr_optim <- function(f, mu = NULL, sigma = NULL, target = NULL) {
  
  if (!is.null(mu) && !is.null(sigma) && !is.null(target)) {

    mu_attr <- attributes(mu)
    mu <- as.matrix(zoo::coredata(mu))

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
    mu <- as.matrix(zoo::coredata(mu))

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
    mu <- as.matrix(zoo::coredata(mu))

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

rollapplyr_xy <- function(f, x, y, width, min_obs) {
  
  x <- as.matrix(zoo::coredata(x))
  y <- as.matrix(zoo::coredata(y))
  
  n_rows <- nrow(x)
  n_cols <- ncol(x)
  result <- matrix(as.numeric(NA), n_rows, n_cols)
  
  for (i in 1:n_rows) {
    
    x_subset <- x[max(1, i - width + 1):i, , drop = FALSE]
    y_subset <- y[max(1, i - width + 1):i, , drop = FALSE]
    
    # roll::roll_crossprod: 'min_obs' complete observations for 'xx' and 'xy'
    n_obs_xx <- sum(complete.cases(x_subset))
    n_obs_xy <- sum(complete.cases(x_subset, y_subset))
    
    if ((n_obs_xx >= min_obs) && (n_obs_xy >= min_obs)) {
      
      status_solve <- tryCatch(f(x_subset, y_subset),
                               
                               error = function(x) {
                                 rep(NA, n_cols)
                               })
      
      result[i, ] <- status_solve
      
    }
    
  }
  
  return(result)
  
}
