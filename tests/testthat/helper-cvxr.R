cvxr_min_var <- function(sigma, mu = NULL, target = NULL,
                         total = 1, lower = 0, upper = 1) {
  
  n_cols <- dim(sigma)[1]
  
  w <- CVXR::Variable(n_cols)
  
  objective <- CVXR::Minimize(CVXR::quad_form(w, sigma))
  constraints <- list(sum(w) == total, w >= lower, w <= upper)
  
  if (!is.null(mu) && !is.null(target)) {
    constraints <- c(constraints, list(t(as.numeric(mu)) %*% w >= target))
  }

  problem <- CVXR::Problem(objective, constraints)
  
  result <- CVXR::solve(problem)$getValue(w)
  
  return(result)
  
}

cvxr_max_mean <- function(mu, sigma = NULL, target = NULL,
                          total = 1, lower = 0, upper = 1) {
  
  mu <- zoo::coredata(mu)
  
  n_cols <- ncol(mu)
  
  w <- CVXR::Variable(n_cols)
  
  objective <- CVXR::Maximize(mu %*% w)
  constraints <- list(sum(w) == total, w >= lower, w <= upper)

  if (!is.null(sigma) && !is.null(target)) {
    constraints <- c(constraints, list(CVXR::quad_form(w, sigma) <= target))
  }

  problem <- CVXR::Problem(objective, constraints)
  
  result <- CVXR::solve(problem)$getValue(w)
  
  return(result)
  
}

cvxr_max_utility <- function(mu, sigma, lambda = 1, total = 1, lower = 0, upper = 1) {
  
  mu <- zoo::coredata(mu)
  
  n_cols <- ncol(mu)
  
  w <- CVXR::Variable(n_cols)
  
  objective <- CVXR::Minimize(0.5 * lambda * CVXR::quad_form(w, sigma) - mu %*% w)
  constraints <- list(sum(w) == total, w >= lower, w <= upper)
  
  problem <- CVXR::Problem(objective, constraints)
  
  result <- CVXR::solve(problem)$getValue(w)
  
  return(result)
  
}

cvxr_min_rss <- function(x, y, total = 1, lower = 0, upper = 1) {
  
  x <- zoo::coredata(x)
  y <- zoo::coredata(y)
  
  n_cols <- ncol(x)
  
  w <- CVXR::Variable(n_cols)
  
  objective <- CVXR::Minimize(CVXR::sum_squares(y - x %*% w))
  constraints <- list(sum(w) == total, w >= lower, w <= upper)
  
  problem <- CVXR::Problem(objective, constraints)
  
  result <- CVXR::solve(problem)$getValue(w)
  
  return(result)
  
}
