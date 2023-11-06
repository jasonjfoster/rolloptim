cvxr_min_var <- function(sigma, total = 1, lower = 0, upper = 1) {
  
  n_cols_sigma <- dim(sigma)[1]
  
  w <- CVXR::Variable(n_cols_sigma)
  
  objective <- CVXR::Minimize(CVXR::quad_form(w, sigma))
  constraints <- list(sum(w) == total, w >= lower, w <= upper)
  
  problem <- CVXR::Problem(objective, constraints)
  
  result <- CVXR::solve(problem)$getValue(w)
  
  return(result)
  
}

cvxr_max_mean <- function(mu, total = 1, lower = 0, upper = 1) {
  
  mu <- zoo::coredata(mu)
  
  n_cols_mu <- ncol(mu)
  
  w <- CVXR::Variable(n_cols_mu)
  
  objective <- CVXR::Maximize(mu %*% w)
  constraints <- list(sum(w) == total, w >= lower, w <= upper)
  
  problem <- CVXR::Problem(objective, constraints)
  
  result <- CVXR::solve(problem)$getValue(w)
  
  return(result)
  
}

cvxr_max_utility <- function(mu, sigma, lambda = 1, total = 1, lower = 0, upper = 1) {
  
  mu <- zoo::coredata(mu)
  
  n_cols_mu <- ncol(mu)
  
  w <- CVXR::Variable(n_cols_mu)
  
  objective <- CVXR::Minimize(0.5 * lambda * CVXR::quad_form(w, sigma) - mu %*% w)
  constraints <- list(sum(w) == total, w >= lower, w <= upper)
  
  problem <- CVXR::Problem(objective, constraints)
  
  result <- CVXR::solve(problem)$getValue(w)
  
  return(result)
  
}

cvxr_min_rss <- function(x, y, total = 1, lower = 0, upper = 1) {
  
  x <- zoo::coredata(x)
  y <- zoo::coredata(y)
  
  n_cols_x <- ncol(x)
  
  w <- CVXR::Variable(n_cols_x)
  
  objective <- CVXR::Minimize(CVXR::sum_squares(y - x %*% w))
  constraints <- list(sum(w) == total, w >= lower, w <= upper)
  
  problem <- CVXR::Problem(objective, constraints)
  
  result <- CVXR::solve(problem)$getValue(w)
  
  return(result)
  
}
