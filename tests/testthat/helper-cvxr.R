cvxr_min_var <- function(sigma, total = 1, lower = 0, upper = 1) {
  
  n_cols_sigma <- dim(sigma)[1]
  
  w <- CVXR::Variable(n_cols_sigma)
  
  objective <- CVXR::Minimize(CVXR::quad_form(w, sigma))
  constraints <- list(sum(w) == total, w >= lower, w <= upper)
  
  problem <- CVXR::Problem(objective, constraints)
  
  result <- CVXR::solve(problem)$getValue(w)
  
  return(result)
  
}

cvxr_max_return <- function(mu, sigma, total = 1, lower = 0, upper = 1) {
  
  mu <- zoo::coredata(mu)
  
  n_cols_mu <- ncol(mu)
  
  w <- CVXR::Variable(n_cols_mu)
  
  objective <- CVXR::Maximize(mu %*% w)
  constraints <- list(sum(w) == total, w >= lower, w <= upper)
  
  problem <- CVXR::Problem(objective, constraints)
  
  result <- CVXR::solve(problem)$getValue(w)
  
  return(result)
  
}

cvxr_max_utility <- function(mu, sigma, gamma = 1, total = 1, lower = 0, upper = 1) {
  
  mu <- zoo::coredata(mu)
  
  n_cols_mu <- ncol(mu)
  
  w <- CVXR::Variable(n_cols_mu)
  
  objective <- CVXR::Minimize(0.5 * gamma * CVXR::quad_form(w, sigma) - mu %*% w)
  constraints <- list(sum(w) == total, w >= lower, w <= upper)
  
  problem <- CVXR::Problem(objective, constraints)
  
  result <- CVXR::solve(problem)$getValue(w)
  
  return(result)
  
}