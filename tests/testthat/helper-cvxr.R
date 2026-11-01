cvxr_min_var <- function(sigma, mu = NULL, target = NULL,
                         total = 1, lower = 0, upper = 1) {
  
  n_cols <- dim(sigma)[1]
  
  w <- CVXR::Variable(n_cols)
  
  objective <- CVXR::Minimize(CVXR::quad_form(w, sigma))
  constraints <- list(sum(w) == total, w >= lower, w <= upper)
  
  if (!is.null(mu) && !is.null(target)) {
    mu <- as.numeric(zoo::coredata(mu))
    constraints <- c(constraints, list(t(mu) %*% w >= target))
  }

  problem <- CVXR::Problem(objective, constraints)
  
  CVXR::psolve(problem)
  
  # CVXR: solutions are returned even if the status is not optimal
  if (CVXR::status(problem) != "optimal") {
    stop("status is not optimal")
  }
  
  result <- CVXR::value(w)
  
  return(result)
  
}

cvxr_max_mean <- function(mu, sigma = NULL, target = NULL,
                          total = 1, lower = 0, upper = 1) {
  
  mu <- as.numeric(zoo::coredata(mu))
  
  n_cols <- length(mu)
  
  w <- CVXR::Variable(n_cols)
  
  objective <- CVXR::Maximize(t(mu) %*% w)
  constraints <- list(sum(w) == total, w >= lower, w <= upper)

  if (!is.null(sigma) && !is.null(target)) {
    constraints <- c(constraints, list(CVXR::quad_form(w, sigma) <= target))
  }

  problem <- CVXR::Problem(objective, constraints)
  
  CVXR::psolve(problem)
  
  # CVXR: solutions are returned even if the status is not optimal
  if (CVXR::status(problem) != "optimal") {
    stop("status is not optimal")
  }
  
  result <- CVXR::value(w)
  
  return(result)
  
}

cvxr_max_utility <- function(mu, sigma, lambda = 1, total = 1, lower = 0, upper = 1) {
  
  mu <- as.numeric(zoo::coredata(mu))
  
  n_cols <- length(mu)
  
  w <- CVXR::Variable(n_cols)
  
  objective <- CVXR::Minimize(0.5 * lambda * CVXR::quad_form(w, sigma) - t(mu) %*% w)
  constraints <- list(sum(w) == total, w >= lower, w <= upper)
  
  problem <- CVXR::Problem(objective, constraints)
  
  CVXR::psolve(problem)
  
  # CVXR: solutions are returned even if the status is not optimal
  if (CVXR::status(problem) != "optimal") {
    stop("status is not optimal")
  }
  
  result <- CVXR::value(w)
  
  return(result)
  
}

cvxr_min_rss <- function(x, y, total = 1, lower = 0, upper = 1) {
  
  x <- as.matrix(zoo::coredata(x))
  y <- as.matrix(zoo::coredata(y))
  
  n_cols <- ncol(x)
  
  # roll::roll_crossprod: complete observations of 'x' for 'xx' and of 'x' and 'y' for 'xy'
  ix_xx <- complete.cases(x)
  ix_xy <- complete.cases(x, y)
  xx <- crossprod(x[ix_xx, , drop = FALSE])
  xy <- crossprod(x[ix_xy, , drop = FALSE], y[ix_xy, , drop = FALSE])
  
  w <- CVXR::Variable(n_cols)
  
  objective <- CVXR::Minimize(CVXR::quad_form(w, xx) - 2 * t(xy) %*% w)
  constraints <- list(sum(w) == total, w >= lower, w <= upper)
  
  problem <- CVXR::Problem(objective, constraints)
  
  CVXR::psolve(problem)
  
  # CVXR: solutions are returned even if the status is not optimal
  if (CVXR::status(problem) != "optimal") {
    stop("status is not optimal")
  }
  
  result <- CVXR::value(w)
  
  return(result)
  
}
