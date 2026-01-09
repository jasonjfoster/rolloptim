# install.packages("ROI.plugin.quadprog")
roi_min_var <- function(sigma, mu = NULL, target = NULL,
                        total = 1, lower = 0, upper = 1) {
  
  n_cols <- dim(sigma)[1]
  
  objective <- ROI::Q_objective(Q = 2 * sigma,
                                L = rep(0, n_cols))
  
  L <- rbind(rep(1, n_cols), rep(1, n_cols), rep(1, n_cols))
  dir <- c("==", ">=", "<=")
  rhs <- c(total, lower, upper)
  
  if (!is.null(mu) && !is.null(target)) {
    
    L <- rbind(L, as.numeric(mu))
    dir <- c(dir, ">=")
    rhs <- c(rhs, target)
  
  }

  constraints <- ROI::L_constraint(L = L, dir = dir, rhs = rhs)
  
  problem <- ROI::OP(objective, constraints)
  
  result <- ROI::ROI_solve(problem, solver = "quadprog")
  
  return(ROI::solution(result))
  
}

# install.packages("ROI.plugin.glpk")
roi_max_mean <- function(mu, total = 1, lower = 0, upper = 1) {
  
  # mu <- zoo::coredata(mu)
  mu <- as.numeric(zoo::coredata(mu))
  
  n_cols <- length(mu)
  
  objective <- ROI::L_objective(L = -mu)
  constraints <- ROI::L_constraint(L = rbind(rep(1, n_cols), rep(1, n_cols), rep(1, n_cols)),
                                   dir = c("==", ">=", "<="),
                                   rhs = c(total, lower, upper))
  
  problem <- ROI::OP(objective, constraints)
  
  result <- ROI::ROI_solve(problem, solver = "glpk")
  
  return(ROI::solution(result))
  
}

roi_max_utility <- function(mu, sigma, lambda = 1, total = 1, lower = 0, upper = 1) {
  
  mu <- zoo::coredata(mu)
  
  n_cols <- ncol(mu)
  
  objective <- ROI::Q_objective(Q = lambda * sigma,
                                L = -mu)
  constraints <- ROI::L_constraint(L = rbind(rep(1, n_cols), rep(1, n_cols), rep(1, n_cols)),
                                   dir = c("==", ">=", "<="),
                                   rhs = c(total, lower, upper))
  
  problem <- ROI::OP(objective, constraints)
  
  result <- ROI::ROI_solve(problem, solver = "quadprog")
  
  return(ROI::solution(result))
  
}

# install.packages("ROI.plugin.qpoases")
roi_min_rss <- function(x, y, total = 1, lower = 0, upper = 1) {
  
  x <- zoo::coredata(x)
  y <- zoo::coredata(y)
  
  n_cols <- ncol(x)
  
  objective <- ROI::Q_objective(Q = 2 * t(x) %*% x,
                                L = -2 * t(y) %*% x)
  constraints <- ROI::L_constraint(L = rbind(rep(1, n_cols), rep(1, n_cols), rep(1, n_cols)),
                                   dir = c("==", ">=", "<="),
                                   rhs = c(total, lower, upper))
  
  problem <- ROI::OP(objective, constraints)
  
  result <- ROI::ROI_solve(problem, solver = "qpoases")
  
  return(ROI::solution(result))
  
}