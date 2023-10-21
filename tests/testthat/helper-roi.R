roi_min_var <- function(sigma, total = 1, lower = 0, upper = 1) {
  
  n_cols_sigma <- dim(sigma)[1]
  
  objective <- ROI::Q_objective(Q = 2 * sigma,
                                L = rep(0, n_cols_sigma))
  constraints <- ROI::L_constraint(L = rbind(rep(1, n_cols_sigma), rep(1, n_cols_sigma), rep(1, n_cols_sigma)),
                                   dir = c("==", ">=", "<="),
                                   rhs = c(total, lower, upper))
  
  problem <- ROI::OP(objective, constraints)
  
  result <- ROI::ROI_solve(problem, solver = "quadprog")
  
  return(ROI::solution(result))
  
}

roi_max_return <- function(mu, sigma, total = 1, lower = 0, upper = 1) {
  
  mu <- zoo::coredata(mu)
  
  n_cols_mu <- ncol(mu)
  
  objective <- ROI::L_objective(L = -mu)
  constraints <- ROI::L_constraint(L = rbind(rep(1, n_cols_mu), rep(1, n_cols_mu), rep(1, n_cols_mu)),
                                   dir = c("==", ">=", "<="),
                                   rhs = c(total, lower, upper))
  
  problem <- ROI::OP(objective, constraints)
  
  result <- ROI::ROI_solve(problem, solver = "glpk")
  
  return(ROI::solution(result))
  
}

roi_max_utility <- function(mu, sigma, gamma = 1, total = 1, lower = 0, upper = 1) {
  
  mu <- zoo::coredata(mu)
  
  n_cols_mu <- ncol(mu)
  
  objective <- ROI::Q_objective(Q = gamma * sigma,
                                L = -mu)
  constraints <- ROI::L_constraint(L = rbind(rep(1, n_cols_mu), rep(1, n_cols_mu), rep(1, n_cols_mu)),
                                   dir = c("==", ">=", "<="),
                                   rhs = c(total, lower, upper))
  
  problem <- ROI::OP(objective, constraints)
  
  result <- ROI::ROI_solve(problem, solver = "quadprog")
  
  return(ROI::solution(result))
  
}