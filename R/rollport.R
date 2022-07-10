##' Rolling Portfolio Optimization to Minimize Risks
##'
##' A function for computing rolling portfolio optimization to minimize risks.
##' 
##' @param mu matrix. Rows are returns and columns are variables.
##' @param sigma cube. Slices are covariance matrices.
##' @param total numeric. Sum of the weights.
##' @param lower numeric. Lower bound of the weights.
##' @param upper numeric. Upper bound of the weights.
##' @examples
##' if (requireNamespace("roll", quietly = TRUE)) {
##' 
##' n_vars <- 3
##' n_obs <- 15
##' x <- matrix(rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars)
##'
##' mu <- roll::roll_mean(x, 5)
##' sigma <- roll::roll_cov(x, width = 5)
##' 
##' # rolling portfolio optimization to minimize risks
##' roll_min_risk(mu, sigma)
##' 
##' }
##' @export
roll_min_risk <- function(mu, sigma, total = 1, lower = 0, upper = 1) {
  return(.Call(`_rollport_roll_min_risk`,
               mu,
               sigma,
               as.numeric(total),
               as.numeric(lower),
               as.numeric(upper)
  ))
}

##' Rolling Portfolio Optimization to Maximize Returns
##'
##' A function for computing rolling portfolio optimization to maximize returns.
##' 
##' @param mu matrix. Rows are returns and columns are variables.
##' @param sigma cube. Slices are covariance matrices.
##' @param total numeric. Sum of the weights.
##' @param lower numeric. Lower bound of the weights.
##' @param upper numeric. Upper bound of the weights.
##' @examples
##' if (requireNamespace("roll", quietly = TRUE)) {
##' 
##' n_vars <- 3
##' n_obs <- 15
##' x <- matrix(rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars)
##'
##' mu <- roll::roll_mean(x, 5)
##' sigma <- roll::roll_cov(x, width = 5)
##' 
##' # rolling portfolio optimization to maximize returns
##' roll_max_return(mu, sigma)
##' 
##' }
##' @export
roll_max_return <- function(mu, sigma, total = 1, lower = 0, upper = 1) {
  return(.Call(`_rollport_roll_max_return`,
               mu,
               sigma,
               as.numeric(total),
               as.numeric(lower),
               as.numeric(upper)
  ))
}

##' Rolling Portfolio Optimization to Maximize Ratios
##'
##' A function for computing rolling portfolio optimization to maximize ratios.
##' 
##' @param mu matrix. Rows are returns and columns are variables.
##' @param sigma cube. Slices are covariance matrices.
##' @param gamma numeric. Risk aversion parameter.
##' @param total numeric. Sum of the weights.
##' @param lower numeric. Lower bound of the weights.
##' @param upper numeric. Upper bound of the weights.
##' @examples
##' if (requireNamespace("roll", quietly = TRUE)) {
##' 
##' n_vars <- 3
##' n_obs <- 15
##' x <- matrix(rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars)
##'
##' mu <- roll::roll_mean(x, 5)
##' sigma <- roll::roll_cov(x, width = 5)
##' 
##' # rolling portfolio optimization to maximize ratios
##' roll_max_ratio(mu, sigma)
##' 
##' }
##' @export
roll_max_ratio <- function(mu, sigma, gamma = 1, total = 1, lower = 0, upper = 1) {
  return(.Call(`_rollport_roll_max_ratio`,
               mu,
               sigma,
               as.numeric(gamma),
               as.numeric(total),
               as.numeric(lower),
               as.numeric(upper)
  ))
}