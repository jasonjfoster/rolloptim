##' Rolling Optimizations to Minimize Variance
##'
##' A function for computing rolling optimizations to minimize variance.
##' 
##' @param sigma cube. Slices are covariance matrices.
##' @param total numeric. Sum of the weights.
##' @param lower numeric. Lower bound of the weights.
##' @param upper numeric. Upper bound of the weights.
##' @return An object of the same class and dimension as \code{mu} with the rolling
##' optimizations to minimize variance.
##' @examples
##' if (requireNamespace("roll", quietly = TRUE)) {
##' 
##' n_vars <- 3
##' n_obs <- 15
##' x <- matrix(rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars)
##'
##' sigma <- roll::roll_cov(x, width = 5)
##' 
##' # rolling portfolio optimizations to minimize variance
##' roll_min_var(sigma)
##' 
##' }
##' @export
roll_min_var <- function(sigma, total = 1, lower = 0, upper = 1) {
  return(.Call(`_rollport_roll_min_var`,
               sigma,
               as.numeric(total),
               as.numeric(lower),
               as.numeric(upper)
  ))
}

##' Rolling Optimizations to Maximize Mean
##'
##' A function for computing rolling optimizations to maximize mean.
##' 
##' @param mu matrix. Rows are means and columns are variables.
##' @param total numeric. Sum of the weights.
##' @param lower numeric. Lower bound of the weights.
##' @param upper numeric. Upper bound of the weights.
##' @return An object of the same class and dimension as \code{mu} with the rolling
##' optimizations to maximize mean.
##' @examples
##' if (requireNamespace("roll", quietly = TRUE)) {
##' 
##' n_vars <- 3
##' n_obs <- 15
##' x <- matrix(rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars)
##'
##' mu <- roll::roll_mean(x, 5)
##' 
##' # rolling optimizations to maximize mean
##' roll_max_mean(mu)
##' 
##' }
##' @export
roll_max_mean <- function(mu, total = 1, lower = 0, upper = 1) {
  return(.Call(`_rollport_roll_max_mean`,
               mu,
               as.numeric(total),
               as.numeric(lower),
               as.numeric(upper)
  ))
}

##' Rolling Portfolio Optimizations to Maximize Utility
##'
##' A function for computing rolling portfolio optimizations to maximize utility.
##' 
##' @param mu matrix. Rows are returns and columns are variables.
##' @param sigma cube. Slices are covariance matrices.
##' @param gamma numeric. Risk aversion parameter.
##' @param total numeric. Sum of the weights.
##' @param lower numeric. Lower bound of the weights.
##' @param upper numeric. Upper bound of the weights.
##' @return An object of the same class and dimension as \code{mu} with the rolling
##' portfolio optimizations to maximize utility.
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
##' # rolling portfolio optimizations to maximize utility
##' roll_max_utility(mu, sigma, gamma = 1)
##' 
##' }
##' @export
roll_max_utility <- function(mu, sigma, gamma = 1, total = 1, lower = 0, upper = 1) {
  return(.Call(`_rollport_roll_max_utility`,
               mu,
               sigma,
               as.numeric(gamma),
               as.numeric(total),
               as.numeric(lower),
               as.numeric(upper)
  ))
}