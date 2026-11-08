#' Rolling Optimizations to Minimize Variance
#'
#' A function for computing rolling optimizations to minimize variance.
#' 
#' @param sigma cube or vector. Slices are covariance matrices, or rows are variances for one
#' variable if \code{mu} and \code{target} are supplied.
#' @param mu matrix. Rows are means and columns are variables.
#' @param target vector. Rows are target means.
#' @param total numeric. Sum of the weights.
#' @param lower numeric. Lower bound of the weights.
#' @param upper numeric. Upper bound of the weights.
#' @details Every combination of active lower and upper bounds is enumerated within each window,
#' the Karush-Kuhn-Tucker system is solved for each combination, and the feasible solution with the
#' smallest variance is kept. If \code{mu} and \code{target} are supplied then the expected return
#' must be greater than or equal to \code{target}.
#'
#' The number of combinations grows as \code{4^p} in the number of variables \code{p}, so \code{p}
#' is limited to 14. A window with no feasible solution, including one whose moments contain
#' missing values or whose \code{target} is unattainable, returns \code{NA}.
#' @return An object of the same class and dimension as \code{mu} with the rolling
#' optimizations to minimize variance, or a matrix with each row the rolling optimizations
#' if \code{mu} is \code{NULL}.
#' @examples
#' if (requireNamespace("roll", quietly = TRUE)) {
#' 
#' n_vars <- 3
#' n_obs <- 15
#' x <- matrix(rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars)
#'
#' sigma <- roll::roll_cov(x, width = 5)
#' 
#' # rolling optimizations to minimize variance
#' roll_min_var(sigma)
#' 
#' }
#' @export
roll_min_var <- function(sigma, mu = NULL, target = NULL,
                         total = 1, lower = 0, upper = 1) {
  return(.Call(`_rolloptim_roll_min_var`,
               sigma,
               mu,
               target,
               as.numeric(total),
               as.numeric(lower),
               as.numeric(upper)
  ))
}

#' Rolling Optimizations to Maximize Mean
#'
#' A function for computing rolling optimizations to maximize mean.
#' 
#' @param mu matrix. Rows are means and columns are variables.
#' @param sigma cube. Slices are covariance matrices.
#' @param target vector. Rows are target variances.
#' @param total numeric. Sum of the weights.
#' @param lower numeric. Lower bound of the weights.
#' @param upper numeric. Upper bound of the weights.
#' @details Every combination of active lower and upper bounds is enumerated within each window,
#' the Karush-Kuhn-Tucker system is solved for each combination, and the feasible solution with the
#' largest mean is kept.
#'
#' The number of combinations grows as \code{4^p} in the number of variables \code{p}, so \code{p}
#' is limited to 14. \code{sigma} and \code{target} are not supported, so a warning is given and the
#' arguments are ignored. A window with no feasible solution, including one whose moments contain
#' missing values, returns \code{NA}.
#' @return An object of the same class and dimension as \code{mu} with the rolling
#' optimizations to maximize mean.
#' @examples
#' if (requireNamespace("roll", quietly = TRUE)) {
#' 
#' n_vars <- 3
#' n_obs <- 15
#' x <- matrix(rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars)
#'
#' mu <- roll::roll_mean(x, 5)
#' 
#' # rolling optimizations to maximize mean
#' roll_max_mean(mu)
#' 
#' }
#' @export
roll_max_mean <- function(mu, sigma = NULL, target = NULL,
                          total = 1, lower = 0, upper = 1) {
  return(.Call(`_rolloptim_roll_max_mean`,
               mu,
               sigma,
               target,
               as.numeric(total),
               as.numeric(lower),
               as.numeric(upper)
  ))
}

#' Rolling Optimizations to Maximize Utility
#'
#' A function for computing rolling optimizations to maximize utility.
#' 
#' @param mu matrix. Rows are means and columns are variables.
#' @param sigma cube or vector. Slices are covariance matrices, or rows are variances for one
#' variable.
#' @param lambda numeric. Risk aversion parameter.
#' @param total numeric. Sum of the weights.
#' @param lower numeric. Lower bound of the weights.
#' @param upper numeric. Upper bound of the weights.
#' @details Every combination of active lower and upper bounds is enumerated within each window,
#' the Karush-Kuhn-Tucker system is solved for each combination, and the feasible solution with the
#' largest quadratic utility, the expected return less one half of \code{lambda} times the
#' variance, is kept.
#'
#' The number of combinations grows as \code{4^p} in the number of variables \code{p}, so \code{p}
#' is limited to 14. A window with no feasible solution, including one whose moments contain
#' missing values, returns \code{NA}.
#' @return An object of the same class and dimension as \code{mu} with the rolling
#' optimizations to maximize utility.
#' @examples
#' if (requireNamespace("roll", quietly = TRUE)) {
#' 
#' n_vars <- 3
#' n_obs <- 15
#' x <- matrix(rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars)
#'
#' mu <- roll::roll_mean(x, 5)
#' sigma <- roll::roll_cov(x, width = 5)
#' 
#' # rolling optimizations to maximize utility
#' roll_max_utility(mu, sigma, lambda = 1)
#' 
#' }
#' @export
roll_max_utility <- function(mu, sigma, lambda = 1, total = 1, lower = 0, upper = 1) {
  return(.Call(`_rolloptim_roll_max_utility`,
               mu,
               sigma,
               as.numeric(lambda),
               as.numeric(total),
               as.numeric(lower),
               as.numeric(upper)
  ))
}

#' Rolling Optimizations to Minimize Residual Sum of Squares
#'
#' A function for computing rolling optimizations to minimize residual sum of squares.
#' 
#' @param xx cube or vector. Slices are crossproducts of \code{x} and \code{x}, or rows are
#' crossproducts for one variable.
#' @param xy cube or vector. Slices are crossproducts of \code{x} and \code{y}, or rows are
#' crossproducts for one variable.
#' @param total numeric. Sum of the weights.
#' @param lower numeric. Lower bound of the weights.
#' @param upper numeric. Upper bound of the weights.
#' @details Every combination of active lower and upper bounds is enumerated within each window,
#' the Karush-Kuhn-Tucker system is solved for each combination, and the feasible solution with the
#' smallest residual sum of squares is kept, so the result is the constrained least squares
#' coefficients implied by the crossproducts \code{xx} and \code{xy}.
#'
#' The number of combinations grows as \code{4^p} in the number of variables \code{p}, so \code{p}
#' is limited to 14. A window with no feasible solution, including one whose crossproducts contain
#' missing values, returns \code{NA}.
#' @return A matrix with each row the rolling optimizations to minimize residual sum of squares.
#' @examples
#' if (requireNamespace("roll", quietly = TRUE)) {
#' 
#' n_vars <- 3
#' n_obs <- 15
#' x <- matrix(rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars)
#' y <- rnorm(n_obs)
#' 
#' xx <- roll::roll_crossprod(x, x, 5)
#' xy <- roll::roll_crossprod(x, y, 5)
#' 
#' # rolling optimizations to minimize residual sum of squares
#' roll_min_rss(xx, xy)
#' 
#' }
#' @export
roll_min_rss <- function(xx, xy, total = 1, lower = 0, upper = 1) {
  return(.Call(`_rolloptim_roll_min_rss`,
               xx,
               xy,
               as.numeric(total),
               as.numeric(lower),
               as.numeric(upper)
  ))
}
