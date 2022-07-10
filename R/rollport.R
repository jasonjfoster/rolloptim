##' Rolling Portfolio Optimizations to Minimize Risk
##'
##' A function for computing the rolling portfolio optimizations to minimize risk.
##' 
##' @param mu matrix. Rows are returns and columns are variables.
##' @param sigma cube. Slices are covariance matrices.
##' @param total numeric. Sum of the weights.
##' @param lower numeric. Lower boundry of the weights.
##' @param upper numeric. Upper boundry of the weights.
##' @examples
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

##' Rolling Portfolio Optimizations to Maximize Return
##'
##' A function for computing the rolling portfolio optimizations to maximize return.
##' 
##' @param mu matrix. Rows are returns and columns are variables.
##' @param sigma cube. Slices are covariance matrices.
##' @param total numeric. Sum of the weights.
##' @param lower numeric. Lower boundry of the weights.
##' @param upper numeric. Upper boundry of the weights.
##' @examples
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

##' Rolling Portfolio Optimizations to Maximize Ratio
##'
##' A function for computing the rolling portfolio optimizations to maximize ratio
##' 
##' @param mu matrix. Rows are returns and columns are variables.
##' @param sigma cube. Slices are covariance matrices.
##' @param total numeric. Sum of the weights.
##' @param lower numeric. Lower boundry of the weights.
##' @param upper numeric. Upper boundry of the weights.
##' @examples
##' @export
roll_max_ratio <- function(mu, sigma, total = 1, lower = 0, upper = 1) {
  return(.Call(`_rollport_roll_max_ratio`,
               mu,
               sigma,
               as.numeric(total),
               as.numeric(lower),
               as.numeric(upper)
  ))
}