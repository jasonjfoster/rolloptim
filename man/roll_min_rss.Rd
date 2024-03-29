% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rolloptim.R
\name{roll_min_rss}
\alias{roll_min_rss}
\title{Rolling Optimizations to Minimize Residual Sum of Squares}
\usage{
roll_min_rss(xx, xy, total = 1, lower = 0, upper = 1)
}
\arguments{
\item{xx}{cube. Slices are crossproducts of \code{x} and \code{x}.}

\item{xy}{cube. Slices are crossproducts of \code{x} and \code{y}.}

\item{total}{numeric. Sum of the weights.}

\item{lower}{numeric. Lower bound of the weights.}

\item{upper}{numeric. Upper bound of the weights.}
}
\value{
An object of the same class and dimension as \code{x} with the rolling
optimizations to minimize residual sum of squares.
}
\description{
A function for computing rolling optimizations to minimize residual sum of squares.
}
\examples{
if (requireNamespace("roll", quietly = TRUE)) {

n_vars <- 3
n_obs <- 15
x <- matrix(rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars)
y <- rnorm(n_obs)

xx <- roll::roll_crossprod(x, x, 5)
xy <- roll::roll_crossprod(x, y, 5)

# rolling optimizations to minimize residual sum of squares
roll_min_rss(xx, xy)

}
}
