# rolloptim

[![](https://github.com/jasonjfoster/rolloptim/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/jasonjfoster/rolloptim/actions/workflows/check-standard.yaml)
[![](https://codecov.io/gh/jasonjfoster/rolloptim/graph/badge.svg)](https://app.codecov.io/github/jasonjfoster/rolloptim)

## Overview

'rolloptim' provides analytical computation of rolling optimization for time-series data.

The 'rolloptim' package solves constrained quadratic and linear programs in closed form by applying Lagrangian multipliers and the Karush-Kuhn-Tucker conditions to perform mean-variance portfolio optimization (Markowitz, 1952, <doi:10.1111/j.1540-6261.1952.tb01525.x>) over rolling windows. For each window, the analytical solution computes the optimal weights that minimize variance, maximize expected return, minimize residual sum of squares, or maximize quadratic utility, subject to a total-weight equality constraint and box bounds on each weight. Use cases include:

* **Mean-variance optimization**: constructing portfolios that minimize variance or maximize utility subject to a budget and bounds
* **Expected-return maximization**: allocating to highest-mean assets under linear constraints
* **Constrained regression**: fitting a linear regression with sum-to-target and box constraints on the coefficients

The package supports rolling optimizations with constraints via the total, lower, and upper arguments. The implementation accepts rolling moments computed via the 'roll' package and uses 'RcppArmadillo' for linear algebra, with parallelism across windows provided by 'RcppParallel'.

## Installation

Install the released version from CRAN:

```r
# install.packages("rolloptim")
```

Or the development version from GitHub:

```r
# install.packages("pak")
pak::pak("jasonjfoster/rolloptim")
```

## Usage

Load the package and supply a dataset:

```r
library(rolloptim) # roll (>= 1.1.7)

n_vars <- 3
n_obs <- 15
x <- matrix(rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars)
y <- rnorm(n_obs)

mu <- roll::roll_mean(x, 5)
xx <- roll::roll_crossprod(x, x, 5)
xy <- roll::roll_crossprod(x, y, 5)
sigma <- roll::roll_cov(x, width = 5)
```

Then, to compute rolling optimizations, use the functions:

```r
# rolling optimization to minimize variance
roll_min_var(sigma)

# rolling optimization to maximize mean
roll_max_mean(mu)

# rolling optimization to minimize residual sum of squares
roll_min_rss(xx, xy)

# rolling optimization to maximize utility
roll_max_utility(mu, sigma, lambda = 1)
```

Constraint handling is implemented by default (see the `total`, `lower`, and `upper` arguments).

## References

Kuhn, H.W. and Tucker, A.W. (1951). "Nonlinear Programming." In *Proceedings of the Second Berkeley Symposium on Mathematical Statistics and Probability*, edited by J. Neyman, 481-492. University of California Press. <doi:10.1525/9780520411586-036>

Markowitz, H.M. (1952). "Portfolio Selection." *The Journal of Finance* 7 (1): 77-91. <doi:10.1111/j.1540-6261.1952.tb01525.x>
