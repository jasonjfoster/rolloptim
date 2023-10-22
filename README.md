# rolloptim

[![](https://github.com/jjf234/rollport/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/jjf234/rollport/actions/workflows/check-standard.yaml)
[![](https://codecov.io/gh/jjf234/rollport/graph/badge.svg)](https://app.codecov.io/github/jjf234/rollport)

## Overview

`rollport` is a package that provides analytical computation of rolling optimization for time-series data.

## Installation

Install the development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("jjf234/rollport")
```

## Usage

Load the package and supply a dataset:

``` r
library(rolloptim)

n_vars <- 3
n_obs <- 15
x <- matrix(rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars)

mu <- roll::roll_mean(x, 5)
sigma <- roll::roll_cov(x, width = 5)
```
Then, to compute rolling optimization, use the functions:

```r
# rolling optimization to minimize variance
roll_min_var(sigma)

# rolling optimization to maximize mean
roll_max_mean(mu)

# rolling optimization to maximize utility
roll_max_utility(mu, sigma, lambda = 1)
```

Note that handling of constraints is implemented by default (see the `total`, `lower`, and `upper` arguments).

## References

Markowitz, H.M. (1952). "Portfolio Selection." *The Journal of Finance*, 7(1), 77–91.

Tam, A. (2021). "Lagrangians and Portfolio Optimization." *https://www.adrian.idv.hk/2021-06-22-kkt/*.