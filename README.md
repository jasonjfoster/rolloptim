# rollport

## Overview

`rollport` is a package that provides fast and efficient computation of rolling portfolio optimization for time-series data.

## Installation

Install the development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("jjf234/rollport")
```

## Usage

Load the package and supply a dataset:

``` r
library(roll)
library(rollport)

n_vars <- 3
n_obs <- 15
x <- matrix(rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars)

mu <- roll_mean(x, 5)
sigma <- roll_cov(x, width = 5)
```
Then, to compute rolling portfolio optimization:

```r
# rolling portfolio optimization to minimize risks
roll_min_risk(mu, sigma)

# rolling portfolio optimization to maximize returns
roll_max_return(mu, sigma)

# rolling portfolio optimization to maximize ratios
roll_max_ratio(mu, sigma)
```
