test_that("input checks", {
  
  packages <- c("roll")
  
  status <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  
  if (!all(status)) {
    
    skip(paste0("'", paste0(names(status[!status]), collapse = "', '"),
                "' package(s) required for this test"))
    
  }
  
  # check_cube
  expect_error(roll_min_var(test_check_sigma[ , , 1]),
               "'sigma' must be cube")
  expect_error(roll_min_var(test_check_sigma[ , , 1], mu = test_check_mu, target = test_target_mu),
               "'sigma' must be cube")
  expect_error(roll_max_utility(test_check_mu, test_check_sigma[ , , 1]),
               "'sigma' must be cube")
  
  # check_cube: vector 'sigma' (roll_min_var only)
  expect_error(roll_min_var(test_check_sigma[1, 1, ]),
               "'sigma' must be cube")
  
  # check_square
  expect_error(roll_min_var(test_check_sigma[ , 1:2, ]),
               "dimensions of 'sigma' must be square")
  expect_error(roll_min_var(test_check_sigma[1:2, , ], mu = test_check_mu, target = test_target_mu),
               "dimensions of 'sigma' must be square")
  expect_error(roll_max_utility(test_check_mu, test_check_sigma[1:2, , ]),
               "dimensions of 'sigma' must be square")
  
  # check_bounds
  expect_error(roll_min_var(test_check_sigma, lower = 1, upper = 0),
               "value of 'lower' must be less than or equal to value of 'upper'")
  expect_error(roll_max_mean(test_check_mu, lower = 1, upper = 0),
               "value of 'lower' must be less than or equal to value of 'upper'")
  expect_error(roll_max_utility(test_check_mu, test_check_sigma, lower = 1, upper = 0),
               "value of 'lower' must be less than or equal to value of 'upper'")
  expect_error(roll_min_rss(test_check_xx, test_check_xy, lower = 1, upper = 0),
               "value of 'lower' must be less than or equal to value of 'upper'")
  
  # check_finite: 'lower'
  expect_error(roll_min_var(test_check_sigma, lower = -Inf),
               "'lower' must be finite value")
  
  # check_finite: 'upper'
  expect_error(roll_min_var(test_check_sigma, upper = Inf),
               "'upper' must be finite value")
  
  # check_finite: 'total'
  expect_error(roll_min_var(test_check_sigma, total = Inf),
               "'total' must be finite value")
  
  # check_total
  expect_error(roll_min_var(test_check_sigma, total = n_vars + 1),
               "'total' must be between")
  expect_error(roll_max_mean(test_check_mu, total = n_vars + 1),
               "'total' must be between")
  expect_error(roll_max_utility(test_check_mu, test_check_sigma, total = n_vars + 1),
               "'total' must be between")
  expect_error(roll_min_rss(test_check_xx, test_check_xy, total = n_vars + 1),
               "'total' must be between")
  
  # check_target
  expect_error(roll_min_var(test_check_sigma, mu = test_check_mu),
               "requires both 'mu' and 'target'")
  
  # check_rows_equal
  expect_error(roll_min_var(test_check_sigma, mu = test_check_mu[-1, ], target = test_target_mu[-1]),
               "number of rows in 'mu' must equal the number of rows in 'sigma'")
  expect_error(roll_max_utility(test_check_mu[-1, ], test_check_sigma),
               "number of rows in 'mu' must equal the number of rows in 'sigma'")
  
  # check_univariate
  expect_error(roll_min_var(test_check_sigma[1, 1, ], mu = test_check_mu, target = test_target_mu),
               "number of columns in 'mu' must be one for univariate 'sigma'")
  expect_error(roll_max_utility(test_check_mu, test_check_sigma[1, 1, ]),
               "number of columns in 'mu' must be one for univariate 'sigma'")
  
  # check_rows_equal: vector 'sigma'
  expect_error(roll_min_var(test_check_sigma[1, 1, ], mu = test_check_mu[-1, 1, drop = FALSE],
                            target = test_target_mu[-1]),
               "number of rows in 'mu' must equal the number of rows in 'sigma'")
  expect_error(roll_max_utility(test_check_mu[-1, 1, drop = FALSE], test_check_sigma[1, 1, ]),
               "number of rows in 'mu' must equal the number of rows in 'sigma'")
  
  # check_cols_equal
  expect_error(roll_min_var(test_check_sigma, mu = test_check_mu[ , 1:2], target = test_target_mu),
               "number of columns in 'mu' must equal the number of columns in 'sigma'")
  expect_error(roll_max_utility(test_check_mu[ , 1:2], test_check_sigma),
               "number of columns in 'mu' must equal the number of columns in 'sigma'")
  
  # check_nonneg_finite: 'lambda'
  expect_error(roll_max_utility(test_check_mu, test_check_sigma, lambda = -1),
               "'lambda' must be non-negative finite value")
  
  # check_cube: 'xx'
  expect_error(roll_min_rss(test_check_xx[ , , 1], test_check_xy),
               "'xx' must be cube")
  
  # check_cube: 'xy'
  expect_error(roll_min_rss(test_check_xx, test_check_xy[ , 1, ]),
               "'xy' must be cube")
  
  # check_cube: vector 'xy' with matrix 'xx'
  expect_error(roll_min_rss(test_check_xx[ , , 1], test_check_xy[1, 1, ]),
               "'xx' must be cube")
  
  # check_cube: vector 'xx' with matrix 'xy'
  expect_error(roll_min_rss(test_check_xx[1, 1, ], test_check_xy[ , 1, ]),
               "'xy' must be cube")
  
  # check_univariate: vector 'xy' with cube 'xx'
  expect_error(roll_min_rss(test_check_xx, test_check_xy[1, 1, ]),
               "number of columns in 'xx' must be one for univariate 'xy'")
  
  # check_univariate: vector 'xx' with cube 'xy'
  expect_error(roll_min_rss(test_check_xx[1, 1, ], test_check_xy),
               "number of columns in 'xy' must be one for univariate 'xx'")
  
  # check_square: vector 'xy' with cube 'xx'
  expect_error(roll_min_rss(test_check_xx[ , 1, , drop = FALSE], test_check_xy[1, 1, ]),
               "dimensions of 'xx' must be square")
  
  # check_rows_equal: vector 'xx' with cube 'xy'
  expect_error(roll_min_rss(test_check_xx[1, 1, ], test_check_xy[1, 1, -1, drop = FALSE]),
               "number of rows in 'xy' must equal the number of rows in 'xx'")
  
  # check_rows_equal: vector 'xx' and 'xy'
  expect_error(roll_min_rss(test_check_xx[1, 1, ], test_check_xy[1, 1, -1]),
               "number of rows in 'xy' must equal the number of rows in 'xx'")
  
  # check_square: 'xx'
  expect_error(roll_min_rss(test_check_xx[ , 1, , drop = FALSE], test_check_xy),
               "dimensions of 'xx' must be square")
  
  # check_rows_equal: 'xy'
  expect_error(roll_min_rss(test_check_xx, test_check_xy[ , , -1, drop = FALSE]),
               "number of rows in 'xy' must equal the number of rows in 'xx'")
  
  # check_cols_equal: 'xy'
  expect_error(roll_min_rss(test_check_xx, test_check_xy[1, , , drop = FALSE]),
               "number of columns in 'xy' must equal the number of columns in 'xx'")
  
  # check_bounds_int: 'n_cols'
  expect_error(roll_min_var(array(0, c(15, 15, n_obs))),
               "value of 'n_cols' must be between 1 and 14")
  expect_error(roll_min_var(array(0, c(15, 15, n_obs)), mu = matrix(0, n_obs, 15), target = test_target_mu),
               "value of 'n_cols' must be between 1 and 14")
  expect_error(roll_max_mean(matrix(0, n_obs, 15)),
               "value of 'n_cols' must be between 1 and 14")
  expect_error(roll_max_utility(matrix(0, n_obs, 15), array(0, c(15, 15, n_obs))),
               "value of 'n_cols' must be between 1 and 14")
  expect_error(roll_min_rss(array(0, c(15, 15, n_obs)), array(0, c(15, 1, n_obs))),
               "value of 'n_cols' must be between 1 and 14")
  
})

test_that("input warnings", {
  
  packages <- c("roll")
  
  status <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  
  if (!all(status)) {
    
    skip(paste0("'", paste0(names(status[!status]), collapse = "', '"),
                "' package(s) required for this test"))
    
  }
  
  # 'sigma' (roll_max_mean only)
  expect_warning(roll_max_mean(test_check_mu, sigma = test_check_sigma),
                 "'sigma' is only supported for linear constraints")
  
  # 'target' (roll_max_mean only)
  expect_warning(roll_max_mean(test_check_mu, target = test_target_mu),
                 "'target' is only supported for linear constraints")
  
})
