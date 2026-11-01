test_that("input checks", {
  
  packages <- c("roll")
  
  status <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  
  if (!all(status)) {
    
    skip(paste0(paste0("'", paste0(names(status[!status]), collapse = ", '", sep = "'")),
                "' package(s) required for this test"))
    
  }
  
  # test data
  width <- test_width[2]
  test_check_mu <- roll::roll_mean(test_ls[[1]], width)
  test_check_sigma <- roll::roll_cov(test_ls[[1]], width = width)
  test_check_xx <- roll::roll_crossprod(test_ls[[1]][ , 1:2], test_ls[[1]][ , 1:2], width)
  test_check_xy <- roll::roll_crossprod(test_ls[[1]][ , 1:2], test_ls[[1]][ , 3], width)
  
  # check_square
  expect_error(roll_min_var(test_check_sigma[ , 1:2, ]),
               "dimensions of 'sigma' must be square")
  
  # check_bounds
  expect_error(roll_min_var(test_check_sigma, lower = 1, upper = 0),
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
  
  # check_target
  expect_error(roll_min_var(test_check_sigma, mu = test_check_mu),
               "requires both 'mu' and 'target'")
  
  # check_rows_equal
  expect_error(roll_min_var(test_check_sigma, mu = test_check_mu[-1, ], target = test_target_mu[-1]),
               "number of rows in 'mu' must equal the number of rows in 'sigma'")
  expect_error(roll_max_utility(test_check_mu[-1, ], test_check_sigma),
               "number of rows in 'mu' must equal the number of rows in 'sigma'")
  
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
  
  # check_square: 'xx' (named 'sigma' in roll_min_rss)
  expect_error(roll_min_rss(test_check_xx[ , 1, , drop = FALSE], test_check_xy),
               "dimensions of 'sigma' must be square")
  
  # check_rows_equal: 'xy' (named 'mu' in roll_min_rss)
  expect_error(roll_min_rss(test_check_xx, test_check_xy[ , , -1, drop = FALSE]),
               "number of rows in 'mu' must equal the number of rows in 'sigma'")
  
  # check_cols_equal: 'xy' (named 'mu' in roll_min_rss)
  expect_error(roll_min_rss(test_check_xx, test_check_xy[1, , , drop = FALSE]),
               "number of columns in 'mu' must equal the number of columns in 'sigma'")
  
})

test_that("input warnings", {
  
  packages <- c("roll")
  
  status <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  
  if (!all(status)) {
    
    skip(paste0(paste0("'", paste0(names(status[!status]), collapse = ", '", sep = "'")),
                "' package(s) required for this test"))
    
  }
  
  # test data
  width <- test_width[2]
  test_check_mu <- roll::roll_mean(test_ls[[1]], width)
  test_check_sigma <- roll::roll_cov(test_ls[[1]], width = width)
  
  # 'sigma' (roll_max_mean only)
  expect_warning(roll_max_mean(test_check_mu, sigma = test_check_sigma),
                 "'sigma' is only supported for linear constraints")
  
  # 'target' (roll_max_mean only)
  expect_warning(roll_max_mean(test_check_mu, target = test_target_mu),
                 "'target' is only supported for linear constraints")
  
})
