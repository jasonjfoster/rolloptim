test_that("solutions on the bounds", {
  
  packages <- c("roll", "zoo")
  
  status <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  
  if (!all(status)) {
    
    skip(paste0("'", paste0(names(status[!status]), collapse = "', '"),
                "' package(s) required for this test"))
    
  }
  
  for (b in 1:length(test_width)) {
    
    width <- test_width[b]
    
    for (ax in 1:length(test_ls)) {
      
      # test data
      test_mu <- roll::roll_mean(test_ls[[ax]], width)
      test_sigma <- roll::roll_cov(test_ls[[ax]], width = width)
      
      # 'target' on the boundary: the only feasible solution is the maximum of 'mu'
      test_target <- apply(zoo::coredata(test_mu), 1, max)
      test_weights <- t(apply(zoo::coredata(test_mu), 1, function(x) as.numeric(x == max(x))))
      test_weights[is.na(test_sigma[1, 1, ]), ] <- NA
      attributes(test_weights) <- attributes(test_mu)
      
      # rolling optimizations to minimize variance
      expect_equal(roll_min_var(test_sigma, mu = test_mu, target = test_target),
                   test_weights, check.attributes = FALSE)
      
    }
    
  }
  
})

test_that("missing values", {
  
  packages <- c("roll")
  
  status <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  
  if (!all(status)) {
    
    skip(paste0("'", paste0(names(status[!status]), collapse = "', '"),
                "' package(s) required for this test"))
    
  }
  
  # test data
  width <- test_width[2]
  test_sigma <- roll::roll_cov(test_ls[[2]], width = width)
  
  # a slice with missing values returns a row of missing values
  expect_true(all(is.na(roll_min_var(test_sigma)[is.na(test_sigma[1, 1, ]), ])))
  
})

test_that("univariate", {
  
  packages <- c("roll")
  
  status <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  
  if (!all(status)) {
    
    skip(paste0("'", paste0(names(status[!status]), collapse = "', '"),
                "' package(s) required for this test"))
    
  }
  
  # test data
  width <- test_width[2]
  test_mu <- roll::roll_mean(test_ls[[2]][ , 1, drop = FALSE], width)
  test_sigma <- roll::roll_cov(test_ls[[2]][ , 1], width = width)
  test_xx <- roll::roll_crossprod(test_ls[[2]][ , 1], test_ls[[2]][ , 1], width)
  test_xy <- roll::roll_crossprod(test_ls[[2]][ , 1], test_ls[[2]][ , 3], width)
  test_xy_cube <- roll::roll_crossprod(test_ls[[2]][ , 1], test_ls[[2]][ , 3, drop = FALSE], width)
  
  # univariate 'sigma': the only feasible solution is 'total'
  test_target <- test_mu[ , 1]
  test_weights <- matrix(1, n_obs, 1)
  test_weights[is.na(test_sigma), ] <- NA
  
  # rolling optimizations to minimize variance
  expect_equal(roll_min_var(test_sigma, mu = test_mu, target = test_target),
               test_weights, check.attributes = FALSE)
  
  # rolling optimizations to maximize utility
  expect_equal(roll_max_utility(test_mu, test_sigma),
               test_weights, check.attributes = FALSE)
  
  # univariate 'xx' and 'xy': the only feasible solution is 'total' (or missing values in 'xy')
  test_weights[is.na(test_xy), ] <- NA
  
  # rolling optimizations to minimize residual sum of squares
  expect_equal(roll_min_rss(test_xx, test_xy),
               test_weights, check.attributes = FALSE)
  expect_equal(roll_min_rss(test_xx, test_xy_cube),
               test_weights, check.attributes = FALSE)
  
})
