test_that("equivalent to CVXR::solve", {
  
  skip("long-running test")
  
  packages <- c("roll", "zoo", "CVXR")
  
  status <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  
  if (!all(status)) {
    
    skip(paste0(paste0("'", paste0(names(status[!status]), collapse = ", '", sep = "'")),
                "' package(s) required for this test"))
    
  }
  
  # test data
  test_zoo_x <- c(lapply(test_ls, function(x){x[ , 1:2]}),
                  list("random vector with 0's" = test_ls[[1]][ , 1]))
  test_zoo_y <- list("random vector with 0's" = test_ls[[1]][ , 3])
  
  for (a in 1:length(test_ls)) {
    
    # test data
    test_mu <- roll::roll_mean(test_ls[[a]], test_width)
    test_sigma <- roll::roll_cov(test_ls[[a]], width = test_width)
    
    # rolling optimizations to minimize variance
    # cvxr_min_var(test_sigma[ , , n_obs])
    expect_equal(roll_min_var(test_sigma),
                 rollapplyr_optim(cvxr_min_var, sigma = test_sigma),
                 check.attributes = FALSE)
    
    # cvxr_min_var(test_sigma[ , , n_obs], test_mu[n_obs, ], test_target_mu[n_obs])
    expect_equal(roll_min_var(test_sigma, mu = test_mu,
                              target = test_target_mu),
                 rollapplyr_optim(cvxr_min_var, sigma = test_sigma,
                                  mu = test_mu, target = test_target_mu),
                 check.attributes = FALSE)
    
    # rolling optimizations to maximize mean
    # cvxr_max_mean(test_mu[n_obs, ])
    expect_equal(roll_max_mean(test_mu),
                 rollapplyr_optim(cvxr_max_mean, mu = test_mu),
                 check.attributes = FALSE, tolerance = 1e-5)
    
    # rolling optimizations to maximize utility
    # cvxr_max_utility(test_mu[n_obs, ], test_sigma[ , , n_obs])
    expect_equal(roll_max_utility(test_mu, test_sigma),
                 rollapplyr_optim(cvxr_max_utility, mu = test_mu,
                                  sigma = test_sigma),
                 check.attributes = FALSE)
    
  }
  
  for (ax in 1:length(test_zoo_x)) {
    for (ay in 1:length(test_zoo_y)) {
      
      # test data
      test_xx <- roll::roll_crossprod(test_zoo_x[[ax]], test_zoo_x[[ax]],
                                      test_width, min_obs = 1)
      test_xy <- roll::roll_crossprod(test_zoo_x[[ax]], test_zoo_y[[ay]],
                                      test_width, min_obs = 1)
      
      # rolling optimizations to minimize residual sum of squares
      # cvxr_min_rss(test_zoo_x[[ax]], test_zoo_y[[ay]])
      expect_equal(roll_min_rss(test_xx, test_xy),
                   rollapplyr_xy(cvxr_min_rss, test_zoo_x[[ax]],
                                 test_zoo_y[[ay]], test_width),
                   check.attributes = FALSE)
      
    }
  }
  
})