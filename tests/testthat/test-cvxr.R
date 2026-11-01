test_that("equivalent to CVXR::psolve", {
  
  packages <- c("roll", "zoo", "CVXR")
  
  status <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  
  if (!all(status)) {
    
    skip(paste0(paste0("'", paste0(names(status[!status]), collapse = ", '", sep = "'")),
                "' package(s) required for this test"))
    
  }
  
  # test data
  test_data_x <- c(lapply(test_ls, function(x){x[ , 1:2]}),
                   list("random vector with 0's and NA's" = test_ls[[2]][ , 1]))
  test_data_y <- c(lapply(test_ls, function(x){x[ , 3, drop = FALSE]}), # univariate 'y' for cvxr_min_rss
                   list("random vector with 0's and NA's" = test_ls[[2]][ , 3]))
  
  for (b in 1:length(test_width)) {

    width <- test_width[b]

    for (ax in 1:length(test_ls)) {
    
      # test data
      test_mu <- roll::roll_mean(test_ls[[ax]], width)
      test_sigma <- roll::roll_cov(test_ls[[ax]], width = width)
      
      # 'target' is feasible, on the boundary, or infeasible
      test_target <- list(test_target_mu,
                               apply(zoo::coredata(test_mu), 1, max),
                               apply(zoo::coredata(test_mu), 1, max) + 1)
      
      # rolling optimizations to minimize variance
      # cvxr_min_var(test_sigma[ , , n_obs])
      expect_equal(roll_min_var(test_sigma),
                  rollapplyr_optim(cvxr_min_var, sigma = test_sigma),
                  check.attributes = FALSE, tolerance = test_tol)
      
      for (g in 1:length(test_target)) {
        
        # cvxr_min_var(test_sigma[ , , n_obs], test_mu[n_obs, ], test_target[[g]][n_obs])
        expect_equal(roll_min_var(test_sigma, mu = test_mu,
                                  target = test_target[[g]]),
                     rollapplyr_optim(cvxr_min_var, sigma = test_sigma,
                                      mu = test_mu, target = test_target[[g]]),
                     check.attributes = FALSE, tolerance = test_tol)
        
      }
      
      # rolling optimizations to maximize mean
      # cvxr_max_mean(test_mu[n_obs, ])
      expect_equal(roll_max_mean(test_mu),
                   rollapplyr_optim(cvxr_max_mean, mu = test_mu),
                   check.attributes = FALSE, tolerance = test_tol)
      
      # rolling optimizations to maximize utility
      # cvxr_max_utility(test_mu[n_obs, ], test_sigma[ , , n_obs])
      expect_equal(roll_max_utility(test_mu, test_sigma),
                   rollapplyr_optim(cvxr_max_utility, mu = test_mu,
                                    sigma = test_sigma),
                   check.attributes = FALSE, tolerance = test_tol)
      
    }
    
    for (ax in 1:length(test_data_x)) {
      for (ay in 1:length(test_data_y)) {
        
        # test data
        test_xx <- roll::roll_crossprod(test_data_x[[ax]], test_data_x[[ax]],
                                        width, min_obs = test_min_obs[1])
        test_xy <- roll::roll_crossprod(test_data_x[[ax]], test_data_y[[ay]],
                                        width, min_obs = test_min_obs[1])
        
        # rolling optimizations to minimize residual sum of squares
        # cvxr_min_rss(test_data_x[[ax]], test_data_y[[ay]])
        expect_equal(roll_min_rss(test_xx, test_xy),
                    rollapplyr_xy(cvxr_min_rss, test_data_x[[ax]],
                                  test_data_y[[ay]], width, test_min_obs[1]),
                    check.attributes = FALSE, tolerance = test_tol)
        
      }
    }

  }
  
})
