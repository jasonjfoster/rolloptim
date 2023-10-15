test_that("equivalent to CVXR::solve", {
  
  skip("long-running test")
  
  packages <- c("roll", "zoo", "CVXR")
  
  status <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  
  if (!all(status)) {
    
    skip(paste0(paste0("'", paste0(names(status[!status]), collapse = ", '", sep = "'")),
                "' package(s) required for this test"))
    
  }
  
  for (a in 1:(length(test_ls))) {
    
    # test data
    test_mu <- roll::roll_mean(test_ls[[a]], test_width)
    test_sigma <- roll::roll_cov(test_ls[[a]], width = test_width)
    
    # rolling portfolio optimizations to minimize risk
    expect_equal(roll_min_risk(test_mu, test_sigma),
                 rollapplyr_port(cvxr_min_risk, test_mu, test_sigma),
                 check.attributes = FALSE)
    
    # rolling portfolio optimizations to maximize return
    expect_equal(roll_max_return(test_mu, test_sigma),
                 rollapplyr_port(cvxr_max_return, test_mu, test_sigma),
                 check.attributes = FALSE)
    
    # rolling portfolio optimizations to maximize utility
    expect_equal(roll_max_utility(test_mu, test_sigma, test_gamma),
                 rollapplyr_port(cvxr_max_utility, test_mu, test_sigma, test_gamma),
                 check.attributes = FALSE)
    
  }
  
})