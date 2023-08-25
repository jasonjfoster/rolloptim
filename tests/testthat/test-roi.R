# https://cran.r-project.org/web/packages/PortfolioAnalytics/vignettes/ROI_vignette.pdf
test_that("equivalent to PortfolioAnalytics::optimize.portfolio", {
  
  packages <- c("roll", "zoo", "PortfolioAnalytics", "PerformanceAnalytics", "ROI", "ROI.plugin.quadprog", "ROI.plugin.glpk")
  
  status <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  
  if (!all(status)) {
    
    skip(paste0(paste0("'", paste0(names(status[!status]), collapse = ", '", sep = "'")),
         "' package(s) required for this test"))
    
  }
  
  # test arguments
  n_vars <- 3
  n_obs <- 15
  test_width <- 5
  test_gamma <- 0.5
  
  # test data
  data(edhec, package = "PerformanceAnalytics") # "object 'edhec' not found" in helper-data.R
  test_x <- edhec[1:n_obs, 1:n_vars]
  test_mu <- roll::roll_mean(test_x, test_width)
  test_sigma <- roll::roll_cov(test_x, width = test_width)
  
  # rolling portfolio optimizations to minimize risk
  expect_equal(roll_min_risk(test_mu, test_sigma),
               zoo::rollapplyr(test_x, width = test_width, roi_min_risk, by.column = FALSE),
               check.attributes = FALSE)
  
  # rolling portfolio optimizations to maximize return
  expect_equal(roll_max_return(test_mu, test_sigma),
               zoo::rollapplyr(test_x, width = test_width, roi_max_return, by.column = FALSE),
               check.attributes = FALSE)
  
  # rolling portfolio optimizations to maximize utility
  expect_equal(roll_max_utility(test_mu, test_sigma, test_gamma),
               zoo::rollapplyr(test_x, width = test_width, roi_max_utility, gamma = test_gamma / 2, by.column = FALSE),
               check.attributes = FALSE)
  
})