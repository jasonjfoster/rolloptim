packages <- c("PerformanceAnalytics", "roll")

if (all(vapply(packages, requireNamespace, logical(1), quietly = TRUE))) {
  
  # test arguments
  n_vars <- 3
  n_obs <- 15
  test_width <- 5
  test_gamma <- 0.5
  
  # test data
  data(edhec, package = "PerformanceAnalytics")
  test_x <- edhec[1:n_obs, 1:n_vars]
  test_mu <- roll::roll_mean(test_x, test_width)
  test_sigma <- roll::roll_cov(test_x, width = test_width)
  
}
