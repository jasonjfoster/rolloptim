set.seed(5640)
n_vars <- 3
n_obs <- 15
n_size <- n_obs * n_vars
dates <- rev(seq(Sys.Date(), length.out = n_obs, by = "-1 day"))

# test arguments
test_width <- c(1, 5, 10, 15)
test_target_mu <- rep(0, n_obs)
test_min_obs <- 1
# test_lambda <- 0.5
test_tol <- 1e-4

# test data
test_ls <- list("random matrix with 0's" =
                  matrix(sample(c(0, rnorm(n_size - 1)), n_size, replace = TRUE,
                                prob = rep(1 / n_size, n_size)),
                         nrow = n_obs, ncol = n_vars),
                
                "random matrix with 0's and NA's" =
                  matrix(sample(c(NA, 0, rnorm(n_size - 2)), n_size, replace = TRUE,
                                prob = c(1 / n_vars, rep(1 / n_size, n_size - 1))),
                         nrow = n_obs, ncol = n_vars))

if (requireNamespace("zoo", quietly = TRUE)) {
  test_ls[[1]] <- zoo::zoo(test_ls[[1]], dates)
}

test_ls[[1]] <- setNames(test_ls[[1]], paste0("x", rep(1:n_vars)))

# colnames(test_ls[[1]]) <- paste0("x", rep(1:n_vars))

if (requireNamespace("roll", quietly = TRUE)) {
  test_check_mu <- roll::roll_mean(test_ls[[1]], test_width[2])
  test_check_sigma <- roll::roll_cov(test_ls[[1]], width = test_width[2])
  test_check_xx <- roll::roll_crossprod(test_ls[[1]][ , 1:2], test_ls[[1]][ , 1:2], test_width[2])
  test_check_xy <- roll::roll_crossprod(test_ls[[1]][ , 1:2], test_ls[[1]][ , 3], test_width[2])
}
