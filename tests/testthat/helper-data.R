set.seed(5640)
n_vars <- 3
n_obs <- 15
n_size <- n_obs * n_vars
dates <- rev(seq(Sys.Date(), length.out = n_obs, by = "-1 day"))

# test arguments
test_width <- 5
test_target_mu <- rep(0, n_obs)
# test_lambda <- 0.5

# test data
test_ls <- list("random matrix with 0's" =
                  matrix(sample(c(0, rnorm(n_size - 1)), n_size, replace = TRUE,
                                prob = rep(1 / n_size, n_size)),
                         nrow = n_obs, ncol = n_vars))

if (requireNamespace("zoo", quietly = TRUE)) {
  test_ls[[1]] <- zoo::zoo(test_ls[[1]], dates)
}

test_ls[[1]] <- setNames(test_ls[[1]], paste0("x", rep(1:n_vars)))

# colnames(test_ls[[1]]) <- paste0("x", rep(1:n_vars))