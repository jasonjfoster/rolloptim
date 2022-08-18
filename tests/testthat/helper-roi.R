# Error in `get(as.character(FUN), mode = "function", envir = envir)`:
# object 'set.portfolio.moments' of mode 'function' was not found
port_moments <- function(R) {
  
  result <- list()
  
  result$mu <- colMeans(R)
  result$sigma <- cov(R)
  
  return(result)
  
}

roi_min_risk <- function(x, total = 1, lower = 0, upper = 1) {
  
  port_spec <- PortfolioAnalytics::portfolio.spec(colnames(x))
  port_spec <- PortfolioAnalytics::add.constraint(portfolio = port_spec, type = "weight_sum", min_sum = total, max_sum = total)
  port_spec <- PortfolioAnalytics::add.constraint(portfolio = port_spec, type = "box", min = lower, max = upper)
  port_spec <- PortfolioAnalytics::add.objective(portfolio = port_spec, type = "risk", name = "var")
  
  result <- PortfolioAnalytics::optimize.portfolio(x, portfolio = port_spec, optimize_method = "ROI",
                                                   momentFUN = "port_moments")$weights
  
  return(result)
  
}

roi_max_return <- function(x, total = 1, lower = 0, upper = 1) {
  
  port_spec <- PortfolioAnalytics::portfolio.spec(colnames(x))
  port_spec <- PortfolioAnalytics::add.constraint(portfolio = port_spec, type = "weight_sum", min_sum = total, max_sum = total)
  port_spec <- PortfolioAnalytics::add.constraint(portfolio = port_spec, type = "box", min = lower, max = upper)
  port_spec <- PortfolioAnalytics::add.objective(portfolio = port_spec, type = "return", name = "mean")
  
  result <- PortfolioAnalytics::optimize.portfolio(x, portfolio = port_spec, optimize_method = "ROI",
                                                   momentFUN = "port_moments")$weights
  
  return(result)
  
}

roi_max_utility <- function(x, gamma = 1, total = 1, lower = 0, upper = 1) {
  
  port_spec <- PortfolioAnalytics::portfolio.spec(colnames(x))
  port_spec <- PortfolioAnalytics::add.constraint(portfolio = port_spec, type = "weight_sum", min_sum = total, max_sum = total)
  port_spec <- PortfolioAnalytics::add.constraint(portfolio = port_spec, type = "box", min = lower, max = upper)
  port_spec <- PortfolioAnalytics::add.objective(portfolio = port_spec, type = "return", name = "mean")
  port_spec <- PortfolioAnalytics::add.objective(portfolio = port_spec, type = "risk", name = "var", risk_aversion = gamma)
  
  result <- PortfolioAnalytics::optimize.portfolio(x, portfolio = port_spec, optimize_method = "ROI",
                                                   momentFUN = "port_moments")$weights
  
  return(result)
  
}