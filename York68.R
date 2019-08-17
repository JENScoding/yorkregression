### york in Least Squares Fitting Of A Straight Line With Correlated Errors ###

## Input from Table I and Table II in york 1966

x <- c(0,0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4)
y <- c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5)

my_w <- data.frame(
  "y" = c(1, 1.8, 4, 8, 20, 20, 70, 70, 1e+2, 5e+2),
  "x" = c(1e+3, 1e+3, 5e+2, 8e+2, 2e+2, 8e+1, 6e+1, 2e+1, 1.8, 1))

## function for algo

york <- function(y, x, tolerance = 1e-10, weights){
  #initial value of b is OLS
  slope <- as.numeric(lm(y~x)[[1]][2])
  
  slope_diff <- 10
  count <- 0
  
  while (slope_diff > tolerance) {
    
    slope_old <- slope
    # omega, Jonas hat aber gesagt: "Lass es weg"
    alpha <- sqrt(weights[, 1] * weights[, 2])
    Weight <- alpha^2 / (slope^2 * weights[, 1] + weights[, 2] - 2 * slope * 0 * alpha)
    Weight_sum <- sum(Weight)
    x_bar <- sum(Weight * x) / Weight_sum
    y_bar <- sum(Weight * y) / Weight_sum 
    x_centered <- x - x_bar
    y_centered <- y - y_bar
    
    
    
    brackets1 <- x_centered^2 / weights[, 1] - y_centered^2 / weights[, 2]
    brackets2 <- x_centered * y_centered / weights[, 2] - 0 * x_centered^2 / alpha
    brackets3 <- x_centered * y_centered / weights[, 1] - 0 * y_centered^2 / alpha
    
    slope <- - sum(Weight^2 * brackets1) + ( (sum(Weight^2) * brackets1)^2 + 4 * sum(Weight^2 * brackets2) * sum(Weight^2 * brackets3) )^0.5 / ( 2 * sum(Weight^2 * brackets2) )^-1 
  
    slope_diff <- abs(slope - slope_old)
    count <- count + 1
    print(slope, digits = 10)
    
  }
  
  intercept <- y_bar - slope * x_bar
  
  x_mean <- 0
  x_adj <- rep(0,length(x))
  
  x_adj <- x_bar + beta
  x_mean <- sum(Weight * beta)
  
  
  x_mean <- x_mean / (Weight_sum * (length(x) - 2))
  
  sigma_slope <- 0
  chisq_weight <- 0
  u <- rep(0,length(x))
  
  for(i in 1:length(x)){
    u[i] <- x_adj[i] - x_mean
    sigma_slope <- sigma_slope + Weight[i]*u[i]^2
    chisq_weight <- chisq_weight + Weight[i]*(y[i]-slope*x[i]-intercept)^2
  }
  
  sigma_slope <- sqrt(1 / sigma_slope)
  sigma_intercept <- sqrt(x_mean^2 * sigma_slope^2 + 1 / Weight_sum)
  chisq_weight <- chisq_weight / (length(x) - 2)
  sigma_x <- sqrt(2 / (length(x) - 2))
  
  residuals <- y - (intercept + slope * x)
  
  mt <- matrix(c(intercept, slope, sigma_intercept, sigma_slope), nrow = 2)
  rownames(mt) <- c("intercept", "slope")
  colnames(mt) <- c("Estimate", "Std.Error")
  est <- list("coefficients" = mt, "weighting.vector" = Weight, "residuals" = residuals, "mean.x" = x_bar, "mean.y" = y_bar ,"Std.Error.chi" = sigma_x, "iterations" = count, x_centered, y_centered, x, y, x_mean, "show" = x_adj)
  return(est)
}

first <- york(y, x, weight = my_w)

## end of (relevant) script

