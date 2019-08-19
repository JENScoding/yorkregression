### york in Least Squares Fitting Of A Straight Line With Correlated Errors ###

## Input from Table I and Table II in york 1966
rm(list=ls())
x1 <- c(0, 0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4)
y1 <- c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5)

x2 <- c(0, 0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4)
y2 <- c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5)

x3 <- c(0, 0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4)
y3 <- c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5)

x4 <- c(0, 0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4)
y4 <- c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5)

x5 <- c(0, 0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4)
y5 <- c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5)

weights.y = c(1, 1.8, 4, 8, 20, 20, 70, 70, 1e+2, 5e+2)
weights.x = c(1e+3, 1e+3, 5e+2, 8e+2, 2e+2, 8e+1, 6e+1, 2e+1, 1.8, 1)

## function for algo

york <- function(x, y, tolerance = 1e-10, weights.x, weights.y){
  #initial value of b is OLS
  lm.OLS <- lm(y~x)
  slope <- as.numeric(lm.OLS[[1]][2])
  
  slope.diff <- 10
  count <- 0
  slope.per.iteration <- NULL
  
  while (slope.diff > tolerance) {
    
    slope.old <- slope
    # omega, Jonas hat aber gesagt: "Lass es weg"
    alpha <- sqrt(weights.x * weights.y)
    Weight <- alpha^2 / (slope^2 * weights.y + weights.x - 2 * slope * 0 * alpha)
    Weight.sum <- sum(Weight)
    x.bar <- sum(Weight * x, na.rm = T) / Weight.sum
    y.bar <- sum(Weight * y, na.rm = T) / Weight.sum
    x.centered <- x - x.bar
    y.centered <- y - y.bar
    
    beta <- Weight * ((x.centered / weights.y) + (slope * y.centered / weights.x) - (slope * x.centered + y.centered) * 0 / alpha)
    Q1 <- sum(Weight * beta * y.centered, na.rm = T)
    Q2 <- sum(Weight * beta * x.centered, na.rm = T)
    slope <- Q1 / Q2
    slope.diff <- abs(slope - slope.old)
    count <- count + 1
    slope.per.iteration <- append(slope.per.iteration, slope)
    
    if (count > tolerance^-1)  stop("\nThe slope coefficient does not converge after ",
                                    count,
                                    " iterations. \nHint: You may reduce the tolerance level.",
                                    cat("Slope coefficient for last 5 iterations:"),
                                    for (i in 4:0){
                                      cat("\n\t", count - i, "\t", slope.per.iteration[count - i])
                                    },
                                    cat("\n")
    )
  }
  slope.per.iteration <- data.frame("slope.per.iteration" = slope.per.iteration)
  
  intercept <- y.bar - slope * x.bar
  
  x.adj <- x.bar + beta
  
  x.mean <- sum(Weight * beta, na.rm = T) / (Weight.sum * (length(x) - 2))
  u <- x.adj - x.mean
  
  sigma.slope <- sqrt(1 / sum(Weight * u^2, na.rm = T))
  sigma.intercept <- sqrt(x.mean^2 * sigma.slope^2 + 1 / Weight.sum)
  
  chisq.weight <- sum(Weight * (y - slope * x - intercept)^2, na.rm = T) / (length(x) - 2)
  sigma.x <- sqrt(2 / (length(x) - 2))
  
  fitted.y <- intercept + slope * x
  
  residuals <- y - fitted.y
  
  c <- 0*alpha
  x.residuals <- (Weight * (intercept + slope * x - y) * (c - slope * weights.y)) / (weights.y * weights.x)
  y.residuals <- (Weight * (intercept + slope * x - y) * (weights.x - slope * c))/ (weights.y * weights.x)
  
  mt <- matrix(c(intercept, slope, sigma.intercept, sigma.slope), nrow = 2)
  rownames(mt) <- c("intercept", "slope")
  colnames(mt) <- c("Estimate", "Std.Error")
  est <- list("coefficients" = mt,
              "weighting.vector" = Weight,
              "x.residuals" = x.residuals,
              "y.residuals"=y.residuals,
              "fitted.y"=fitted.y,
              "mean.x" = x.bar,
              "mean.y" = y.bar ,
              "Std.Error.chi" = sigma.x,
              "Weighted_chi^2" = chisq.weight,
              "number.of.iterations" = count,
              "slope.after.each.iteration" = slope.per.iteration,
              x.centered, y.centered, x, y, x.mean, "show" = x.adj,
              "original.x.values" = x,
              "original.y.values" = y,
              "fitted.ols" = lm.OLS$fitted.values)
  return(est)
}

york.output <- york(x, y, weights.x = weights.x, weights.y = weights.y)

york.output$slope.after.each.iteration


## end of (relevant) script


#### Testing

