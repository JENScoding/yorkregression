### york in Least Squares Fitting Of A Straight Line With Correlated Errors ###
#setwd()
load("multiple_samples.RData")


## function for algo

york <- function(x, y, tolerance = 1e-10, weights.x, weights.y){
  #initial value of b is OLS
  lm.OLS <- list()
  slope <- NULL
  for (i in 1:5){
    lm.OLS[[i]] <- lm(y[, i]~x[, i])
    slope[i] <- as.numeric(lm.OLS[[i]][[1]][2])
  }
  slope <- mean(slope)
  
  slope.diff <- 10
  count <- 0
  slope.per.iteration <- NULL
  
  while (slope.diff > tolerance) {
    slope.old <- slope
    x.bar <- 0
    y.bar <- 0
    W.sum <- 0
    alpha <- rep(0,length(x))
    Weight <- rep(0,length(x))
    omega.x <- rep(0,length(x))
    omega.y <- rep(0,length(x))
    error.correlation <- rep(0,length(x))
    for (i in 1:length(x)) {
      omega.x[i] <- 1 / var(x[i, ])
      omega.y[i] <- 1 / var(y[i, ])
      alpha[i] <- sqrt(omega.x[i] * omega.y[i])
      Weight[i] <- alpha[i]^2 / (slope^2 * omega.y[i] + omega.x[i] - 2 * slope * error.correlation[i] * alpha[i])
      x.bar <- x.bar + Weight[i] * x[i, ]
      y.bar <- y.bar + Weight[i] * y[i, ]
      W.sum <- W.sum + Weight[i]
    }
    X_bar <- X_bar / W_sum
    Y_bar <- Y_bar / W_sum
    
    Q1 <- 0
    Q2 <- 0
    U <- rep(0,length(x))
    V <- rep(0,length(x))
    beta <- rep(0,length(x))
    for (i in 1:length(x)) {
      U[i] <- x[i] - X_bar
      V[i] <- y[i] - Y_bar
      beta[i] <- W[i] * ((U[i] / weights[i, 1]) + (b * V[i] / weights[i, 2]) - (b * U[i] + V[i]) * 0 / alpha[i])
      Q1 <- Q1 + W[i] * beta[i] * V[i]
      Q2 <- Q2 + W[i] * beta[i] * U[i]
    }
    #Q1 <- sum(W*beta*V)
    #Q2 <- sum(W*beta*U)
    b <- round(Q1 / Q2, 100)
    b_diff <- abs(b - b_old)
    count <- count + 1
    
  }
    
    slope.old <- slope
    omega, Jonas hat aber gesagt: "Lass es weg"
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

