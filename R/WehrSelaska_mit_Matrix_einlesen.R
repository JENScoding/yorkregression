### york in Least Squares Fitting Of A Straight Line With Correlated Errors ###

## Input from Table I and Table II in york 1966
rm(list=ls())
x <- c(0,0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4)
y <- c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5)

#weights_y = c(1, 1.8, 4, 8, 20, 20, 70, 70, 1e+2, 5e+2)
#weights_x = c(1e+3, 1e+3, 5e+2, 8e+2, 2e+2, 8e+1, 6e+1, 2e+1, 1.8, 1)

sX <- x*0.01
sY <- y*0.005
n <- length(x)
rXY <- rep(0.8,n)
dat <- cbind(x,sX,y,sY,rXY)
dat
## function for algo

lm_OLS <- lm(y~x)
slope <- as.numeric(lm_OLS[[1]][2])

york <- function(x, tolerance = 1e-10){
  #initial value of b is OLS
  lm_OLS <- lm(x[,'y']~x[,'x'])
  slope <- as.numeric(lm_OLS[[1]][2])

  slope_diff <- 10
  count <- 0

  while (slope_diff > tolerance) {

    slope_old <- slope
    # omega, Jonas hat aber gesagt: "Lass es weg"
    weights_x <- 1 / x[,'sX']^2
    weights_y <- 1 / x[,'sY']^2
    alpha <- sqrt(weights_x * weights_y)
    Weight <- alpha^2 / (slope^2 * weights_y + weights_x - 2 * slope * x[,'rXY'] * alpha)
    Weight_sum <- sum(Weight)
    x_bar <- sum(Weight * x[,'x']) / Weight_sum
    y_bar <- sum(Weight * x[,'y']) / Weight_sum
    x_centered <- x[,'x'] - x_bar
    y_centered <- x[,'y'] - y_bar

    beta <- Weight * ((x_centered / weights_y) + (slope * y_centered / weights_x) - (slope * x_centered + y_centered) * x[,'rXY'] / alpha)
    Q1 <- sum(Weight * beta * y_centered)
    Q2 <- sum(Weight * beta * x_centered)
    slope <- Q1 / Q2
    slope_diff <- abs(slope - slope_old)
    count <- count + 1
    print(slope, digits = 10)

  }

  intercept <- y_bar - slope * x_bar

  x_mean <- 0
  x_adj <- rep(0,length(x[,'x']))

  x_adj <- x_bar + beta
  x_mean <- sum(Weight * beta)


  x_mean <- x_mean / (Weight_sum * (length(x[,'x']) - 2))

  sigma_slope <- 0
  chisq_weight <- 0
  u <- rep(0,length(x[,'x']))

  for(i in 1:length(x[,'x'])){
    u[i] <- x_adj[i] - x_mean
    sigma_slope <- sigma_slope + Weight[i]*u[i]^2
    chisq_weight <- chisq_weight + Weight[i]*(x[,'y'][i]-slope*x[,'x'][i]-intercept)^2
  }

  sigma_slope <- sqrt(1 / sigma_slope)
  sigma_intercept <- sqrt(x_mean^2 * sigma_slope^2 + 1 / Weight_sum)
  chisq_weight <- chisq_weight / (length(x[,'x']) - 2)
  sigma_x <- sqrt(2 / (length(x[,'x']) - 2))

  fitted_y <- intercept + slope * x[,'x']

  residuals <- x[,'y'] - fitted_y

  c <- x[,'rXY']*alpha

  x_resdiduals <- (Weight * (intercept + slope * x[,'x'] - x[,'y']) * (c - slope * weights_y)) / (weights_y * weights_x)

  y_resdiduals = (Weight * (intercept + slope * x[,'x'] - x[,'y']) * (weights_x - slope * c))/ (weights_y * weights_x)



  mt <- matrix(c(intercept, slope, sigma_intercept, sigma_slope), nrow = 2)
  rownames(mt) <- c("intercept", "slope")
  colnames(mt) <- c("Estimate", "Std.Error")
  est <- list("coefficients" = mt, "weighting.vector" = Weight, "x.residuals" = x_resdiduals,"y.residuals"=y_resdiduals,"fitted_y"=fitted_y, "mean.x" = x_bar, "mean.y" = y_bar ,"Std.Error.chi" = sigma_x, "iterations" = count, x_centered, y_centered, x, y, x_mean, "show" = x_adj)
  return(est)
}

first <- york(x = dat)
first

plot(dat[,'x'],y, pch = 16)
lines(dat[,'x'],first$fitted_y , col = "red",lwd = 2)
lines(dat[,'x'], lm_OLS$fitted.values, col = "blue", lty = "dashed",lwd = 2)
legend("topright",legend = c("OLS","York"), fill = c("blue","red"))
#the York regression line and the OLS regression line both go to the "center of gravity"
abline(v= first$mean.x,h= first$mean.y ,lty = "dashed", col = "red")
points(x= first$mean.x,y= first$mean.y,col="red", pch = 16)
abline(v= mean(dat[,'x']),h= mean(dat[,'y']) ,lty = "dashed", col = "blue")
points(x= mean(dat[,'x']),y= mean(dat[,'y']),col="blue", pch = 16)

plot(first$x.residuals,first$y.residuals)

## end of (relevant) script

