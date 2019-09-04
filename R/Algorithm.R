### york in Least Squares Fitting Of A Straight Line With Correlated Errors ###
## Input from Table I and Table II in york 1966
#setwd("/Users/jonascedrodelgado/Desktop/York-Regression/York/R")
load("R/original_data.RData")
# here you can also load other data and weights

## function for algo
york <- function(x, y, tolerance = 1e-10, weights.x = NULL, weights.y = NULL,
                 rxy = NULL, sd.x = NULL, sd.y = NULL, mult.samples = FALSE) {
  if (mult.samples == FALSE) {
    if (is.null(c(sd.x, sd.y, weights.x, weights.y))) {
      stop("Specify either standard errors or weights")
    }
    if (all(sapply(list(sd.x, sd.y, weights.x, weights.y),
                   function(x) !is.null(x)))) {
      stop("You can't specify weights and standard errors at the same time!")
    }
    if (length(sd.x) == 1) {
      sd.x = rep(sd.x, length(x))
    }
    if (length(sd.y) == 1) {
      sd.y = rep(sd.y, length(y))
    }
    if(length(x) != length(y)) {
      stop("x and y must have same length!")
    }
    if (length(rxy) == 1) {
      rxy = rep(rxy, length(x))
    }
    if (length(rxy) != length(x)) {
      stop("Length of correlation vector must equal length of x")
    }
    #delete rows with NA values
    to.delete  <- c(which(is.na(x)), which(is.na(y)),
                    which(is.na(weights.x)), which(is.na(weights.y)),
                    which(is.na(sd.x)), which(is.na(sd.y)))
    rm.share <- length(to.delete) / length(x)
    if (rm.share > 0.1) {
      warning(rm.share * 100, "% of the data were removed due to missing values!")
    }
    if (length(to.delete) > 0){
      y <- y[-to.delete]
      x <- x[-to.delete]
      weights.x <- weights.x[-to.delete]
      weights.y <- weights.y[-to.delete]
      sd.x <- sd.x[-to.delete]
      sd.y <- sd.y[-to.delete]
      rxy <- rxy[-to.delete]
    }
    if (is.null(weights.x) & is.null(weights.y)) {
      weights.x <- 1/sd.x^2
      weights.y <- 1/sd.y^2
    }
    if (is.null(sd.x) & is.null(sd.y)) {
      sd.x <- 1/ sqrt(weights.x)
      sd.y <- 1/sqrt(weights.y)
    }
    if(length(sd.x) != length(x) | length(sd.y) != length(y)) {
      stop("Sd.x and sd.y must have the same length of x resp. y")
    }
    #initial value of b is OLS
    x.input <- matrix(c(rep(1, length(x)), x), ncol =2)
    lm.ols <- solve(t(x.input) %*% x.input) %*% t(x.input) %*% y
    fitted.y.ols <- x.input %*% lm.ols
    residuals <- y - fitted.y.ols
    slope <- as.numeric(lm.ols[2])
    intercept.ols <- as.numeric(lm.ols[1])
    sigma.squared.hat <- (1 / (length(x) - 2)) * sum((residuals)^2)
    se.of.reg.ols <- sqrt(sigma.squared.hat)
    mean.x <- mean(x)
    mean.y <- mean(y)
    centered.x <- x - mean.x
    centered.y <- y - mean.y
    SS.x <- sum(centered.x^2)
    SS.y <- sum(centered.y^2)
    S.x <- sum(x^2)
    SS.xy <- sum((centered.x) * (centered.y))
    se.intercept.ols <- sqrt(sigma.squared.hat * (S.x / (length(x) * SS.x)))
    se.slope.ols <- sqrt(sigma.squared.hat / SS.x)
  } else {
    stop("Not ready yet")
  }


  slope.diff <- 10
  count <- 0
  slope.per.iteration <- NULL
  while (slope.diff > tolerance) {
    slope.old <- slope
    alpha <- sqrt(weights.x * weights.y)
    Weight <- alpha^2 / (slope^2 * weights.y + weights.x -
                           2 * slope * rxy * alpha)
    Weight.sum <- sum(Weight)
    x.bar <- sum(Weight * x, na.rm = T) / Weight.sum
    y.bar <- sum(Weight * y, na.rm = T) / Weight.sum
    x.centered <- x - x.bar
    y.centered <- y - y.bar
    beta <- Weight * ((x.centered / weights.y) + (slope * y.centered /
                                                    weights.x) -
                        (slope * x.centered + y.centered) * rxy / alpha)
    Q1 <- sum(Weight * beta * y.centered, na.rm = T)
    Q2 <- sum(Weight * beta * x.centered, na.rm = T)
    slope <- Q1 / Q2
    slope.diff <- abs(slope - slope.old)
    count <- count + 1
    slope.per.iteration <- append(slope.per.iteration, slope)
    if (count > tolerance^-1)
      stop("\nThe slope coefficient does not converge after ",
           count," iterations. \nHint: You may reduce the tolerance level.",
           cat("Slope coefficient for last 5 iterations:"),
           for (i in 4:0){
             cat("\n\t", count - i, "\t", slope.per.iteration[count - i])},
           cat("\n"))}
  slope.per.iteration <- data.frame("slope.per.iteration" =
                                      slope.per.iteration)
  intercept <- y.bar - slope * x.bar
  x.adj <- x.bar + beta
  x.mean <- sum(Weight * beta, na.rm = T) / (Weight.sum * (length(x) - 2))
  u <- x.adj - x.mean
  sigma.slope <- sqrt(1 / sum(Weight * u^2, na.rm = T))
  sigma.intercept <- sqrt(x.mean^2 * sigma.slope^2 + 1 / Weight.sum)
  sigma.slope.intercept <- -x.mean*sigma.slope^2
  reduced.chisq <- sum(Weight * (y - slope * x - intercept)^2, na.rm = T) /
    (length(x) - 2)
  sigma.chisq <- sqrt(2 / (length(x) - 2))
  fitted.y <- intercept + slope * x
  residuals <- y - fitted.y
  df.regression <- 2*(length(x)-1)
  c <- rxy*alpha
  x.residuals <- (Weight * (intercept + slope * x - y)
                  * (c - slope * weights.y)) /
    (weights.y * weights.x)
  y.residuals <- (Weight * (intercept + slope * x - y) *
                    (weights.x - slope * c)) / (weights.y * weights.x)
  #total least squares/ simple major axis regression

  slope.mayor <- (SS.y - SS.x + sqrt((SS.y - SS.x)^2 + 4*(SS.xy)^2)) / (2*SS.xy)
  intercept.mayor <- mean.y - slope.mayor * mean.x
  fitted.y.orthogonal <- intercept.mayor + slope.mayor * x
  r <- SS.xy / (sqrt(SS.x) * sqrt(SS.y))
  se.slope.mayor <- (slope.mayor/r) * sqrt((1 - r^2) / (length(x)))
  se.intercept.mayor <- ((1 / length(x)) * (sqrt(var(y)) - sqrt(var(x)) * slope.mayor)^2 +
                   (1 - r) * slope.mayor * (2 * sqrt(var(x)) * sqrt(var(y)) +
                                          ((mean.x  *slope.mayor*(1+r)) / (r^2))))

  york.reg <- matrix(c(intercept, slope, sigma.intercept, sigma.slope), nrow = 2)
  rownames(york.reg) <- c("intercept", "slope")
  colnames(york.reg) <- c("Estimate", "Std.Error")
  ols.reg <- matrix(c(intercept.ols, lm.ols[2], se.intercept.ols, se.slope.ols), nrow = 2)
  rownames(ols.reg) <- c("intercept", "slope")
  colnames(ols.reg) <- c("Estimate", "Std.Error")
  mayor.reg <- matrix(c(intercept.mayor, slope.mayor, se.intercept.mayor, se.slope.mayor), nrow = 2)
  rownames(mayor.reg) <- c("intercept", "slope")
  colnames(mayor.reg) <- c("Estimate", "Std.Error")
  data <- matrix(c(x, y, sd.x, sd.y, rxy), ncol = 5)
  colnames(data) <- c("x", "y", "sd.x", "sd.y", "rxy")
  est <- list("coefficients.york" = york.reg,
              "coefficients.mayor" = mayor.reg,
              "coefficients.ols" = ols.reg,
              "weighting.vector" = Weight,
              "x.residuals" = x.residuals,
              "y.residuals"= y.residuals,
              "fitted.y"=fitted.y,
              "df.regression" = df.regression,
              "mean.x" = x.bar,
              "mean.y" = y.bar ,
              "reduced.chisq" = reduced.chisq,
              "std.Error.chisq" = sigma.chisq,
              "number.of.iterations" = count,
              "slope.after.each.iteration" = slope.per.iteration,
              x.centered, y.centered, x, y, x.mean, "show" = x.adj,
              "original.x.values" = x,
              "original.y.values" = y,
              "fitted.y.ols" = fitted.y.ols,
              "se.of.reg.ols" = se.of.reg.ols,
              "fitted.y.orthogonal" = fitted.y.orthogonal,
              "data" = data)
  return(est)
}

(york.output <- york(x, y, weights.x = weights.x, weights.y = weights.y, rxy = 0, mult.samples = T))
#york.output$slope.after.each.iteration

#york.plots(x = x, y = y, rxy = 0.3)

## end of (relevant) script

#### Testing

