### york in Least Squares Fitting Of A Straight Line With Correlated Errors ###
## Input from Table I and Table II in york 1966
#setwd("/Users/jonascedrodelgado/Desktop/York-Regression/York/R")
load("original_data.RData")
# here you can also load other data and weights

#' Simple linear regression of X and Y-variables with correlated errors.
#' @description
#' Implements the algorithm for the problem of the best-ﬁt straight line to
#' independent points with errors in both x and y general York 1969 solution
#' according to the algorithm of Wehr & Saleska (2017)
#' @details
#' Given \eqn{n} pairs of \eqn{(X_i, Y_i), i = 1, \ldots, n}, their weights
#' \eqn{(\omega(X_i), \omega(Y_i)), i = 1, \ldots, n} or their standard errors
#' \eqn{SE(X_i)} and \eqn{SE(Y_i)}, the \code{york} function finds the best fit
#' straight line using the algorithm of York et al. (1966)/ York et al. (1969)
#' as presented in Wehr & Saleska (2017). If the data contains NA values then
#' the share of NA values in the total values is calculated und the rows with NA
#'  values will be deleted.


#' @param x A 1 times n numeric row vector.of the \code{X} Variable
#' @param y A 1 times n numeric row vector.of the \code{Y} Variable
#' @param tolerance The tolerance for convergence which is set a priori to
#' \code{1e-10}
#' @param weights.x The prespecified 1 times n weights vector for x values
#' @param weights.y The prespecified 1 times n weights vector for y values
#' @param xy.error.correlation The prespecified correlation coefficient between X and Y
#' @param x.errors The standard error of the x values
#' @param y.errors The standard error of the y values
#' @return yorks Returns an object of class "york" the York Regression for the
#' \code{x} and \code{y} data for either specified weights \code{weights.x}
#' and \code{weights.y} or specified standard errors \code{x.errors} abd \code{y.errors}
#' An object of class "york" containing the following components:
#' coefficients.york (a matrix which contains the york estimates for intercept
#' and slope with the respective standard errors), coefficients.orthogonal
#' (a matrix which contains the deming estimates for intercept and slope with
#' the respective standard errors), coefficients.ols (a matrix which contains
#' the ols estimates for intercept and slope with the respective
#' standard errors), weighting.vector (the speciefied or calculated weights),
#' x.residuals (the York x-residuals), y.residuals (the York y-residuals),
#' fitted.y (the fitted York y-values), df.regression (the number of degress
#' of freedom of the York regression), mean.x (the ordinary mean of X), mean.y
#' (the ordinary mean of Y), reduced.chisq (the reduced chi-squared statistic
#' i.e. the goodness of fit measure of York's regression), std.Error.chisq
#' (the standard error of the chi-squared statistic), number.of.iterations
#' (the total number of iterations), slope.after.each.iteration (the york slope
#' after each interation), original.x.values (the original x values),
#' original.y.values (the original y values), fitted.y.ols (the fitted values
#' for OLS), se.of.reg.ols (the standard error of the regression for OLS),
#' fitted.y.orthogonal (the fitted values for deming regression), data
#' (a data matrix which contains as colums the x-, y-,x.errors- and y.errors-values)
#' @examples
#' york(x, y, tolerance = 1e-10, weights.x = NULL, weights.y = NULL,
#'xy.error.correlation = NULL, x.errors = NULL, y.errors = NULL)
#' \dontrun{
#' # Example: No standard errors or weights specified
#' york(x, y, xy.error.correlation = 0)
#' # Example: You can't specify weights and standard errors at the same time
#' york(x , y, x.errors, y.errors, weights.x, weights.y, xy.error.correlation = 0)
#' # Example: x and y must have same length
#' york(x = c(0.0, 0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5),
#' y = c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5),
#' weights.x = weights.x, weights.y = weights.y, xy.error.correlation = 0)
#' }
#' @name york
#' @export
?lm

help(York)

york <- function(x, y, tolerance = 1e-10, weights.x = NULL, weights.y = NULL,
                 xy.error.correlation = NULL, x.errors = NULL, y.errors = NULL,
                 mult.samples = FALSE) {

  if (mult.samples == FALSE) {
    if (all(sapply(list(x.errors, y.errors, weights.x, weights.y), is.null))) {
      stop("Specify either standard errors or weights")
    }
    if (all(sapply(list(x.errors, y.errors, weights.x, weights.y),
                   function(x) !is.null(x)))) {
      stop("You can't specify weights and standard errors at the same time!")
    }
    if (any(is.na(x))) {
      stop("There is at least one NA value. Please specify this value(s)!")
    }
    if (length(x.errors) == 1) {
      x.errors = rep(x.errors, length(x))
    }
    if (length(y.errors) == 1) {
      y.errors = rep(y.errors, length(y))
    }
    if(length(x) != length(y)) {
      stop("x and y must have same length!")
    }
    if (is.null(c(x.errors, y.errors, weights.x, weights.y))) {
      stop("Specify either standard errors or weights")
    }
    if (all(sapply(list(x.errors, y.errors, weights.x, weights.y),
                   function(x) !is.null(x)))) {
      stop("You can't specify weights and standard errors at the same time!")
    }
    if (length(x.errors) == 1) {
      x.errors = rep(x.errors, length(x))
    }
    if (length(y.errors) == 1) {
      y.errors = rep(y.errors, length(y))
    }
    if(length(x) != length(y)) {
      stop("x and y must have same length!")
    }
    if (length(xy.error.correlation) == 1) {
      xy.error.correlation = rep(xy.error.correlation, length(x))
    }
    if (length(xy.error.correlation) != length(x)) {
      stop("Length of correlation vector must equal length of x")
    }
    #delete rows with NA values
    to.delete  <- c(which(is.na(x)), which(is.na(y)),
                    which(is.na(weights.x)), which(is.na(weights.y)),
                    which(is.na(x.errors)), which(is.na(y.errors)))
    rm.share <- length(to.delete) / length(x)
    if (rm.share > 0.1) {
      warning(rm.share * 100,
              "% of the data were removed due to missing values!")
    }
    if (length(to.delete) > 0){
      y <- y[-to.delete]
      x <- x[-to.delete]
      weights.x <- weights.x[-to.delete]
      weights.y <- weights.y[-to.delete]
      x.errors <- x.errors[-to.delete]
      y.errors <- y.errors[-to.delete]
      xy.error.correlation <- xy.error.correlation[-to.delete]
    }
    if (is.null(weights.x) & is.null(weights.y)) {
      weights.x <- 1/x.errors^2
      weights.y <- 1/y.errors^2
    }
    if (is.null(x.errors) & is.null(y.errors)) {
      x.errors <- 1/ sqrt(weights.x)
      y.errors <- 1/sqrt(weights.y)
    }
    if(length(x.errors) != length(x) | length(y.errors) != length(y)) {
      stop("x.errors and y.errors must have the same length of x resp. y!")
    }
    #initial value of b is OLS
    x.input <- matrix(c(rep(1, length(x)), x), ncol =2)
    lm.ols <- solve(t(x.input) %*% x.input) %*% t(x.input) %*% y
    fitted.y.ols <- x.input %*% lm.ols
    residuals <- y - fitted.y.ols
    slope <- as.numeric(lm.ols[2])
    intercept.ols <- as.numeric(lm.ols[1])
    if (any(is.na(c(slope,intercept.ols)))){
      stop("Cannot fit a line through these data!")
    }
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
    # For Jonas
    x.input <- 1
    lm.ols <- 1
    fitted.y.ols <- 1
    residuals <- 1
    slope <- 1
    intercept.ols <- 1
    sigma.squared.hat <- 1
    se.of.reg.ols <- 1
    mean.x <- 1
    mean.y <- 1
    SS.x <- 1
    SS.y <- 1
    S.x <- 1
    SS.xy <- 1
    se.intercept.ols <- 1
    se.slope.ols <- 1

    lm.OLS <- list()
    slope <- NULL
    for (i in 1:5) {
      lm.OLS[[i]] <- lm(y[, i]~x[, i])
      slope[i] <- as.numeric(lm.OLS[[i]][[1]][2])
    }
    slope <- mean(slope)

    mean.xi <- apply(x, 1, mean)
    mean.yi <- apply(y, 1, mean)
    x.errors <- x - mean.xi
    y.errors <- y - mean.yi
    xy.error.correlation <- diag(cor(t(x.errors), t(y.errors)))
    weights.x <- apply(x, 1, var)
    weights.y <- apply(y, 1, var)
    x.original <- x
    y.original <- y
    x <- mean.xi
    y <- mean.yi
  }

  slope.diff <- 10
  count <- 0
  slope.per.iteration <- NULL
  alpha <- sqrt(weights.x * weights.y)
  while (slope.diff > tolerance) {
    slope.old <- slope
    Weight <- alpha^2 / (slope^2 * weights.y + weights.x -
                           2 * slope * xy.error.correlation * alpha)
    Weight.sum <- sum(Weight)
    x.bar <- sum(Weight * x, na.rm = T) / Weight.sum
    y.bar <- sum(Weight * y, na.rm = T) / Weight.sum
    x.centered <- x - x.bar
    y.centered <- y - y.bar
    beta <- Weight * ((x.centered / weights.y) + (slope * y.centered /
                                                    weights.x) -
                        (slope * x.centered + y.centered) * xy.error.correlation / alpha)
    Q1 <- sum(Weight * beta * y.centered, na.rm = T)
    Q2 <- sum(Weight * beta * x.centered, na.rm = T)
    slope <- Q1 / Q2
    slope.diff <- abs(slope - slope.old)
    count <- count + 1
    slope.per.iteration <- append(slope.per.iteration, slope)
    if (count > tolerance^-1) {
      stop("\nThe slope coefficient does not converge after ",
           count," iterations. \nHint: You may reduce the tolerance level.",
           cat("Slope coefficient for last 5 iterations:"),
           for (i in 4:0){
             cat("\n\t", count - i, "\t", slope.per.iteration[count - i])},
           cat("\n"))
    }
  }
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
  c <- xy.error.correlation*alpha
  x.residuals <- (Weight * (intercept + slope * x - y)
                  * (c - slope * weights.y)) /
    (weights.y * weights.x)
  y.residuals <- (Weight * (intercept + slope * x - y) *
                    (weights.x - slope * c)) / (weights.y * weights.x)
  #orthogonal regression

  slope.orthogonal <- (SS.y - SS.x + sqrt((SS.y - SS.x)^2 + 4*(SS.xy)^2)) / (2*SS.xy)
  intercept.orthogonal <- mean.y - slope.orthogonal * mean.x
  fitted.y.orthogonal <- intercept.orthogonal + slope.orthogonal * x
  r <- SS.xy / (sqrt(SS.x) * sqrt(SS.y))
  se.slope.orthogonal <- (slope.orthogonal/r) * sqrt((1 - r^2) / (length(x)))
  se.intercept.orthogonal <- ((1 / length(x)) * (sqrt(var(y)) - sqrt(var(x)) * slope.orthogonal)^2 +
                           (1 - r) * slope.orthogonal * (2 * sqrt(var(x)) * sqrt(var(y)) +
                                                      ((mean.x  *slope.orthogonal*(1+r)) / (r^2))))

  york.reg <- matrix(c(intercept, slope, sigma.intercept, sigma.slope), nrow = 2)
  rownames(york.reg) <- c("intercept", "slope")
  colnames(york.reg) <- c("Estimate", "Std.Error")
  ols.reg <- matrix(c(intercept.ols, lm.ols[2], se.intercept.ols, se.slope.ols), nrow = 2)
  rownames(ols.reg) <- c("intercept", "slope")
  colnames(ols.reg) <- c("Estimate", "Std.Error")
  orthogonal.reg <- matrix(c(intercept.orthogonal, slope.orthogonal, se.intercept.orthogonal, se.slope.orthogonal), nrow = 2)
  rownames(orthogonal.reg) <- c("intercept", "slope")
  colnames(orthogonal.reg) <- c("Estimate", "Std.Error")
  if (mult.samples == T) {
    data <- list("x" = x.original, "y" = y.original, "x.errors" = x.errors, "y.errors" = y.errors,
                 "xy.error.correlation" = xy.error.correlation)
  } else {
    data <- matrix(c(x, y, x.errors, y.errors, xy.error.correlation), ncol = 5)
    colnames(data) <- c("x", "y", "x.errors", "y.errors", "xy.error.correlation")
  }
  est <- list("coefficients.york" = york.reg,
              "coefficients.orthogonal" = orthogonal.reg,
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


(york.output <- york(x, y, weights.x = weights.x, weights.y = weights.y, xy.error.correlation = 0, mult.samples = F))
