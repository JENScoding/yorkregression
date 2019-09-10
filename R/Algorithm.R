### york in Least Squares Fitting Of A Straight Line With Correlated Errors ###
## Input from Table I and Table II in york 1966
#setwd("/Users/jonascedrodelgado/Desktop/York-Regression/York/R")
'load("original_data.RData")
# here you can also load other data and weights
#' @title
#' Simple linear regression of X- and Y-variables with correlated errors.
#'
#' @description
#' Implements the algorithm for the problem of the best-ﬁt straight line to
#' independent points with errors in both x and y general York 1969 solution
#' according to the algorithm of Wehr & Saleska (2017)
#' @details
#' Given \eqn{n} pairs of \eqn{(X_i, Y_i), i = 1, \ldots, n}, their weights
#' \eqn{(\omega(X_i), \omega(Y_i)), i = 1, \ldots, n} or their standard errors
#' \eqn{SE(X_i)} and \eqn{SE(Y_i)}, the \code{york} function finds the best-fit
#' straight line using the algorithm of York et al. (1966)/ York et al. (1969)
#' as presented in Wehr & Saleska (2017). In addition, the function provides
#' numerous statistics, parameters and goodness of fit criteria. If the data
#' contains NA values then the share of NA values in the total values is
#' calculated and the rows with NA values will be deleted.
#' @param x A 1 times n numeric row vector of the \code{X}-variable
#' @param y A 1 times n numeric row vector of the \code{Y}-variable
#' @param tolerance The tolerance for convergence which is set a priori to
#' \code{1e-10}
#' @param weights.x The prespecified 1 times n weights vector for X-values
#' @param weights.y The prespecified 1 times n weights vector for Y-values
#' @param r.xy The prespecified correlation coefficient between X and Y
#' @param sd.x The standard error of the X-values
#' @param sd.y The standard error of the Y-values
#' @return York Returns an object of class "York" the York regression for the
#' \code{x} and \code{y} data for either specified weights \code{weights.x}
#' and \code{weights.y} or specified standard errors \code{sd.x} and \code{sd.y}
#' An object of class "York" containing the following components:
#'
#' \describe{
#'
#' \item{coefficients.york}{a matrix which contains the York estimates for
#' intercept and slope with the respective standard errors}
#' \item{coefficients.orthogonal}{a matrix which contains the Deming estimates
#' for intercept and slope with the respective standard errors}
#' \item{coefficients.ols}{a matrix which contains the OLS estimates for
#' intercept and slope with the respective standard errors}
#' \item{weighting.vector}{the prespecified or calculated weights}
#' \item{x.residuals}{the York X-residuals}
#' \item{y.residuals}{the York Y-residuals}
#' \item{fitted.y}{the fitted York Y-values}
#' \item{df.regression}{the number of degrees of freedom of the York regression}
#' \item{mean.x}{the ordinary mean of X}
#' \item{mean.y}{the ordinary mean of Y}
#' \item{reduced.chisq}{the reduced chi-squared statistic, i.e. the goodness
#' of fit measure of York's regression}
#' \item{std.Error.chisq}{the standard error of the chi-squared statistic}
#' \item{number.of.iterations}{the total number of iterations}
#' \item{slope.after.each.iteration}{the York slope after each iteration}
#' \item{original.x.values}{the original X-values}
#' \item{original.y.values}{the original Y-values}
#' \item{fitted.y.ols}{the fitted values for OLS}
#' \item{se.of.reg.ols}{the standard error of the regression for OLS}
#' \item{fitted.y.orthogonal}{the fitted values for Deming regression}
#' \item{data}{a data matrix which contains as columns the X-, Y-, sd.X- and
#' sd.Y-values}
#' }
#'
#' @references
#' Wehr, Richard, and Scott R. Saleska.
#' "The long-solved problem of the best-fit straight line: Application to
#'  isotopic mixing lines." Biogeosciences 14.1 (2017). pp. 17-29.
#'
#' York, Derek. "Least squares fitting of a straight line with correlated
#' errors." Earth and planetary science letters 5 (1968), pp. 320-324.
#'
#' York, Derek. "Least-squares fitting of a straight line.",
#' Canadian Journal of Physics 44.5 (1966), pp. 1079-1086.
#'
#' @examples
#' x <- c(0.0, 0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4)
#' y <- c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5)
#' weights.x <- c(1e+3, 1e+3, 5e+2, 8e+2, 2e+2, 8e+1, 6e+1, 2e+1, 1.8, 1)
#' weights.y <- c(1, 1.8, 4, 8, 20, 20, 70, 70, 1e+2, 5e+2)
#' r.xy <- 0
#' york(x = x, y = y, tolerance = 1e-10, weights.x = weights.x,
#' weights.y = weights.y, r.xy = r.xy, mult.samples = F)
#'
#' \dontrun{
#' # Example: No standard errors or weights specified
#' york(x, y, r.xy = 0)
#' # Example: You can't specify weights and standard errors at the same time
#' york(x , y, sd.x, sd.y, weights.x, weights.y, r.xy = 0)
#' # Example: x and y must have same length
#' york(x = c(0.0, 0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5),
#' y = c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5),
#' weights.x = weights.x, weights.y = weights.y, r.xy = 0)
#' }
#' @name york
#' @export
york <- function(x, y, tolerance = 1e-10, weights.x = NULL, weights.y = NULL,
                 sd.x = NULL, sd.y = NULL, r.xy = NULL, mult.samples = FALSE,
                 exact.solution = FALSE) {

  if (mult.samples == FALSE) {
    if (all(sapply(list(sd.x, sd.y, weights.x, weights.y), is.null))) {
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
    if (length(r.xy) == 1) {
      r.xy = rep(r.xy, length(x))
    }
    if (length(r.xy) != length(x)) {
      stop("Length of correlation vector must equal length of x")
    }
    #delete rows with NA values
    to.delete  <- c(which(is.na(x)), which(is.na(y)),
                    which(is.na(weights.x)), which(is.na(weights.y)),
                    which(is.na(sd.x)), which(is.na(sd.y)))
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
      sd.x <- sd.x[-to.delete]
      sd.y <- sd.y[-to.delete]
      r.xy <- r.xy[-to.delete]
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
      stop("sd.x and sd.y must have the same length of x resp. y!")
    }
    #initial value of the slope is the OLSE
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
    r.squared.ols <- 1 - sum((residuals)^2) / SS.y
  } else {
    if(exact.solution == T) {
      stop("There is no exact solution in case of multiple samples!")
    }
    # For Jonas
    x.input <- matrix(1)
    lm.ols <- matrix(1)
    fitted.y.ols <- matrix(1)
    residuals <- matrix(1)
    slope <- matrix(1)
    intercept.ols <- matrix(1)
    sigma.squared.hat <- matrix(1)
    se.of.reg.ols <- 1
    mean.x <- 1
    mean.y <- 1
    SS.x <- matrix(1)
    SS.y <- matrix(1)
    S.x <- matrix(1)
    SS.xy <- matrix(1)
    se.intercept.ols <- 1
    se.slope.ols <- 1
    r.squared.ols <- 1 - sum((residuals)^2) / SS.y

    lm.OLS <- list()
    slope <- NULL
    for (i in 1:5) {
      lm.OLS[[i]] <- lm(y[, i]~x[, i])
      slope[i] <- as.numeric(lm.OLS[[i]][[1]][2])
    }
    slope <- mean(slope)

    mean.xi <- apply(x, 1, mean)
    mean.yi <- apply(y, 1, mean)
    x.errors<- x - mean.xi
    y.errors <- y - mean.yi
    r.xy <- diag(cor(t(x.errors), t(y.errors)))
    weights.x <- apply(x, 1, var)
    weights.y <- apply(y, 1, var)
    x.original <- x
    y.original <- y
    x <- mean.xi
    y <- mean.yi
  }
  if (exact.solution == F) {
    slope.diff <- 10
    count <- 0
    slope.per.iteration <- NULL
    alpha <- sqrt(weights.x * weights.y)
    while (slope.diff > tolerance) {
      slope.old <- slope
      Weight <- alpha^2 / (slope^2 * weights.y + weights.x -
                             2 * slope * r.xy * alpha)
      Weight.sum <- sum(Weight)
      x.bar <- sum(Weight * x, na.rm = T) / Weight.sum
      y.bar <- sum(Weight * y, na.rm = T) / Weight.sum
      x.centered <- x - x.bar
      y.centered <- y - y.bar
      beta <- Weight * ((x.centered / weights.y) + (slope * y.centered /
                                                      weights.x) -
                          (slope * x.centered + y.centered) * r.xy / alpha)
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

  } else {
    if (any(r.xy != 0)) {
      stop("There is no exact solution in case of correlation between x and y errors!")
    }
    ## Apply formula and use lm estimate as intitial value for b
    alpha <- sqrt(weights.x * weights.y)
    Weight <- alpha^2 / (slope^2 * weights.y + weights.x)

    # see formula 19 and 20 for the following:
    Weight.sum <- sum(Weight)
    x.bar <- sum(Weight * x, na.rm = T) / Weight.sum
    y.bar <- sum(Weight * y, na.rm = T) / Weight.sum
    x.centered <- x - x.bar
    y.centered <- y - y.bar
    # calculate alpha.exact, beta.exact and gamma.exact. See York 66 page 1084
    alpha.exact <- 2 * sum(x.centered * y.centered * Weight^2 / weights.x) /
                  (3 * sum(x.centered^2 * Weight^2 / weights.x))
    beta.exact <- (sum(y.centered^2 * Weight^2 / weights.x) - sum(Weight * x.centered^2)) /
                  (3 * sum(x.centered^2 * Weight^2 / weights.x))
    gamma.exact <- - sum(x.centered * y.centered * Weight) / (sum(x.centered^2 * Weight^2 / weights.x))

    # use formula of York 66 to find slope, given on page 1084
    phi <- acos((alpha.exact^3 - 3 /2 * alpha.exact * beta.exact + 0.5 * gamma.exact) /
           (alpha.exact^2 - beta.exact)^(3 / 2))

    sol.cubic2 <- alpha.exact + 2 * (alpha.exact^2 - beta.exact)^0.5 * cos( 1 / 3 *(phi + 2 * pi * c(0:2)))
    sol.cubic2
    slope <- sol.cubic2[3]

    beta <- Weight * (x.centered / weights.y) + (slope * y.centered / weights.x)
    count <- 0
    slope.per.iteration <- data.frame("slope.per.iteration" =
                                        slope)
  }

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
  c <- r.xy*alpha
  x.residuals <- (Weight * (intercept + slope * x - y)
                  * (c - slope * weights.y)) /
    (weights.y * weights.x)
  y.residuals <- (Weight * (intercept + slope * x - y) *
                    (weights.x - slope * c)) / (weights.y * weights.x)
  #orthogonal regression

  slope.orthogonal <- (SS.y - SS.x + sqrt((SS.y - SS.x)^2 + 4*(SS.xy)^2)) /
    (2*SS.xy)
  intercept.orthogonal <- mean.y - slope.orthogonal * mean.x
  fitted.y.orthogonal <- intercept.orthogonal + slope.orthogonal * x
  r <- SS.xy / (sqrt(SS.x) * sqrt(SS.y))
  se.slope.orthogonal <- (slope.orthogonal/r) * sqrt((1 - r^2) / (length(x)))
  se.intercept.orthogonal <- ((1 / length(x)) * (sqrt(var(y)) - sqrt(var(x)) *
                                                   slope.orthogonal)^2 + (1 - r) * slope.orthogonal *
                                (2 * sqrt(var(x)) * sqrt(var(y)) +
                                   ((mean.x  *slope.orthogonal*(1+r)) / (r^2))))

  york.reg <- matrix(c(intercept, slope, sigma.intercept, sigma.slope),
                     nrow = 2)
  rownames(york.reg) <- c("intercept", "slope")
  colnames(york.reg) <- c("Estimate", "Std.Error")
  ols.reg <- matrix(c(intercept.ols, lm.ols[2], se.intercept.ols, se.slope.ols),
                    nrow = 2)
  rownames(ols.reg) <- c("intercept", "slope")
  colnames(ols.reg) <- c("Estimate", "Std.Error")
  orthogonal.reg <- matrix(c(intercept.orthogonal, slope.orthogonal,
                             se.intercept.orthogonal, se.slope.orthogonal), nrow = 2)
  rownames(orthogonal.reg) <- c("intercept", "slope")
  colnames(orthogonal.reg) <- c("Estimate", "Std.Error")
  if (mult.samples == F) {
    data <- matrix(c(x, y, sd.x, sd.y, r.xy), ncol = 5)
    colnames(data) <- c("x", "y", "sd.x", "sd.y", "r.xy")
  } else {
    data <- list("x" = x.original, "y" = y.original, "x.errors" = x.errors,
                 "y.errors" = y.errors, "r.xy" = r.xy)
  }

  output <- list("coefficients.york" = york.reg,
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
                  "fitted.y.ols" = fitted.y.ols,
                  "se.of.reg.ols" = se.of.reg.ols,
                  "r.squared.ols" = r.squared.ols,
                  "fitted.y.orthogonal" = fitted.y.orthogonal,
                  "data" = data)
  attr(output, "class") <- "york"

  return(output)
}

(york.output <- york(x, y, weights.x = weights.x, weights.y = weights.y, r.xy = 0, mult.samples = F))

class(york.output)

(york.output <- york(x, y, weights.x = weights.x, weights.y = weights.y, r.xy = 0, mult.samples = F, exact.solution = T))
