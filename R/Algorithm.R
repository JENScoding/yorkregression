### york in Least Squares Fitting Of A Straight Line With Correlated Errors ###

#' @title Simple linear regression of X- and Y-variables with correlated errors.
#'
#' @description Implements the algorithm for the problem of the best-ﬁt straight
#' line to independent points with errors in both x and y general York 1969
#' solution according to the algorithm of Wehr & Saleska (2017)
#' @details Given \eqn{n} pairs of \eqn{(X_i, Y_i), i = 1, \ldots, n}, their
#' weights \eqn{(\omega(X_i), \omega(Y_i)), i = 1, \ldots, n} or their standard
#' errors \eqn{SE(X_i)} and \eqn{SE(Y_i)}, the \code{york} function finds the
#' best-fit straight line using the algorithm of York et al. (1966)/ York et al.
#' (1969) as presented in Wehr & Saleska (2017). In addition, the function
#' provides numerous statistics, parameters and goodness of fit criteria. If the
#' data contains NA values then the share of \code{NA} values in the total
#' values is calculated and the rows with \code{NA} values will be deleted.
#' @param x A 1 times n numeric row vector or a dataframe of the \code{X}-variable
#' @param y A 1 times n numeric row vector or a dataframe of the \code{Y}-variable
#' @param tolerance The tolerance for convergence of the slope coefficent.
#' It is set a priori to \code{1e-5}
#'@param max.iterations The maximum number of iterations for convergence.
#'   Default is \code{50}.
#' @param weights.x The prespecified 1 times n weights vector for
#'   \code{X}-values
#' @param weights.y The prespecified 1 times n weights vector for
#'   \code{Y}-values
#' @param sd.x The standard error of the \code{X}-values
#' @param sd.y The standard error of the \code{Y}-values
#' @param r.xy The prespecified correlation coefficient between the errors in
#'   \code{X} and \code{Y}
#' @param mult.samples An indicator if the multiple samples option is turned on
#' or not. The standard value is \code{FALSE}
#' @param approx.solution An indicator if the approximate solution option is turned on
#' or not. The standard value is \code{FALSE}
#' @return York Returns an object of class "york" the York regression for the
#'   \code{x} and \code{y} data for either specified weights \code{weights.x}
#'   and \code{weights.y} or specified standard errors \code{sd.x} and
#'   \code{sd.y} An object of class "york" containing the following components:
#'
#'   \describe{
#'
#'   \item{coefficients}{a matrix which contains the York estimates for
#'   intercept and slope of the best-fit straight line with the respective
#'   standard errors}
#'   \item{coefficients.orthogonal}{a matrix which contains the
#'   orthogonal estimates for intercept and slope with the respective standard
#'   errors}
#'   \item{coefficients.ols}{a matrix which contains the OLS estimates
#'   for the intercept and slope with the respective standard errors}
#'   \item{weights}{a matrix representation of the prespecified or calculated
#'   weights for the X- and Y-observations}
#'   \item{x.residuals}{a vector of the
#'   York X-residuals}
#'   \item{y.residuals}{a vector of the York Y-residuals}
#'   \item{fitted.y}{a vector of the fitted York Y-values}
#'   \item{df.regression}{the number of degrees of freedom
#'   (See \url{https://en.wikipedia.org/wiki/Degrees_of_freedom_(statistics)})
#'   of York's regression}
#'   \item{weighted.mean.x}{the weighted.mean of X}
#'   \item{weighted.mean.y}{the weighted.mean of Y}
#'   \item{reduced.chisq}{the reduced chi-squared statistic
#'   (See \url{https://en.wikipedia.org/wiki/Reduced_chi-squared_statistic}),
#'   i.e. the goodness of fit measure of York's regression}
#'   \item{std.error.chisq}{the standard error of the chi-squared statistic}
#'   \item{number.of.iterations}{the total number of iterations}
#'   \item{slope.after.each.iteration}{the York slope after each iteration}
#'   \item{fitted.y.ols}{the fitted values for OLS}
#'   \item{residuals.ols}{the OLS residuals}
#'   \item{residual.sum.of.squares}{the residual sum of squares (RSS) of OLS}
#'   \item{total.sum.of.squares}{the total sum of squares of OLS}
#'   \item{se.of.reg.ols}{the standard error of the regression (SER) for OLS}
#'   \item{r.squared.ols}{the R squared of the OLS regression}
#'   \item{r.squared.adjusted.ols}{the adjusted R squared of OLS}
#'   \item{f.statistic.ols}{the F statistic of OLS}
#'   \item{fitted.y.orthogonal}{the fitted values for orthogonal regression (See
#'   \url{https://en.wikipedia.org/wiki/Deming_regression})}
#'   \item{data}{a data matrix which contains as columns the observed points
#'   X-, Y-, sd.X- and sd.Y-values}}
#'
#' @references Wehr, Richard, and Scott R. Saleska. "The long-solved problem of
#' the best-fit straight line: Application to isotopic mixing lines."
#' Biogeosciences 14.1 (2017). pp. 17-29.
#'
#' York, Derek. "Least squares fitting of a straight line with correlated
#' errors." Earth and planetary science letters 5 (1968), pp. 320-324.
#'
#' York, Derek. "Least-squares fitting of a straight line.", Canadian Journal of
#' Physics 44.5 (1966), pp. 1079-1086.
#'
#' York, Derek, et al. "Unified equations for the slope, intercept, and standard
#' errors of the best straight line." American Journal of Physics 72.3 (2004),
#' pp. 367-375.
#'
#' @examples
#' # Example: York's regression with weight data taken from Pearson (1901):
#' x <- c(0.0, 0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4)
#' y <- c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5)
#' weights.x <- c(1e+3, 1e+3, 5e+2, 8e+2, 2e+2, 8e+1, 6e+1, 2e+1, 1.8, 1)
#' weights.y <- c(1, 1.8, 4, 8, 20, 20, 70, 70, 1e+2, 5e+2)
#' r.xy <- 0
#' york(x = x, y = y, tolerance = 1e-10, weights.x = weights.x,
#' weights.y = weights.y, r.xy = r.xy)
#'
#' # Example: York's regression arbitrary values for sd.x and sd.y:
#' x <- c(0.0, 0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4)
#' y <- c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5)
#' sd.x <- 0.2
#' sd.y <- 0.4
#' r.xy <- 0.3
#' york(x = x, y = y, tolerance = 1e-10, sd.x = sd.x,
#' sd.y = sd.y, r.xy = r.xy)
#'
#' \dontrun{
#' # Example: No standard errors or weights specified
#' york(x, y, r.xy = 0)
#' # Example: You can't specify weights and standard errors at the same time
#' york(x , y, sd.x, sd.y, weights.x, weights.y, r.xy = 0)
#' # Example: x and y must have same length
#' york(x = c(0.0, 0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5),
#' y = c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5),
#' weights.x = weights.x, weights.y = weights.y,  r.xy = 0)
#' }
#' @name york
#' @export
#' @importFrom stats pchisq
#' @importFrom utils stack
york <- function(x, y, weights.x = NULL, weights.y = NULL, tolerance = 1e-5,
                 max.iterations = 50, sd.x = NULL, sd.y = NULL, r.xy = NULL,
                 mult.samples = FALSE, approx.solution = FALSE) {

  if (mult.samples == FALSE) {

    # rewrite input and delete rows with NA values
    input <- f_rewrite(x, y, weights.x = weights.x, weights.y = weights.y,
                     sd.x = sd.x, sd.y = sd.y, r.xy = r.xy)
    x <- input$x
    y <- input$y
    weights.x <- input$weights.x
    weights.y <- input$weights.y
    sd.x <- input$sd.x
    sd.y <- input$sd.y
    r.xy <- input$r.xy

    # expected errors for wrongly specified input
    exp_error_simple(x, y, weights.x = weights.x, weights.y = weights.y,
                       sd.x = sd.x, sd.y = sd.y, r.xy = r.xy)

  } else { # mult.sample = TRUE

    # expected errors for wrongly specified multiple sample input
    exp_error_multiple(x, y, weights.x = weights.x, weights.y = weights.y,
                       sd.x = sd.x, sd.y = sd.y, r.xy = r.xy,
                       approx.solution = approx.solution)

    # Define errors, the error correlation and weights
    mean.xi <- apply(x, 1, mean)
    mean.yi <- apply(y, 1, mean)
    x_errors<- x - mean.xi
    y_errors <- y - mean.yi
    for (i in 1:nrow(x_errors)) {
      r.xy[i] <- f_corr_row(x_errors[i, ],
                           y_errors[i, ])
      weights.x[i] <- 1 / f_var_row(x[i, ])
      weights.y[i] <- 1 / f_var_row(y[i, ])
    }
    x_original <- x
    y_original <- y
    x <- as.matrix(stack(data.frame(t(x_original)))[1])
    y <- as.matrix(stack(data.frame(t(y_original)))[1])
  }

  # initial value of the slope is the olse
  ols_reg <- f_ols_reg(x, y)
  slope <- ols_reg$slope


  if (approx.solution == F) {

    # algorithm to find york slope
    slope_diff <- 10
    count <- 0
    slope_per_iteration <- NULL
    alpha <- sqrt(weights.x * weights.y)
    while (slope_diff > tolerance) {
      slope_old <- slope
      Weight <- alpha^2 / (slope^2 * weights.y + weights.x -
                             2 * slope * r.xy * alpha)
      if (mult.samples == TRUE) {
        Weight <- rep(Weight, each = ncol(x_original))
      }
      Weight_sum <- sum(Weight)
      x_bar <- sum(Weight * x) / Weight_sum
      y_bar <- sum(Weight * y) / Weight_sum
      x_centered <- x - x_bar
      y_centered <- y - y_bar
      beta <- Weight * ((x_centered / weights.y) + (slope * y_centered /
                                                      weights.x) -
                          (slope * x_centered + y_centered) * r.xy / alpha)
      Q1 <- sum(Weight * beta * y_centered)
      Q2 <- sum(Weight * beta * x_centered)
      slope <- Q1 / Q2
      slope_diff <- abs(slope - slope_old)
      count <- count + 1
      slope_per_iteration <- append(slope_per_iteration, slope)
      if (count > max.iterations) {
        stop("\nThe slope coefficient does not converge after ",
             count, paste(" iterations. \nHint: You may reduce the tolerance level",
             "or increase the maximum number of iterations.", sep = " "),
             cat("Slope coefficient for last 5 iterations:"),
             for (i in 4:0){
               cat("\n\t", count - i, "\t", slope_per_iteration[count - i])},
             cat("\n"))
      }
    }
    slope_per_iteration <- data.frame("slope" = slope_per_iteration)

  } else {
    if (any(r.xy != 0)) {
      stop(paste("There is no approximate solution in case of correlation",
                 "between x and y errors!", sep = " "))
    }
    ## Apply formula and use lm estimate as intitial value for b
    alpha <- sqrt(weights.x * weights.y)
    Weight <- alpha^2 / (slope^2 * weights.y + weights.x)

    # see formula 19 and 20 for the following:
    Weight_sum <- sum(Weight)
    x_bar <- sum(Weight * x) / Weight_sum
    y_bar <- sum(Weight * y) / Weight_sum
    x_centered <- x - x_bar
    y_centered <- y - y_bar
    # calculate alpha_cubic, beta_cubic and gamma_cubic. See York 66 page 1084
    xy <- x_centered * y_centered
    xW_w <- x_centered^2 * Weight^2 / weights.x
    yW_w <- y_centered^2 * Weight^2 / weights.x
    alpha_cubic <- 2 * sum(xy * Weight^2 / weights.x) /
                    (3 * sum(xW_w))
    beta_cubic <- (sum(yW_w) - sum(Weight * x_centered^2)) /
                    (3 * sum(xW_w))
    gamma_cubic <- - sum(xy * Weight) / (sum(x_centered^2 *
                    Weight^2 / weights.x))

    # use formula of York 66 to find slope, given on page 1084
    phi <- acos((alpha_cubic^3 - 3 /2 * alpha_cubic * beta_cubic + 0.5 *
                        gamma_cubic) /
                  (alpha_cubic^2 - beta_cubic)^(3 / 2))

    sol_cubic2 <- alpha_cubic + 2 * (alpha_cubic^2 - beta_cubic)^0.5 *
                        cos( 1 / 3 *(phi + 2 * pi * c(0:2)))
    ols_range <- c(ols_reg$slope - 4 * ols_reg$se_slope,
                   ols_reg$slope + 4 * ols_reg$se_slope)
    pick_right_root <- which(sol_cubic2 >= ols_range[1] &
                               sol_cubic2 <= ols_range[2])
    slope <- sol_cubic2[pick_right_root]
    if (length(slope) == 0 | length(slope) > 1) {
      stop("An approximate solution does not exist!")
    }

    beta <- Weight * (x_centered / weights.y) + (slope * y_centered / weights.x)
    count <- 0
    slope_per_iteration <- data.frame("slope" = slope)
  }
  # York intercept
  intercept <- y_bar - slope * x_bar

  # SE of coefficients
  x_adj <- x_bar + beta
  x_mean <- sum(Weight * beta) / (Weight_sum * (length(x) - 2))
  u <- x_adj - x_mean
  sigma_slope <- sqrt(1 / sum(Weight * u^2))
  sigma_intercept <- sqrt(x_mean^2 * sigma_slope^2 + 1 / Weight_sum)

  # Goodness of fit + Test (H0: S <= df)
  S <- sum(Weight * (y - slope * x - intercept)^2)
  chisq_df <- (length(x) - 2)
  reduced_chisq <- S / chisq_df
  sigma_chisq <- sqrt(2 / chisq_df)
  if (mult.samples == F) {
    p_value <- 1 - pchisq(S, df = chisq_df)
    test_result <- if (p_value > 0.1 ) {
      "The assumption of a good fit cannot be rejected."
    } else if (p_value < 0.01) {
      paste("The assumption of a good fit can be rejected",
            "at a significance level of 1%.", sep = " ")
    } else if (p_value < 0.05) {
      paste("The assumption of a good fit can be rejected",
            "at a significance level of 5%.", sep = " ")
    } else {
      paste("The assumption of a good fit can be rejected",
            "at a significance level of 10%.", sep = " ")
    }
  } else {
    p_value <- "No p-value for multiple sample case."
    test_result <- "No test results for multiple sample case."
  }



  df_regression <- 2 * (length(x) - 1)

  # fitted values
  fitted_y <- intercept + slope * x

  # residuals
  c <- r.xy * alpha
  x_residuals <- (Weight * (intercept + slope * x - y)
                  * (c - slope * weights.y)) /
    (weights.y * weights.x)
  y_residuals <- (Weight * (intercept + slope * x - y) *
                    (weights.x - slope * c)) / (weights.y * weights.x)

  # define output
  york_reg <- matrix(c(intercept, slope, sigma_intercept, sigma_slope),
                     nrow = 2)
  rownames(york_reg) <- c("intercept", "slope")
  colnames(york_reg) <- c("Estimate", "Std_Error")

  weights_matrix <- matrix(c(weights.x, weights.y), ncol = 2)
  colnames(weights_matrix) <- c("weights of X_i", "weights of Y_i")
  if (mult.samples == F) {
    data <- matrix(c(x, y, sd.x, sd.y, r.xy), ncol = 5)
    colnames(data) <- c("x", "y", "sd.x", "sd.y", "r.xy")
  } else {
    data <- list("x" = x_original, "y" = y_original, "x_errors" = x_errors,
                 "y_errors" = y_errors, "r.xy" = r.xy, "mean.x.i" = x,
                 "mean.y.i" = y)
  }
  york_arguments <- list("tolerance" = tolerance, "max.iterations" = max.iterations,
                         "mult.samples" = mult.samples, "approx.solution" =
                           approx.solution)
  chisq_test_results <- list("test.result" = test_result, "p.value" = p_value)

  output <- list("coefficients" = york_reg,
                 "weights" = weights_matrix,
                 "x.residuals" = x_residuals,
                 "y.residuals"= y_residuals,
                 "fitted.y"=fitted_y,
                 "df.regression" = df_regression,
                 "weighted.mean.x" = x_bar,
                 "weighted.mean.y" = y_bar ,
                 "reduced.chisq" = reduced_chisq,
                 "std.error.chisq" = sigma_chisq,
                 "Overall.significance.of.fit" = chisq_test_results,
                 "number.of.iterations" = count,
                 "slope.after.each.iteration" = slope_per_iteration,
                 "ols_summary" = ols_reg[-c(1:2)],
                 "york.arguments" = york_arguments,
                 "data" = data)
  attr(output, "class") <- "york"

  return(output)
}

