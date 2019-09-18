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
#'@param max_iterations The maximum number of iterations for convergence.
#'   Default is \code{50}.
#' @param weights_x The prespecified 1 times n weights vector for
#'   \code{X}-values
#' @param weights_y The prespecified 1 times n weights vector for
#'   \code{Y}-values
#' @param sd_x The standard error of the \code{X}-values
#' @param sd_y The standard error of the \code{Y}-values
#' @param r_xy_errors The prespecified correlation coefficient between the errors in
#'   \code{X} and \code{Y}
#' @param mult_samples An indicator if the multiple samples option is turned on
#' or not. The standard value is \code{FALSE}
#' @param approx_solution An indicator if the approximate solution option is turned on
#' or not. The standard value is \code{FALSE}
#' @return York Returns an object of class "york" the York regression for the
#'   \code{x} and \code{y} data for either specified weights \code{weights_x}
#'   and \code{weights_y} or specified standard errors \code{sd_x} and
#'   \code{sd_y} An object of class "york" containing the following components:
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
#'   \item{x_residuals}{a vector of the
#'   York X-residuals}
#'   \item{y_residuals}{a vector of the York Y-residuals}
#'   \item{fitted_y}{a vector of the fitted York Y-values}
#'   \item{df_regression}{the number of degrees of freedom
#'   (See \url{https://en.wikipedia.org/wiki/Degrees_of_freedom_(statistics)})
#'   of York's regression}
#'   \item{weighted_mean_x}{the weighted.mean of X}
#'   \item{weighted_mean_y}{the weighted.mean of Y}
#'   \item{reduced_chisq}{the reduced chi-squared statistic
#'   (See \url{https://en.wikipedia.org/wiki/Reduced_chi-squared_statistic}),
#'   i.e. the goodness of fit measure of York's regression}
#'   \item{se_chisq}{the standard error of the chi-squared statistic}
#'   \item{n_iterations}{the total number of iterations}
#'   \item{slope_per_iteration}{the York slope after each iteration}
#'   \item{fitted_y.ols}{the fitted values for OLS}
#'   \item{residuals.ols}{the OLS residuals}
#'   \item{residual.sum.of.squares}{the residual sum of squares (RSS) of OLS}
#'   \item{total.sum.of.squares}{the total sum of squares of OLS}
#'   \item{se.of.reg.ols}{the standard error of the regression (SER) for OLS}
#'   \item{r.squared.ols}{the R squared of the OLS regression}
#'   \item{r.squared.adjusted.ols}{the adjusted R squared of OLS}
#'   \item{f.statistic.ols}{the F statistic of OLS}
#'   \item{fitted_y.orthogonal}{the fitted values for orthogonal regression (See
#'   \url{https://en.wikipedia.org/wiki/Deming_regression})}
#'   \item{data}{a data matrix which contains as columns the observed points
#'   X-, Y-, sd_x- and sd_y-values}}
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
#' weights_x <- c(1e+3, 1e+3, 5e+2, 8e+2, 2e+2, 8e+1, 6e+1, 2e+1, 1.8, 1)
#' weights_y <- c(1, 1.8, 4, 8, 20, 20, 70, 70, 1e+2, 5e+2)
#' r_xy_errors <- 0
#' york(x = x, y = y, tolerance = 1e-10, weights_x = weights_x,
#' weights_y = weights_y, r_xy_errors = r_xy_errors)
#'
#' # Example: York's regression arbitrary values for sd_x and sd_y:
#' x <- c(0.0, 0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4)
#' y <- c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5)
#' sd_x <- 0.2
#' sd_y <- 0.4
#' r_xy_errors <- 0.3
#' york(x = x, y = y, tolerance = 1e-10, sd_x = sd_x,
#' sd_y = sd_y, r_xy_errors = r_xy_errors)
#'
#' \dontrun{
#' # Example: No standard errors or weights specified
#' york(x, y, r_xy_errors = 0)
#' # Example: You can't specify weights and standard errors at the same time
#' york(x , y, sd_x, sd_y, weights_x, weights_y, r_xy_errors = 0)
#' # Example: x and y must have same length
#' york(x = c(0.0, 0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5),
#' y = c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5),
#' weights_x = weights_x, weights_y = weights_y,  r_xy_errors = 0)
#' }
#' @name york
#' @export
#' @importFrom stats pchisq
#' @importFrom utils stack
york <- function(x, y, weights_x = NULL, weights_y = NULL, tolerance = 1e-5,
                 max_iterations = 50, sd_x = NULL, sd_y = NULL, r_xy_errors = NULL,
                 mult_samples = FALSE, approx_solution = FALSE) {

  if (mult_samples == FALSE) {

    # rewrite input and delete rows with NA values
    input <- f_rewrite(x, y, weights_x, weights_y,
                       sd_x, sd_y, r_xy_errors)
    x <- input$x
    y <- input$y
    weights_x <- input$weights_x
    weights_y <- input$weights_y
    sd_x <- input$sd_x
    sd_y <- input$sd_y
    r_xy_errors <- input$r_xy_errors
    x_original <- NULL
    y_original <- NULL
    x_errors <- NULL
    y_errors <- NULL
    mean_x_i <- NULL
    mean_y_i <- NULL


    # expected errors for wrongly specified input
    exp_error_simple(x, y, weights_x, weights_y,
                       sd_x, sd_y, r_xy_errors)

  } else { # mult_samples = TRUE

    # expected errors for wrongly specified multiple sample input
    exp_error_multiple(x, y, weights_x, weights_y,
                       sd_x, sd_y, r_xy_errors, approx_solution)

    # Define errors, error correlation and weights in multiple sample case
    mean_x_i <- apply(x, 1, mean)
    mean_y_i <- apply(y, 1, mean)
    x_errors<- x - mean_x_i
    y_errors <- y - mean_y_i
    for (i in 1:nrow(x_errors)) {
      r_xy_errors[i] <- f_corr_row(x_errors[i, ],
                           y_errors[i, ])
      sd_x[i] <- f_var_row(x[i, ])
      sd_y[i] <- f_var_row(y[i, ])
      weights_x[i] <- 1 / sd_x[i]
      weights_y[i] <- 1 / sd_y[i]
    }
    x_original <- x
    y_original <- y
    x <- as.matrix(stack(data.frame(t(x_original)))[1])
    y <- as.matrix(stack(data.frame(t(y_original)))[1])
  }

  # initial value of the slope is olse slope
  ols_reg <- f_ols_reg(x, y)
  slope <- ols_reg$slope


  if (approx_solution == FALSE) {

    # algorithm to find york slope
    slope_diff <- 10
    count <- 0
    slope_per_iteration <- NULL
    alpha <- sqrt(weights_x * weights_y)
    while (slope_diff > tolerance) {
      slope_old <- slope
      Weight <- alpha^2 / (slope^2 * weights_y + weights_x -
                             2 * slope * r_xy_errors * alpha)
      if (mult_samples == TRUE) {
        Weight <- rep(Weight, each = ncol(x_original))
      }
      Weight_sum <- sum(Weight)
      x_bar <- sum(Weight * x) / Weight_sum
      y_bar <- sum(Weight * y) / Weight_sum
      x_centered <- x - x_bar
      y_centered <- y - y_bar
      beta <- Weight * ((x_centered / weights_y) + (slope * y_centered /
                                                      weights_x) -
                          (slope * x_centered + y_centered) * r_xy_errors / alpha)
      Q1 <- sum(Weight * beta * y_centered)
      Q2 <- sum(Weight * beta * x_centered)

      ### the slope
      slope <- Q1 / Q2

      # save values and start if slope_diff < tolerance next iteration
      slope_diff <- abs(slope - slope_old)
      count <- count + 1
      slope_per_iteration <- append(slope_per_iteration, slope)

      # error if it does not converge
      exp_error_convergence(count, max_iterations, slope_per_iteration)
    }

  } else { # approx_solution = TRUE

    # solve cubic problem and use estimates as approximation
    approx <- f_cubic_root(x, y, weights_x, weights_y,
                                           r_xy_errors, slope, ols_reg$se_slope)
    slope <- approx$slope
    Weight <- approx$Weight
    Weight_sum <- approx$Weight_sum
    x_bar <- approx$x_bar
    y_bar <- approx$y_bar
    alpha <- approx$alpha
    beta <- approx$beta
    count <- 0
    slope_per_iteration <- slope
  }

  # York intercept
  intercept <- y_bar - slope * x_bar

  # SE of coefficients
  x_adj <- x_bar + beta
  x_mean <- sum(Weight * beta) / (Weight_sum * (length(x) - 2))
  u <- x_adj - x_mean
  sigma_slope <- sqrt(1 / sum(Weight * u^2))
  sigma_intercept <- sqrt(x_mean^2 * sigma_slope^2 + 1 / Weight_sum)

  # Goodness of fit + Test (H0: S <= chisq_df)
  S <- sum(Weight * (y - slope * x - intercept)^2)
  chisq_df <- (length(x) - 2)
  reduced_chisq <- S / chisq_df
  sigma_chisq <- sqrt(2 / chisq_df)
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

  # fitted values
  fitted_y <- intercept + slope * x

  # residuals
  c <- r_xy_errors * alpha
  x_residuals <- (Weight * (intercept + slope * x - y)
                  * (c - slope * weights_y)) /
    (weights_y * weights_x)
  y_residuals <- (Weight * (intercept + slope * x - y) *
                    (weights_x - slope * c)) / (weights_y * weights_x)

  # define output
  def_output <- f_define_output(intercept, slope, sigma_intercept, sigma_slope,
                                weights_x, weights_y, mult_samples, x, y, sd_x,
                                sd_y, r_xy_errors, x_original, y_original,
                                x_errors, y_errors, mean_x_i, mean_y_i,
                                slope_per_iteration, tolerance, max_iterations,
                                approx_solution, S, chisq_df, p_value,
                                test_result)

  output <- list("coefficients" = def_output$york_reg,
                 "weights" = def_output$weights_matrix,
                 "x_residuals" = x_residuals,
                 "y_residuals" = y_residuals,
                 "fitted_y" = fitted_y,
                 "weighted_mean_x" = x_bar,
                 "weighted_mean_y" = y_bar ,
                 "reduced_chisq" = reduced_chisq,
                 "se_chisq" = sigma_chisq,
                 "goodness_of_fit" = def_output$chisq_test_results,
                 "n_iterations" = count,
                 "slope_per_iteration" = def_output$slope_per_iteration,
                 "ols_summary" = ols_reg[-c(1:2)],
                 "york_arguments" = def_output$york_arguments,
                 "data" = def_output$data)
  attr(output, "class") <- "york"

  return(output)
}

