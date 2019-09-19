#' @title
#'  Fitting Linear Models With York's Method.
#'
#' @description
#'  The york function is used to fit a model when both x and y variables are subject
#'  to measurement errors. The york function returnes an object of class "york", which is a fit of a
#'  York's regression. The model can also take into account correlations between
#'  x and y errors.
#'
#' @param x
#'  a 1 times n vector or in case of multiple samples a dataframe, where
#'  the first column represents the first sample, of the \code{x}-variable(s).
#' @param y
#'  a 1 times n vector or in case of multiple samples a a dataframe, where
#'  the first column represents the first sample, of the \code{y}-variable(s).
#' @param weights_x
#'  the prespecified 1 times n weighting vector for \code{x}-values.
#' @param weights_y
#'  the prespecified 1 times n weighting vector for \code{y}-values.
#' @param r_xy_errors
#'  the prespecified correlation coefficient between the errors
#'  in \code{x} and \code{y}. Either a 1 times n vector or a single value.
#' @param tolerance
#'  the tolerance for convergence of the slope coefficent. The default
#'  is \code{1e-5}.
#' @param max_iterations
#'  the maximum number of iterations for convergence. The default is \code{50}.
#' @param sd_x
#'  the standard error of the \code{x}-values. If the true errors
#'  in x are known.
#' @param sd_y
#'  the standard error of the \code{y}-values. If the true errors
#'  in y are known.
#' @param mult_samples
#'  \code{logical}. Default is \code{FALSE}. Change to TRUE, if the input is
#'  a data frame with multiple samples of both x and y.
#' @param approx_solution
#'  \code{logical}. Default is \code{FALSE}. Change to TRUE, if you want an
#'  approximate solution of the slope coefficients. No iteration is needed.
#'  TRUE is not recommended when the accuracy for the slope coefficient shall be
#'  more than one decimal point.
#'
#' @details
#'  The york function implements the algorithm for the problem of the
#'  best-ﬁt straight line to independent points with errors in both x and y
#'  variables. General York (1969) solution according to the algorithm of Wehr & Saleska
#'  (2017).
#'
#'  Given \eqn{n} pairs of \eqn{(x_i, y_i), i = 1, \ldots, n}, their
#'  weights \eqn{(\omega(x_i), \omega(y_i)), i = 1, \ldots, n} or their standard
#'  errors \eqn{sd(x_i)} and \eqn{sd(y_i)}, the \code{york} function finds the
#'  best-fit straight line using the algorithm of York et al. (1966)/ York et
#'  al. (1969) as presented in Wehr & Saleska (2017). In addition, the function
#'  provides numerous statistics, parameters and goodness of fit criteria. If
#'  the data contains NA values they will be omitted.
#'
#' @return
#'  York Returns an object of class "york". An
#'  object of class "york" is a list containing the following components:
#'
#'  \describe{
#'
#'  \item{coefficients}{a matrix which contains the York estimates for intercept
#'    and slope of the best-fit straight line with their respective standard errors.}
#'  \item{x_residuals}{a vector of the York x-residuals.}
#'  \item{y_residuals}{a vector of the York y-residuals.}
#'  \item{fitted_y}{a vector of the fitted York y-values.}
#'  \item{weights}{a matrix representation of the prespecified or calculated
#'    weights for the x- and y-observations.}
#'  \item{data}{a data matrix which contains as columns the observed points x-,
#'    y-values, the errors sd_x- and sd_y and the correlation of the errors
#'    (error_correlation). If the input are multiple samples the data element will be
#'    a list containing the observed points x-, y-values, the errors sd_x- and sd_y
#'    and the correlation of the errors (error_correlation), the errors in x and y
#'    (x_errors// y_errors) and the mean of each obersvation i for variable x
#'    and y, respectively (mean_x_i// mean).}
#'  \item{reduced_chisq}{the reduced chi-squared statistic
#'    (See \url{https://en.wikipedia.org/wiki/Reduced_chi-squared_statistic}), i.e.
#'    the goodness of fit measure of York's regression.}
#'  \item{se_chisq}{the standard error of the chi-squared statistic.}
#'  \item{goodness_of_fit}{a list with the test results of a \eqn{\chi^2}--test,
#'    containing the test-statistic, the degrees of freedom, the p-value and
#'    a string saying whether \eqn{H0} (the assumption of a good fit) can be
#'    rejected or not for \eqn{\alpha = 0.01}.}
#'  \item{n_iterations}{the total number of iterations.}
#'  \item{slope_per_iteration}{the York slope after each iteration.}
#'  \item{weighted_mean_x}{the weighted.mean of x.}
#'  \item{weighted_mean_y}{the weighted.mean of y.}
#'  \item{ols_summary}{a list containing ols statistics.}
#'  \item{york_arguments}{a list containing the specfied input arguments.}
#'  }
#'
#' @references
#'  Wehr, Richard, and Scott R. Saleska. "The long-solved problem of
#'  the best-fit straight line: Application to isotopic mixing lines."
#'  Biogeosciences 14.1 (2017). pp. 17-29.
#'
#'  York, Derek. "Least squares fitting of a straight line with correlated
#'  errors." Earth and planetary science letters 5 (1968), pp. 320-324.
#'
#'  York, Derek. "Least-squares fitting of a straight line.", Canadian Journal
#'  of Physics 44.5 (1966), pp. 1079-1086.
#'
#'  York, Derek, et al. "Unified equations for the slope, intercept, and
#'  standard errors of the best straight line." American Journal of Physics 72.3
#'  (2004), pp. 367-375.
#'
#' @examples
#'  # Example: York's regression with weight data taken from Pearson (1901):
#'  x <- c(0.0, 0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4)
#'  y <- c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5)
#'  weights_x <- c(1e+3, 1e+3, 5e+2, 8e+2, 2e+2, 8e+1, 6e+1, 2e+1, 1.8, 1)
#'  weights_y <- c(1, 1.8, 4, 8, 20, 20, 70, 70, 1e+2, 5e+2)
#'  r_xy_errors <- 0
#'  york(x, y, weights_x, weights_y, r_xy_errors)
#'
#'  # Example: York's regression arbitrary values for sd_x and sd_y:
#'  x <- c(0.0, 0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4)
#'  y <- c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5)
#'  sd_x <- 0.2
#'  sd_y <- 0.4
#'  r_xy_errors <- 0.3
#'
#'  # fit york model
#'  york(x, y, sd_x = sd_x, sd_y = sd_y, r_xy_errors = r_xy_errors)
#'
#' \dontrun{
#'  # Example: No standard errors or weights specified
#'  york(x, y, r_xy_errors = 0)
#'
#'  # Example: You can't specify weights and standard errors at the same time
#'  york(x , y, sd_x, sd_y, weights_x, weights_y, r_xy_errors = 0)
#' }
#'
#' @name york
#'
#' @importFrom stats pchisq
#' @importFrom utils stack
#'
#' @export
#'
york <- function(x, y, weights_x = NULL, weights_y = NULL, r_xy_errors = NULL,
                 tolerance = 1e-5, max_iterations = 50, sd_x = NULL, sd_y = NULL,
                 mult_samples = FALSE, approx_solution = FALSE) {

  if (mult_samples == FALSE) {

    # rewrite input and delete rows with NA values
    re_input <- f_rewrite(x, y, weights_x, weights_y,
                       sd_x, sd_y, r_xy_errors)
    x <- re_input$x
    y <- re_input$y
    weights_x <- re_input$weights_x
    weights_y <- re_input$weights_y
    sd_x <- re_input$sd_x
    sd_y <- re_input$sd_y
    r_xy_errors <- re_input$r_xy_errors
    x_data <- NULL
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

    # rewrite input and delete rows with NA values
    re_input <- f_rewrite_mult(x, y)
    x <- re_input$x
    y <- re_input$y

    # Define errors, error correlation and weights in multiple sample case
    mean_x_i <- apply(x, 1, mean)
    mean_y_i <- apply(y, 1, mean)
    x_errors<- x - mean_x_i
    y_errors <- y - mean_y_i
    r_xy_errors <- f_corr_row(x_errors,
                           y_errors)
    sd_x <- f_var_row(x)
    sd_y <- f_var_row(y)
    weights_x <- 1 / sd_x
    weights_y <- 1 / sd_y

    x_data <- x
    y_original <- y
    x <- as.matrix(stack(data.frame(t(x_data)))[1])
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
        Weight <- rep(Weight, each = ncol(x_data))
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
                                sd_y, r_xy_errors, x_data, y_original,
                                x_errors, y_errors, mean_x_i, mean_y_i,
                                slope_per_iteration, tolerance, max_iterations,
                                approx_solution, S, chisq_df, p_value,
                                test_result)

  output <- list("coefficients" = def_output$york_reg,
                 "x_residuals" = x_residuals,
                 "y_residuals" = y_residuals,
                 "fitted_y" = fitted_y,
                 "weights" = def_output$weights_matrix,
                 "data" = def_output$data,
                 "reduced_chisq" = reduced_chisq,
                 "se_chisq" = sigma_chisq,
                 "goodness_of_fit" = def_output$chisq_test_results,
                 "n_iterations" = count,
                 "slope_per_iteration" = def_output$slope_per_iteration,
                 "weighted_mean_x" = x_bar,
                 "weighted_mean_y" = y_bar,
                 "ols_summary" = ols_reg[-c(1:2)],
                 "york_arguments" = def_output$york_arguments)
  attr(output, "class") <- "york"

  return(output)
}

