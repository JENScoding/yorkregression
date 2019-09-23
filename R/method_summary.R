#' @title
#'  Summarizing York Model Fits
#'
#' @description
#'  Summary method for class "york".
#'
#' @param object
#'  a fitted object of class inheriting from "york".
#' @param ...
#'  additional arguments affecting the summary produced.
#'
#' @return
#'  The function summary.york computes and returns a list, of summary
#'  statistics of the fitted york model given in object.
#'
#' \describe{
#'
#'   \item{x_residuals}{summary of the the weighted x residuals.}
#'   \item{y_residuals}{summary of the the weighted y residuals.}
#'   \item{coefficients}{a matrix which contains the York estimates for
#'     intercept and slope of the best-fit straight line with their respective
#'     standard errors.}
#'   \item{Regression_Test}{a vector of strings with the information of
#'     the value of the test statistic the degrees of freedom and the p-value
#'    as first element and the test result whether \eqn{H0} can be rejected
#'    or not as second element.}
#' }
#'
#' @examples
#'  # Example: York's regression with weight data taken from Pearson (1901):
#'  x <- c(0.0, 0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4)
#'  y <- c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5)
#'  weights_x <- c(1e+3, 1e+3, 5e+2, 8e+2, 2e+2, 8e+1, 6e+1, 2e+1, 1.8, 1)
#'  weights_y <- c(1, 1.8, 4, 8, 20, 20, 70, 70, 1e+2, 5e+2)
#'
#'  # fit york model and summarize
#'  york_fit <- york(x, y, weights_x, weights_y, r_xy_errors = 0)
#'  summary(york_fit)
#'
#' @name summary.york
#'
#' @export
#'
summary.york <- function(object, ...) {

  if (class(object) != "york") {
    stop("Input must be of class york (Output of york function)")
  }
  if (object$york_arguments$mult_samples == FALSE) {
    x_resid <- object$x_residuals
    y_resid <- object$y_residuals
  } else {# mult_samples = TRUE
    x_resid <- as.vector(object$x_residuals)
    y_resid <- as.vector(object$y_residuals)
  }

  # to be shown in console
  x_residuals_summary <- round(summary(x_resid), 3)
  y_residuals_summary <- round(summary(y_resid), 3)
  p_value <- round(object$goodness_of_fit$p_value, 5)
  df <- object$goodness_of_fit$chisq_df
  test_statistic <- round(object$goodness_of_fit$chisq_statistic, 3)

  phrase <- object$goodness_of_fit$test_result
  regression_test = c(paste("Chisq-statistic:", test_statistic,
                            "on", df, "degrees of freedom, ",
                            "p-value:", p_value),
                      phrase)

  # define output
  output <- list("x_Residuals" = x_residuals_summary,
                 "y_Residuals" = y_residuals_summary,
                 "Coefficients" = round(object$coefficients, 5),
                 "Regression_Test" = regression_test)

  return(output)
}

