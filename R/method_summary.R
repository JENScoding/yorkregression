### york in Least Squares Fitting Of A Straight Line With Correlated Errors ###
## Prediction function

#' @title Summarizing York Model Fits
#'
#' @description Summary method for class "york".
#'
#' @param object A fitted object of class inheriting from "york".
#' @param ... additional arguments affecting the summary produced.
#'
#' @return The function summary.york computes and returns a list, of summary
#' statistics of the fitted york model given in object.
#'
#'   \describe{
#'
#'   \item{x_residuals}{the weighted x residuals.}
#'   \item{y_residuals}{the weighted y residuals.}
#'   \item{coefficients}{ll.}}
#'
#' @examples
#' # Example: York's regression with weight data taken from Pearson (1901):
#' x <- c(0.0, 0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4)
#' y <- c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5)
#' weights_x <- c(1e+3, 1e+3, 5e+2, 8e+2, 2e+2, 8e+1, 6e+1, 2e+1, 1.8, 1)
#' weights_y <- c(1, 1.8, 4, 8, 20, 20, 70, 70, 1e+2, 5e+2)
#' york_fit <- york(x = x, y = y, weights_x = weights_x, weights_y = weights_y,
#' r_xy_errors = 0)
#' summary(object = york_fit)
#'
#' @name summary.york
#' @export
summary.york <- function(object, ...) {

  if (class(object) != "york") {
    stop("Input must be of class york (Output of york function)")
  }
  x_residuals_summary <- round(summary(object$x_residuals), 3)
  y_residuals_summary <- round(summary(object$y_residuals), 3)
  p_value <- round(object$goodness_of_fit$p_value, 5)
  df <- object$goodness_of_fit$chisq_df
  test_statistic <- round(object$goodness_of_fit$chisq_statistic, 3)

  # define output
  output <- list("x_Residuals" = x_residuals_summary,
                 "y_Residuals" = y_residuals_summary,
                 "Coefficients" = round(object$coefficients, 5),
                 "Regression_Test" = paste("Chisq-statistic:",
                                           test_statistic, "on",
                                           df, "degrees of freedom, ",
                                           "p-value:", p_value))
  attr(output, "class") <- "summary.york"

  return(output)
}




