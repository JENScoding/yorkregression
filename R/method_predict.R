### york in Least Squares Fitting Of A Straight Line With Correlated Errors ###
## Prediction function

#' @title Predict method for York fits
#'
#' @description Obtains predictions from a fitted York's regression.
#'
#' @param york_output A fitted object of class inheriting from "york".
#' @param newdata A vector or a dataframe in which to look for variables
#' with which to predict.
#'
#' @return Returns a list, containing a dataframe with predicted values
#' as first element.
#'
#'   \describe{
#'
#'   \item{prediction}{A dataframe containing the \code{newdata} in the first
#'   column and the corresponding predicted values in the second column.}}
#'
#' @examples
#' # Example: York's regression with weight data taken from Pearson (1901):
#' x <- c(0.0, 0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4)
#' y <- c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5)
#' weights_x <- c(1e+3, 1e+3, 5e+2, 8e+2, 2e+2, 8e+1, 6e+1, 2e+1, 1.8, 1)
#' weights_y <- c(1, 1.8, 4, 8, 20, 20, 70, 70, 1e+2, 5e+2)
#' york_output <- york(x = x, y = y, weights_x = weights_x, weights_y = weights_y,
#' r_xy_errors = 0)
#' data <- c(2,3,4,3.6)
#' york_predict(york_output = york_output, newdata = data)
#'
#' @name predict.york
#' @export
predict.york <- function(york_output, newdata = NULL) {

  if (class(york_output) != "york") {
    stop("Input must be of class york (Output of york function)")
  }
  if (york_output$york_arguments$mult_samples == FALSE) {
    x_data <- york_output$data[, 1]
    y_data <- york_output$data[, 2]
  } else {
    x_data <- stack(york_output$data$x)[, 1]
    y_data <- stack(york_output$data$y)[, 1]
  }
  if (is.null(newdata)) {
    predict_y <- data.frame("x" = x_data, "predicted_y" =
                              york_output$fitted_y)
  } else {
    predict_y <- data.frame("x" = newdata, "predicted_y" =
                              york_output$coefficients[1, 1] +
                              york_output$coefficients[2, 1] * newdata)
    colnames(predict_y) <- c("x", "predicted_y")
  }

  output <- list("prediction" = predict_y)

  return(output)
}






