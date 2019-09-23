#' @title
#'  Predict method for York Model Fits
#'
#' @description
#'  Predictions from a model fitted by \code{york}.
#'
#' @param object
#'  a fitted object of class inheriting from "york".
#' @param newdata
#'  a vector or a dataframe to be used for the prediction.
#'  If not specified the original x-values with their respective
#'  fitted y-values will be returned.
#' @param ...
#'  additional arguments affecting the predictions produced.
#'
#' @return
#'  Returns a list, containing a dataframe with predicted values.
#'
#' \describe{
#'
#'   \item{prediction}{a dataframe containing either the original x-values
#'     or the specified newdata in the first column and the corresponding
#'     predicted values in the second column.}
#' }
#'
#' @examples
#'  # Example: York's regression with weight data taken from Pearson (1901):
#'  x <- c(0.0, 0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4)
#'  y <- c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5)
#'  weights_x <- c(1e+3, 1e+3, 5e+2, 8e+2, 2e+2, 8e+1, 6e+1, 2e+1, 1.8, 1)
#'  weights_y <- c(1, 1.8, 4, 8, 20, 20, 70, 70, 1e+2, 5e+2)
#'
#'  # fit york model
#'  york_fit <- york(x, y, weights_x, weights_y, r_xy_errors = 0)
#'
#'  # values for prediction
#'  new <- c(2,3,4,3.6)
#'  predict(york_fit, new)
#'
#' @name predict.york
#'
#' @export
#'
predict.york <- function(object, newdata = NULL, ...) {

  if (class(object) != "york") {
    stop("Input must be of class york (Output of york function)")
  }
  if (object$york_arguments$mult_samples == FALSE) {
    x_data <- object$data[, 1]
    y_data <- object$data[, 2]
  } else { # mult_sample = TRUE
    x_t <- data.frame(t(object$data$x))
    y_t <- data.frame(t(object$data$y))
    x_data <- stack(x_t)[, 1]
    y_data <- stack(y_t)[, 1]
  }

  # show fitted y's if now newdata is given
  if (is.null(newdata)) {
    predict_y <- data.frame("x" = x_data, "predict_y" =
                              object$fitted_y)
    colnames(predict_y) <- c("x", "predict_y")
  } else {

    # show fitted y for given newdata
    predict_y <- data.frame("x" = newdata, "predict_y" =
                              object$coefficients[1, 1] +
                              object$coefficients[2, 1] * newdata)
    colnames(predict_y) <- c("x", "predict_y")
  }

  output <- list("prediction" = predict_y)

  return(output)
}
