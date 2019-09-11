### york in Least Squares Fitting Of A Straight Line With Correlated Errors ###
## Prediction function

#' @title Predict Method for York Fits
#'
#' @description Obtains predictions from a fitted York's regression.
#'
#' @param x A fitted object of class inheriting from "york".
#' @param newdata A vector or a dataframe in which to look for variables with which
#'   to predict.
#'
#' @return Returns a list, containing a dataframe with predicted values as first element.
#'
#'   \describe{
#'
#'   \item{prediction}{A dataframe containing the \code{newdata} in the first
#'   column and the corresponding predicted values in the second column.}
#'
#' @examples
#' # Example: York's regression with weight data taken from Pearson (1901):
#' x <- c(0.0, 0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4)
#' y <- c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5)
#' weights.x <- c(1e+3, 1e+3, 5e+2, 8e+2, 2e+2, 8e+1, 6e+1, 2e+1, 1.8, 1)
#' weights.y <- c(1, 1.8, 4, 8, 20, 20, 70, 70, 1e+2, 5e+2)
#' r.xy <- 0
#' fit <- york(x = x, y = y, weights.x = weights.x, weights.y = weights.y)
#' # predict with:
#' data<- c(2,3,4,3.6)
#' york.predict(york.output = fit, newdata = data)
#'
#' @name york
#' @export
york.predict <- function(york.output, newdata) {
  if (class(york.output) != "york") {
    stop("Input must be of class york (Output of york function)")
  }
  predict.y <- data.frame("x" = newdata, "predicted.y" = york.output$coefficients.york[1, 1] + york.output$coefficients.york[2, 1] * newdata)
  colnames(predict.y) <- c("x", "predicted.y")

  output <- list("prediction" = predict.y)
  return(output)
}






