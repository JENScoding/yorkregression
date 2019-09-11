
york.predict <- function(york.output, new) {
  if (class(york.output) != "york") {
    stop("Input must be of class york (Output of york function)")
  }
  predict.y <- data.frame("x" = new, "predicted.y" = york.output$coefficients.york[1, 1] + york.output$coefficients.york[2, 1] * new)
  colnames(predict.y) <- c("x", "predicted.y")

  output <- list("prediction" = predict.y)
  attr(output, "class") <- "york"
  return(output)
}






