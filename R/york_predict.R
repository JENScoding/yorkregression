
york.predict <- function(york.output, new) {
  if (class(york.output) != "york") {
    stop("Input must be of class york (Output of york function)")
  }
  york.output$coefficients.york[1, 1] + york.output$coefficients.york[2, 1] * new
}

# new <- c(2,3,4,3.6)
# york.predict(york.output = york.output, new = new)
#
#
# york.predict(york.output = york.output, new = 1.5)




