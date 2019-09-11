stop.mult.sample <- function(exact.solution, x, y){
  if (exact.solution == T) {
    stop("There is no exact solution in case of multiple samples!")
  }
  if (is.data.frame(x) == F|| is.data.frame(y) == F) {
    stop("Inputs x and y must be of class data.frame!")
  }
  if (ncol(x) != ncol(y) || nrow(x) != nrow(y)) {
    stop("x and y must have the same number of columns/ rows")
  }
  if (ncol(x) == 1 || ncol(y) == 1) {
    stop("You need more than one sample of x and y, if you specify mult.samples = T")
  }
}
