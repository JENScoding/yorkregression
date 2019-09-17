### expected errors with only one sample each (simple) (mult.sample = FALSE)
exp_error_simple <- function(x, y, weights.x = NULL, weights.y = NULL,
                               sd.x = NULL, sd.y = NULL, r.xy = NULL) {

  # Specify all wrong inputs if mult.sample = FALSE
  if (is.null(c(sd.x, sd.y, weights.x, weights.y))) {
    stop("Specify either standard errors or weights")
  }
  if(length(x) != length(y)) {
    stop("x and y must have the same length!")
  }
  if (length(r.xy) != length(x)) {
    stop("Length of correlation vector must equal length of x")
  }
  if (any(r.xy <= -1 | r.xy >= 1)) {
    stop("Wrong input for r.xy:
       r.xy must be element of (-1, ... , 1)")
  }
  if (length(weights.x) != length(x) | length(weights.y) != length(y)) {
    stop("weights.x and weights.y must have the same length as x and y resp.!")
  }
  if (length(sd.x) != length(x) | length(sd.y) != length(y)) {
    stop("sd.x and sd.y must have the same length as x and y resp.!")
  }
}

### expected errors when input has multiple samples (mult.sample = TRUE)
exp_error_multiple <- function(x, y, weights.x = NULL, weights.y = NULL,
                               sd.x = NULL, sd.y = NULL, r.xy = NULL,
                               approx.solution = FALSE) {

  # Specify all wrong inputs if mult.sample = FALSE
  if (approx.solution == T) {
    stop("There is no approximate solution in case of multiple samples!")
  }
  if (is.data.frame(x) == F|| is.data.frame(y) == F) {
    stop("Inputs x and y must be of class data.frame!")
  }
  if (ncol(x) != ncol(y) || nrow(x) != nrow(y)) {
    stop("x and y must have the same number of columns/ rows")
  }

  if (ncol(x) == 1 || ncol(y) == 1) {
    stop(paste("You need more than one sample of x and y",
               "if you specify mult.samples = T.", sep = " "))
  }
  if (ncol(x) < 5) {
    stop(paste("You need more than at least 4 samples",
               "of the x and y variables.", sep = " "))
  }
  if (ncol(x) < 10) {
    warning(paste("You have less than 10 samples",
                  "of the x and y variables.",
                  "Increasing the number of samples is recommended",
                  "in order to get accurate estimates.", sep = " "))
  }
}
