
error_wrong_input1 <- function(x, y, weights.x = NULL, weights.y = NULL,
                               sd.x = NULL, sd.y = NULL, r.xy = NULL) {

  # Specify all wrong inputs if mult.sample == FALSE
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
}
