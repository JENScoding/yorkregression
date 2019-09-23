#' @title
#'  Internal Functions II
#'
#' @description
#'  Functions used in york function. Internal Functions II
#'  compromise all possible error messages.
#'
#'  The exp_error_simple function returns all possible error messages
#'  if the input is a simple regression with only one sample of x and y.
#'
#'  The exp_error_multiple function returns all possible error messages if
#'  mult.sample = TRUE is specified
#'
#'  The exp_error_convergence function returns an error message if the slope
#'  coefficient does not converge. It will also print the value if the last
#'  determined slope coefficients.
#'
#' @name internal_functions_II
#'
#' @keywords
#'  internal
#'
# function to print error message if input is wrong for mult.samples = FALSE
exp_error_simple <- function(x, y, weights_x = NULL, weights_y = NULL,
                               sd_x = NULL, sd_y = NULL, r_xy_errors = NULL) {

  # Specify all wrong inputs if mult.sample = FALSE
  if (is.null(c(sd_x, sd_y, weights_x, weights_y))) {
    stop("Specify either standard errors or weights")
  }
  if(length(x) != length(y)) {
    stop("x and y must have the same length!")
  }
  if (length(r_xy_errors) != length(x)) {
    stop("Length of correlation vector must equal length of x")
  }
  if (any(r_xy_errors <= -1 | r_xy_errors >= 1)) {
    stop("Wrong input for r_xy_errors:
       r_xy_errors must be element of (-1, ... , 1)")
  }
  if (length(weights_x) != length(x) | length(weights_y) != length(y)) {
    stop("weights_x and weights_y must have the same length as x and y resp.!")
  }
  if (length(sd_x) != length(x) | length(sd_y) != length(y)) {
    stop("sd_x and sd_y must have the same length as x and y resp.!")
  }
}


#' @title Internal Functions II
#'
#' @name internal_functions_II
#'
#' @keywords internal
#'
# error function for multiple sample case
exp_error_multiple <- function(x, y, weights_x = NULL, weights_y = NULL,
                               sd_x = NULL, sd_y = NULL, r_xy_errors = NULL,
                               approx_solution = FALSE) {

  # Specify all wrong inputs if mult.sample = FALSE
  if (is.null(c(sd_x, sd_y, weights_x, weights_y)) == FALSE) {
    stop(paste("Standard errors and weights cannot be specified",
               "in multiple sample case!", sep = " "))
  }
  if (approx_solution == TRUE) {
    stop("There is no approximate solution in case of multiple samples!")
  }
  if (is.data.frame(x) == FALSE || is.data.frame(y) == FALSE) {
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


#' @title Internal Functions II
#'
#' @name internal_functions_II
#'
#' @keywords internal
#'
# error function for no convergence
exp_error_convergence <- function(count, max_iterations, slope_per_iteration) {

  if (count <= 3) {
    last <- count - 1
  } else {
    last <- 4
  }
  if (count == max_iterations) {
    stop("\nThe slope coefficient does not converge after ",
         count, paste(" iterations. \nHint: You may reduce the tolerance level",
                      "or increase the maximum number of iterations.", sep = " "),
         cat("Slope coefficient for the last", last + 1, "iterations:"),
         for (i in last:0){
           cat("\n\t", count - i, "\t", slope_per_iteration[count - i])},
         cat("\n"))
  }
}
