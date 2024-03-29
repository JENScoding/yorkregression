% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal_functions.R
\name{internal_functions_I}
\alias{internal_functions_I}
\alias{f_rewrite}
\alias{f_rewrite_mult}
\alias{f_var_row}
\alias{f_corr_row}
\alias{f_ols_reg}
\alias{f_cubic_root}
\alias{f_define_output}
\title{Internal Functions I}
\usage{
f_rewrite(x, y, weights_x = NULL, weights_y = NULL, sd_x = NULL,
  sd_y = NULL, r_xy_errors = NULL)

f_rewrite_mult(x, y)

f_var_row(x)

f_corr_row(x, y)

f_ols_reg(x, y)

f_cubic_root(x, y, weights_x, weights_y, r_xy, slope_ols, se_slope_ols)

f_define_output(intercept, slope, sigma_intercept, sigma_slope, weights_x,
  weights_y, mult_samples, x, y, sd_x, sd_y, r_xy_errors, x_data,
  y_original, x_errors, y_errors, mean_x_i, mean_y_i, slope_per_iteration,
  tolerance, max_iterations, approx_solution, S, chisq_df, p_value,
  test_result)
}
\description{
Function used in york function. For Internal Functions I
 the output is used for proceeding calculations or displaying in the
 york function. All the Internal Functions I functions begin with f_ and
 then the name of the function is given.

 f_rewrite rewrites the input data to make it suitable for the algorithm.
 If there are missing values in the input data, they will be omitted.
 If only one value of the weights or sd for x and y respectively, or
 only one value for the error correlation is given, this value is
 repeated as many times as the length of the x and y vector.

 f_var_row calculates the variance of the respective rows of a data frame
 and f_corr_row the correlations of the respective rows of a data frame.

 f_ols_reg returns the fit of a linear model. The OLS slope coefficient
 is then used in the first iteration of the york algorithm.

 f_cubic_root avoids the iterative determination of the slope coefficient.
 It solves the cubic equation for the slope coefficient.
 It is just an approximation, for more accurate results the slope has to be
 determined iteratively. Only possible when the error correlation
 is equal to 0. For further detail see references.

 f_define_output returns a list that is used for the york function
 output. It mainly restructures the information given in the york
 function.
}
\references{
York, Derek. "Least-squares fitting of a straight line.", Canadian Journal
 of Physics 44.5 (1966), pp. 1079-1086.
}
\keyword{internal}
