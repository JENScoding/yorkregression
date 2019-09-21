#' @title
#'  Internal Functions for Plot Method
#'
#' @description
#'  Functions used in plot.york method. For functions in "Internal Functions
#'  for Plot Method" the output is used for proceeding calculations or
#'  displaying in the plot.york method. All the functions in "Internal
#'  Functions for Plot Method" begin with f_p_ and then the name of the
#'  function is given.
#'
#'  f_p_influential determines influential points in a york regression.
#'  An influential point is defined as a point that significantly changes
#'  the slope coefficient if omitted. Assymptotic normality of the slope
#'  coefficient is assumed. The test is based on a significance level
#'  of alpha = 0.01.
#'
#'  f_p_chisq_test is a function that calculates the critical values for
#'  alpha = c(0.1, 0.05, 0.01) and the corresponding probabilities in a
#'  \eqn{\chi^2}--test. Additionally, values for the rejection region
#'  for alpha = 0.01 are given.
#'
#'
#' @name method_plot_internal
#'
#' @keywords
#'  internal
#'

# function to determine influential observations:
f_p_influential <- function(x_data, y_data, mult_samples,
                        sd_x, sd_y, r_xy_errors, coef, coef_se,
                        tolerance, max_iterations, x_mult, y_mult) {

  # first for non multiple case then for multiple
  if (mult_samples == FALSE) {

    # run york regression and omit one observation
    slope_influential <- NULL
    for (i in 1:length(x_data)) {
      slope_influential[i] <- york(x_data[-i], y_data[-i],
                               sd_x = sd_x[-i],
                               sd_y = sd_y[-i],
                               r_xy_errors = r_xy_errors[-i],
                               tolerance = tolerance,
                               max_iterations = max_iterations)[[1]][2, 1]
    }

    # determine significance
    t_statistic <- (coef - slope_influential) /
      coef_se
    p_value <- 2 * (1 - pnorm(abs(t_statistic)))
    detect_influential <- which(p_value <= 0.01)

    # write all influential observations with x and y coordinate in a data frame
    # and define output
    influential <- data.frame("x" = x_data[detect_influential],
                          "y" = y_data[detect_influential],
                          "which_row" = factor(detect_influential))

    # define titles for plot
    if (length(detect_influential) != 0) {
      legend_title = "Influential observation at observation point:"
      plot_title = paste("York's best-fit straight line with encircled ",
                         "potential influential points")
      plot_subtitle = "t-test with significance level alpha = 0.01"
    } else {
      legend_title = "."
      plot_title = paste("York's best-fit straight line with encircled ",
                         "potential influential points")
      plot_subtitle = paste("No influential points were detected for ",
                            "t-test with significance level alpha = 0.01")
    }
  } else { # mult.samples = TRUE

    # run york regression and omit one row of observations
    slope_influential <- NULL
    for (i in 1:nrow(x_mult)) {
      slope_influential[i] <- suppressWarnings(
        york(x_mult[-i,], y_mult[-i,],
             tolerance = tolerance,
             max_iterations =max_iterations,
             mult_samples = mult_samples)[[1]][2, 1])
    }

    # determine significance
    t_statistic <- (coef - slope_influential) /
      coef_se
    p_value <- 2 * (1 - pnorm(abs(t_statistic)))
    detect_influential <- which(p_value <= 0.01)

    # define output
    x_mult_t <- data.frame(t(x_mult[detect_influential,]))
    y_mult_t <- data.frame(t(y_mult[detect_influential,]))

    if (length(detect_influential) == 0) {
      x_mult_st <- NULL
      y_mult_st <- NULL
    } else {
      x_mult_st <- stack(x_mult_t)[, 1]
      y_mult_st <- stack(y_mult_t)[, 1]
    }

    # write all influentials with x and y coordinate in a data frame
    influential <- data.frame("x" = x_mult_st,
                              "y" = y_mult_st,
                              "which_row" = factor(rep(detect_influential,
                                                       each =
                                                         ncol(x_mult))))


    if (length(detect_influential) != 0) {
      colnames(influential) <- c("x", "y", "which_row")
      legend_title = "Influential observation(s) in row(s):"
      plot_title = paste("York's best-fit straight line with encircled ",
                       "potential influential points")
      plot_subtitle = "t-test with significance level alpha = 0.01"
    } else {
      legend_title = "."
      plot_title = paste("York's best-fit straight line with encircled ",
                         "potential influential points")
      plot_subtitle = paste("No influential points were detected for ",
                            "t-test with significance level alpha = 0.01")
    }
  }

  output <- list("detect_influential" = influential,
                 "leg_title" = legend_title,
                 "pl_title" = plot_title,
                 "pl_subtitle" = plot_subtitle)

  return(output)
}


#' @title
#'  Internal Functions for Plot Method
#'
#' @name method_plot_internal
#'
#' @keywords
#'  internal
#'
# function for Chi^2 test with test results
f_p_chisq_test <- function(chisq_df, p_value, chisq_statistic) {

  # set up parameters for visualisation
  df <- chisq_df
  options(scipen = 999)
  p_value <- round(p_value, 5)
  math <- " = "
  if (p_value < 0.00001) {
    p_value <- 0.00001
    math <- " < "
  }

  # calculate corresponding x and y values of the p-values,
  # the critical value and rejection region for the test (alpha = 0.01)
  x_p_value <- qchisq(1 - p_value, df)
  critical_value <- qchisq(c(0.9, 0.95, 0.99,
                             seq(0.9915, 0.995, length.out = 15),
                             seq(0.99501, 0.999, length.out = 14)), df)
  y_critical <- dchisq(critical_value, df)
  y_p_value <- dchisq(x_p_value, df)
  reject_region <- data.frame("x" = critical_value[4:32],
                              "y" = rep(0, length(critical_value) - 3),
                              "yend" = y_critical[4:32])
  test_statistic <- round(chisq_statistic, 3)
  if (test_statistic >= 100){
    adj <- 1.09
  } else {
    adj <- 1.17
  }

  x_limit <- (2 * (1 + log(df^2) / df)) * df

  output <- list("df" = df,
                 "test_statistic" = test_statistic,
                 "p_value" = p_value,
                 "x_p_value" = x_p_value,
                 "y_p_value" = y_p_value,
                 "critical_value" = critical_value,
                 "y_critical" = y_critical,
                 "reject_region" = reject_region,
                 "adj" = adj,
                 "math" = math,
                 "x_limit" = x_limit)

  return(output)
}
