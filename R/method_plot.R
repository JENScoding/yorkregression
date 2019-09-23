#' @title
#'  Plot an york object
#'
#' @description
#'  Plots for york model fit are obtained. The plot function returns
#'  6 different plots.
#'  A plot for the fitted line, diagnostic plots and a trace plot
#'  (if the slope coefficient was determined approximately,
#'  no trace plot will be shown).
#'
#' @param x
#'  an object of class "york"
#' @param ...
#'  arguments to be passed to methods.
#'
#' @details
#'  Given an object of class "york" this function gives several york plots.
#'  The first plot shows York's best-fit straight line only. The
#'  second plot shows York's best fit straight line compared to OLS and
#'  orthogonal regression. Additionally, the different "center of gravity" are
#'  shown. The third plot is a "y-residuals vs. x-residuals plot". The fourth is
#'   a "x-residuals vs. fitted-y plot". The fith plot is a "y-residuals vs.
#'  fitted y" plot and the last plot is a "trace plot" which shows the slope
#'  after each interation, i.e. how the slope coefficient converges.
#'
#' @examples
#'  # Example: York's regression with weight data taken from Pearson (1901):
#'  x <- c(0.0, 0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4)
#'  y <- c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5)
#'  weights_x <- c(1e+3, 1e+3, 5e+2, 8e+2, 2e+2, 8e+1, 6e+1, 2e+1, 1.8, 1)
#'  weights_y <- c(1, 1.8, 4, 8, 20, 20, 70, 70, 1e+2, 5e+2)
#'  r_xy_errors <- 0
#'
#'  # fit york model
#'  york_fit <- york(x, y, weights_x, weights_y, r_xy_errors = 0)
#'
#'  # Obtain plots
#'  plot(york_fit)
#'
#' @name plot.york
#'
#' @importFrom
#'  ggplot2 ggplot aes geom_abline geom_point labs theme element_text
#'  draw_key_rect scale_colour_manual geom_vline geom_hline geom_smooth
#'  geom_line stat_function geom_segment annotate
#' @importFrom
#'  utils stack
#' @importFrom
#'  stats pnorm qchisq dchisq
#'
#' @method
#'  plot york
#' @S3method
#'  plot york
#'

plot.york <- function(x, ...) {
  if (class(x) != "york") {
    stop("Input must be of class york (Output of york function)")
  }

  ### Rename x and make variable y global
  york_fit <- x
  y <- NULL


  ### York's best fit line with data points
  if (york_fit$york_arguments$mult_samples == FALSE) {

    x_data <- york_fit$data[, 1]
    y_data <- york_fit$data[, 2]

    ddf <- data.frame(x = x_data, y = y_data)

    # plot fitted line and data
    plot_1 <- ggplot(data = ddf, aes(x = x,
                                   y = y)) +
      geom_point() +
      geom_abline(aes(slope = york_fit$coefficients[2, 1],
                      intercept = york_fit$coefficients[1, 1]), col = "red") +
      labs(title = "York's best-fit straight line",
           x = "x data", y = "y data") +
      theme(plot.title = element_text(hjust = 0.5))

  } else {# mult.samples = TRUE

    x_data <- stack(york_fit$data$x)[, 1]
    y_data <- stack(york_fit$data$y)[, 1]

    ddf <- data.frame(x = x_data, y = y_data)

    # plot fitted line and data in multiple sample case
    plot_1 <- ggplot(data = ddf, aes(x = x,
                                     y = y)) +
      geom_point() +
      geom_abline(aes(slope = york_fit$coefficients[2, 1],
                      intercept = york_fit$coefficients[1, 1]), col = "red") +
      labs(title = "York's best-fit straight line",
           x = "x data", y = "y data") +
      theme(plot.title = element_text(hjust = 0.5))
  }


  ### Compare York's best fit line with OLS and orthogonal
  mean_x <- mean(x_data)
  mean_y <- mean(y_data)
  ddf2 <- data.frame(x = x_data, y = y_data)

  # calculate orthogonal fit (weights = 1)
  orthogonal <- york(x = x_data, y = y_data, weights_x = 1,
                     weights_y = 1, r_xy_errors = 0,
                     tolerance = york_fit$york_arguments$tolerance,
                     max_iterations = york_fit$york_arguments$max_iterations)

  # plot york, orthogonal and OLS line
  plot_2 <- ggplot(data = ddf2, aes(x = x,
                                    y = y)) +
    geom_point() +
    geom_abline(aes(slope = york_fit$coefficients[2, 1],
                    intercept = york_fit$coefficients[1, 1], colour = "York"),
                key_glyph = draw_key_rect) +
    geom_abline(aes(slope = york_fit$ols_summary$coefficients_ols[2, 1],
                    intercept = york_fit$ols_summary$coefficients_ols[1, 1],
                    colour = "OLS")) +
    geom_abline(aes(slope = orthogonal$coefficients[2, 1], intercept =
                      orthogonal$coefficients[1, 1], colour = "Orthogonal")) +
    geom_vline(xintercept = york_fit$weighted_mean_x, linetype = "dashed",
               color = "red", size = 0.4) +
    geom_hline(yintercept = york_fit$weighted_mean_y, linetype = "dashed",
               color = "red", size = 0.4) +
    geom_vline(xintercept = mean_x,
               linetype = "dashed", color = "blue", size = 0.4) +
    geom_hline(yintercept = mean_y,
               linetype = "dashed", color = "blue", size = 0.4) +
    geom_point(aes(x = mean_x,
                   y = mean_y), col = "blue") +
    geom_point(aes(x = york_fit$weighted_mean_x,
                   y = york_fit$weighted_mean_y), col = "red") +
    scale_colour_manual(values = c("blue", "green", "red")) +
    labs(title = "York' best fit straight line compared to OLS and orthogonal ",
         x = "x data", y = "y data", colour = "") +
    theme(plot.title = element_text(hjust = 0.5))


  ### x y residuals plotted against each other
  ddf3 <- data.frame(x = as.numeric(york_fit$x_residuals),
                     y = as.numeric(york_fit$y_residuals))

  # plot x and y residuals
  plot_3 <- ggplot(aes(x = x, y = y), data = ddf3) +
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed", col = "red") +
    geom_vline(xintercept = 0, linetype = "dashed", col = "red") +
    labs( title = "y residuals vs. x residuals",
          x = "x residuals", y = "y residuals") +
    theme(plot.title = element_text(hjust = 0.5))


  ### x resid vs fitted y
  ddf4 <- data.frame(x = as.numeric(york_fit$fitted_y),
                     y = as.numeric(york_fit$x_residuals))

  # plot x resid vs fitted y
  plot_4 <- ggplot(aes(x = x, y = y), data = ddf4) +
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed", col = "red") +
    labs(title = "x residuals vs. fitted y",
          x = "fitted y", y = "x residuals") +
    theme(plot.title = element_text(hjust = 0.5))


  ### y resid vs fitted y
  ddf5 <- data.frame(x = as.numeric(york_fit$fitted_y),
                     y = as.numeric(york_fit$y_residuals))

  # plot y resid vs fitted y
  plot_5 <- ggplot(aes(x = x, y = y), data = ddf5) +
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed", col = "red") +
    labs(title = "y residuals vs. fitted y",
         x = "fitted y", y = "y residuals") +
    theme(plot.title = element_text(hjust = 0.5))


  ### detect influential observations
  ###  -> assume asymptotic normality of slope coefficients)
  if (york_fit$york_arguments$mult_samples == FALSE) {
    sd_x <- york_fit$data[, 3]
    sd_y <- york_fit$data[, 4]
    r_xy_errors <- york_fit$data[, 5]
    x_mult <- NULL
    y_mult <- NULL
  } else {# mult_sample = TRUE
    sd_x <- york_fit$data$sd_x
    sd_y <- york_fit$data$sd_y
    r_xy_errors <- york_fit$data$error_correlation
    x_mult <- york_fit$data$x
    y_mult <- york_fit$data$y
  }
  influential <- f_p_influential(x_data, y_data,
                                 york_fit$york_arguments$mult_samples,
                                 sd_x, sd_y,
                                 r_xy_errors,
                                 york_fit$coefficients[2, 1],
                                 york_fit$coefficients[2, 2],
                                 york_fit$york_arguments$tolerance,
                                 york_fit$york_arguments$max_iterations,
                                 x_mult, y_mult)


  ddf6 <- data.frame(x = x_data, y = y_data)
  ddf_influential <- data.frame(influential$detect_influential)

  # plot fitted line and observations with encircled influential points
  plot_6 <- ggplot(data = ddf6, aes(x = x,
                                      y = y)) +
    geom_point() +
    geom_abline(aes(slope = york_fit$coefficients[2, 1],
                    intercept = york_fit$coefficients[1, 1]), col = "red") +
    geom_point(aes(x = ddf_influential$x,
                   y = ddf_influential$y,
                   col = ddf_influential$which_row),
               size = 15, shape = 1,
               data = ddf_influential) +
    labs(title = influential$pl_title,
         subtitle = influential$pl_subtitle,
         x = "x data", y = "y data", col = influential$leg_title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.subtitle = element_text(hjust = 0.5))


  ### Chi^2 test graph and test results
  chisq_test <- f_p_chisq_test(york_fit$goodness_of_fit$chisq_df,
                               york_fit$goodness_of_fit$p_value,
                               york_fit$goodness_of_fit$chisq_statistic)
  ddf_reject_region <- data.frame(chisq_test$reject_region)

  # plot chisq graph and test results
  plot_7 <- ggplot(data = data.frame(x = c(0, chisq_test$x_limit)), aes(x)) +
    stat_function(fun = dchisq, n = 1e+3, args = list(df = chisq_test$df),
                  size = 1.1) +
    geom_hline(yintercept = 0, col = "black") +
    geom_vline(xintercept = 0, col = "black") +
    geom_segment(aes(x = chisq_test$x_p_value,
                     y = 0,
                     xend = chisq_test$x_p_value,
                     yend = chisq_test$y_p_value,
                     colour = "test statistic")) +
    geom_segment(aes(x = chisq_test$critical_value[1],
                     y =  0,
                     xend = chisq_test$critical_value[1],
                     yend = chisq_test$y_critical[1],
                     colour = "critical value, alpha = 0.1"),
                 alpha = 0.8) +
    geom_segment(aes(x = chisq_test$critical_value[2],
                     y =  0,
                     xend = chisq_test$critical_value[2],
                     yend = chisq_test$y_critical[2],
                     colour = "critical value, alpha = 0.05")) +
    geom_segment(aes(x = chisq_test$critical_value[3],
                     y =  0,
                     xend = chisq_test$critical_value[3],
                     yend = chisq_test$y_critical[3],
                     colour = "critical value, alpha = 0.01")) +
    geom_segment(aes(x = ddf_reject_region$x,
                     y = ddf_reject_region$y,
                     xend = ddf_reject_region$x,
                     yend = ddf_reject_region$yend),
                 colour = "red", alpha = 0.05,
                 size = (15 - log(chisq_test$df) * 2.5),
                 data = ddf_reject_region) +
    scale_colour_manual("", values = c("red4", "red3", "red", "lightgreen")) +
    annotate("text", label = paste("test statistic = ",
                                   chisq_test$test_statistic),
             x =  Inf, y = Inf,
             vjust = 7.5, hjust = chisq_test$adj, size = 6) +
    annotate("text", label = paste("p-value", chisq_test$math,
                                   chisq_test$p_value ),
              x =  Inf, y = Inf,
             vjust = 10, hjust = 1.3, size = 6) +
    labs(title = paste("Chi-squared test: ",
                       york_fit$goodness_of_fit$test_result),
         x = "x", y = "Density") +
    theme(plot.title = element_text(hjust = 0.5))


  ### Trace plot -> Convergence of slope
  if (york_fit$york_arguments$approx_solution == FALSE) {
    ddf8 <- data.frame(x = c(0, 1:york_fit$n_iterations),
                       y = c(york_fit$ols_summary$coefficients_ols[2, 1],
                             york_fit$slope_per_iteration[,1]))

    # plot trace plot
    plot_8 <- ggplot(aes(x = x, y = y), data = ddf8) +
      geom_line() +
      geom_point() +
      geom_hline(yintercept = york_fit$coefficients[2,1], col = "darkblue") +
      labs(title = "Trace plot",
           x = "Number of iterations", y = "slope coefficient") +
      theme(plot.title = element_text(hjust = 0.5))
  } else {# approx_solution = TRUE
    plot_8 <- NULL
  }
  return(list(plot_1, plot_2, plot_3, plot_4, plot_5, plot_6, plot_7, plot_8))
}

