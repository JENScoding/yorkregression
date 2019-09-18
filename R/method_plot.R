#' @title Plots for York's regression
#'
#' @description Gives several diagnostic plots for York's regression.
#' @details Given an object of class "york" this function gives several
#' dianostic plots. The first plot shows York's best-fit straight line only. The
#' second plot shows ork' best fit straight line compared to OLS and orthogonal
#' regression. Additionally, the different "center of gravity" are shown. The
#' third plot is a "y-residuals vs. x-residuals plot". The fourth is a
#' "x-residuals vs. fitted-y plot". The fith plot is a "y-residuals vs.
#' fitted y" plot and the last plot is a "trace plot" which shows the slope after
#' each interation, i.e. how the slope coefficient converges.
#' @param x An object of class "york"
#' @param ... Arguments to be passed to methods.
#'
#' @return The plot function returns 7 different plots.
#'
#' @examples
#' # Example: York's regression with weight data taken from Pearson (1901):
#' x <- c(0.0, 0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4)
#' y <- c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5)
#' weights_x <- c(1e+3, 1e+3, 5e+2, 8e+2, 2e+2, 8e+1, 6e+1, 2e+1, 1.8, 1)
#' weights_y <- c(1, 1.8, 4, 8, 20, 20, 70, 70, 1e+2, 5e+2)
#' r_xy_errors <- 0
#' york_fit <- york(x = x, y = y, weights_x = weights_x, weights_y = weights_y,
#'                     r_xy_errors = 0)
#' plot(york_fit)
#'
#' @name plot.york
#'
#' @importFrom ggplot2 ggplot aes geom_abline geom_point labs theme element_text
#' draw_key_rect scale_colour_manual geom_vline geom_hline geom_smooth geom_line
#' @importFrom utils stack
#' @importFrom stats pnorm
#'
#' @method plot york
#' @S3method plot york

plot.york <- function(x, ...) {
  if (class(x) != "york") {
    stop("Input must be of class york (Output of york function)")
  }

  # Rename x and make variable y global
  york_fit <- x
  y <- NULL

  # York's best fit line with data points
  if (york_fit$york_arguments$mult_samples == FALSE) {

    x_data <- york_fit$data[, 1]
    y_data <- york_fit$data[, 2]

    ddf <- data.frame(x = x_data, y = y_data)
    plot_1 <- ggplot(data = ddf, aes(x = x,
                                   y = y)) +
      geom_abline(aes(slope = york_fit$coefficients[2, 1],
                      intercept = york_fit$coefficients[1, 1]), col = "red") +
      geom_point() +
      labs(title="York's best-fit straight line",
           x ="x data", y = "y data") +
      theme(plot.title =element_text(hjust = 0.5))
  } else {

    x_data <- stack(york_fit$data$x)[, 1]
    y_data <- stack(york_fit$data$y)[, 1]

    ddf <- data.frame(x = x_data, y = y_data)
    plot_1 <- ggplot(data = ddf, aes(x = x,
                                     y = y)) +
      geom_abline(aes(slope = york_fit$coefficients[2, 1],
                      intercept = york_fit$coefficients[1, 1]), col = "red") +
      geom_point() +
      labs(title="York's best-fit straight line",
           x ="x data", y = "y data") +
      theme(plot.title = element_text(hjust = 0.5))
  }

  # Compare York's best fit line with olse and orthogonal
  mean_x <- mean(x_data)
  mean_y <- mean(y_data)
  orthogonal <- york(x = x_data, y = y_data, weights_x = 1,
                     weights_y = 1, r_xy_errors = 0,
                     tolerance = york_fit$york_arguments$tolerance,
                     max_iterations = york_fit$york_arguments$max_iterations)
  ddf2 <- data.frame(x = x_data, y = y_data)
  plot_2 <- ggplot(data = ddf2, aes(x = x,
                                    y = y)) +
    geom_abline(aes(slope = york_fit$coefficients[2, 1],
                    intercept = york_fit$coefficients[1, 1], colour = "York"),
                key_glyph = draw_key_rect) +
    geom_abline(aes(slope = york_fit$ols_summary$coefficients_ols[2, 1],
                    intercept = york_fit$ols_summary$coefficients_ols[1, 1],
                    colour = "OLS")) +
    geom_abline(aes(slope = orthogonal$coefficients[2, 1], intercept =
                      orthogonal$coefficients[1, 1], colour = "Orthogonal")) +
    labs(colour = "") +
    scale_colour_manual(values=c("blue", "green", "red")) +
    geom_point() +
    labs(title="York' best fit straight line compared to OLS and orthogonal ",
         x ="x data", y = "y data") +
    geom_vline(xintercept = york_fit$weighted_mean_x, linetype = "dashed",
               color = "red", size = 0.4) +
    geom_hline(yintercept = york_fit$weighted_mean_y, linetype = "dashed",
               color = "red", size = 0.4) +
    geom_vline(xintercept = mean_x,
               linetype = "dashed", color = "blue", size = 0.4) +
    geom_hline(yintercept = mean_y,
               linetype= "dashed", color = "blue", size = 0.4) +
    geom_point(aes(x = mean_x,
                   y = mean_y), col = "blue") +
    geom_point(aes(x = york_fit$weighted_mean_x,
                   y = york_fit$weighted_mean_y), col = "red") +
    theme(plot.title = element_text(hjust = 0.5))

  # X Y residuals plotted against each other
  ddf3 <- data.frame(x = as.numeric(york_fit$x_residuals),
                     y = as.numeric(york_fit$y_residuals))
  plot_3 <- ggplot(aes(x = x, y = y), data = ddf3) +
    geom_point() +
    labs( title = "y residuals vs. x residuals",
          x = "x residuals", y = "y residuals") +
    geom_hline(yintercept = 0, linetype = "dashed", col = "red") +
    geom_vline(xintercept = 0, linetype = "dashed", col = "red") +
    theme(plot.title = element_text(hjust = 0.5))

  # x resid vs fitted y
  ddf4 <- data.frame(x = as.numeric(york_fit$fitted_y),
                     y = as.numeric(york_fit$x_residuals))
  plot_4 <- ggplot(aes(x = x, y = y), data = ddf4) +
    geom_point() +
    labs( title = "x residuals vs. fitted y",
          x = "fitted y", y = "x residuals") +
    geom_hline(yintercept = 0, linetype = "dashed", col = "red") +
    theme(plot.title = element_text(hjust = 0.5))

  # y resid vs fitted y
  ddf5 <- data.frame(x = as.numeric(york_fit$fitted_y),
                     y = as.numeric(york_fit$y_residuals))
  plot_5 <- ggplot(aes(x = x, y = y), data = ddf5) +
    geom_point() +
    labs(title = "y residuals vs. fitted y",
         x = "fitted y", y = "y residuals") +
    geom_hline(yintercept = 0, linetype = "dashed", col = "red") +
    theme(plot.title = element_text(hjust = 0.5))

  # detect outliers, assume asymptotic normality of slope coefficients
  if (york_fit$york_arguments$mult_samples == FALSE) {
    slope_outlier <- NULL
    for (i in 1:length(x_data)) {
      slope_outlier[i] <- york(x_data[-i], y_data[-i],
                               sd_x = york_fit$data[-i, 3],
                               sd_y = york_fit$data[-i, 4],
                               r_xy_errors = york_fit$data[-i, 5],
                               tolerance =
                                 york_fit$york_arguments$tolerance,
                               max_iterations =
                                 york_fit$york_arguments$max_iterations)[[1]][2, 1]
    }
    detect_outlier1 <- (york_fit$coefficients[2, 1] - slope_outlier) /
      york_fit$coefficients[2, 2]
    detect_outlier2 <- 2 * (1 - pnorm(abs(detect_outlier1)))
    detect_outlier3 <- which(detect_outlier2 <= 0.01)
    detect_outlier4 <- x_data[detect_outlier3]
  } else {
    # multiple samples
    x_data <- york_fit$data$x
    y_data <- york_fit$data$y
    slope_outlier <- NULL
    for (i in 1:nrow(x_data)) {
      slope_outlier[i] <- suppressWarnings(
        york(x_data[-i,], y_data[-i,],
             tolerance =
               york_fit$york_arguments$tolerance,
             max_iterations =
               york_fit$york_arguments$max_iterations,
             mult_samples =
               york_fit$york_arguments$mult_samples)[[1]][2, 1])
    }
    detect_outlier1 <- (york_fit$coefficients[2, 1] - slope_outlier) /
      york_fit$coefficients[2, 2]
    detect_outlier2 <- 2 * (1 - pnorm(abs(detect_outlier1)))
    detect_outlier3 <- which(detect_outlier2 <= 0.01)
  }


  # ddf6 <- data.frame(x = x_data, y = y_data)
  # plot_6 <- ggplot(aes(x = x, y = y),
  #                       data = ddf6) +
  #   geom_abline(aes(slope = x$coefficients[2, 1],
  #                   intercept = x$coefficients[1, 1]), col = "red") +
  #   geom_point() +
  #   geom_point(aes(x = detect_outlier4[1],
  #                  y = y_data[detect_outlier3][1]), col = "blue", fill = "white",
  #              shape = "O", size = 3) +
  #   labs(title = "York's best-fit straight line with correlated errors",
  #        x = "x data", y = "y data") +
  #   theme(plot.title =element_text(hjust = 0.5))
  plot_6 <- NULL

  # Trace plot -> Convergence of slope
  if (york_fit$york_arguments$approx_solution == FALSE) {
    ddf7 <- data.frame(x = 1:york_fit$n_iterations,
                       y = york_fit$slope_per_iteration[,1])
    plot_7 <- ggplot(aes(x = x, y = y), data = ddf7)+
      geom_line() +
      geom_point() +
      labs(title = "Trace plot",
           x = "Number of iterations", y = "slope coefficient") +
      geom_hline(yintercept = york_fit$coefficients[2,1], col = "darkblue") +
      theme(plot.title = element_text(hjust = 0.5))
  } else {
    plot_7 <- NULL
  }
  return(list(plot_1, plot_2, plot_3, plot_4, plot_5, plot_6, plot_7))
}

