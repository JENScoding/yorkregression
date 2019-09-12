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
#' @param york.output An object of class "york"
#' @return york.plots returns 6 different diagnostic plots.
#'
#' @examples
#' # Example: York's regression with weight data taken from Pearson (1901):
#' x <- c(0.0, 0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4)
#' y <- c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5)
#' weights.x <- c(1e+3, 1e+3, 5e+2, 8e+2, 2e+2, 8e+1, 6e+1, 2e+1, 1.8, 1)
#' weights.y <- c(1, 1.8, 4, 8, 20, 20, 70, 70, 1e+2, 5e+2)
#' r.xy <- 0
#' york.output <- york(x = x, y = y, weights.x = weights.x, weights.y = weights.y,
#'                     r.xy = 0)
#' york.plots(york.output)
#' @name york.plots
#' @export
#' @importFrom ggplot2 ggplot aes geom_abline geom_point labs theme element_text
#' draw_key_rect scale_colour_manual geom_vline geom_hline geom_smooth geom_line
#' @importFrom utils stack
utils::globalVariables(c("x", "y"))
york.plots <- function(york.output) {
  if (class(york.output) != "york") {
    stop("Input must be of class york (Output of york function)")
  }
  if (york.output$york.arguments$mult.samples == F) {

    x.data <- york.output$data[, 1]
    y.data <- york.output$data[, 2]

    ddf <- data.frame(x = x.data, y = y.data)
    plot.1 <- ggplot(data=ddf, aes(x = x,
                                   y = y)) +
      geom_abline(aes(slope = york.output$coefficients.york[2, 1],
                      intercept = york.output$coefficients.york[1, 1]), col = "red") +
      geom_point() +
      labs(title="York's best-fit straight line with correlated errors",
           x ="x data", y = "y data") +
      theme(plot.title =element_text(hjust = 0.5))
  } else {

    x.data <- york.output$data$mean.x.i
    y.data <- york.output$data$mean.y.i

    x.data.1 <- stack(york.output$data$x)[, 1]
    y.data.1 <- stack(york.output$data$y)[, 1]

    ddf <- data.frame(x = x.data.1, y = y.data.1)
    plot.1 <- ggplot(data = ddf, aes(x = x,
                                   y = y)) +
      geom_abline(aes(slope = york.output$coefficients.york[2, 1],
                      intercept = york.output$coefficients.york[1, 1]), col = "red") +
      geom_point() +
      labs(title="York's best-fit straight line with correlated errors",
           x ="x data", y = "y data") +
      theme(plot.title =element_text(hjust = 0.5))
  }

  ddf2 <- data.frame(x = x.data, y = y.data)
  plot.2 <- ggplot(data = ddf2, aes(x = x,
                                 y = y)) +
    geom_abline(aes(slope = york.output$coefficients.york[2, 1],
                    intercept = york.output$coefficients.york[1, 1], colour="York"),
                key_glyph = draw_key_rect) +
    geom_abline(aes(slope = york.output$coefficients.ols[2, 1],
                    intercept = york.output$coefficients.ols[1, 1], colour="OLS")) +
    geom_abline(aes(slope = york.output$coefficients.orthogonal[2, 1], intercept =
                      york.output$coefficients.orthogonal[1, 1], colour = "Orthogonal")) +
    labs(colour="") +
    scale_colour_manual(values=c("blue", "green", "red")) +
    geom_point() +
    labs(title="York' best fit straight line compared to OLS and orthogonal ",
         x ="x data", y = "y data") +
    geom_vline(xintercept = york.output$weighted.mean.x, linetype = "dashed",
               color = "red", size = 0.4) +
    geom_hline(yintercept = york.output$weighted.mean.y, linetype = "dashed",
               color = "red", size = 0.4) +
    geom_vline(xintercept = mean(x.data),
               linetype = "dashed", color = "blue", size = 0.4) +
    geom_hline(yintercept = mean(y.data),
               linetype= "dashed", color = "blue", size = 0.4) +
    geom_smooth(method = "lm") +
    geom_point(aes(x = mean(x.data),
                   y = mean(y.data)), col = "blue") +
    geom_point(aes(x = york.output$weighted.mean.x,
                   y = york.output$weighted.mean.y), col = "red") +
    theme(plot.title = element_text(hjust = 0.5))

  ddf3 <- data.frame(x = york.output$x.residuals,  y = york.output$y.residuals)
  plot.3 <- ggplot(aes(x = x, y = y), data = ddf3) +
    geom_point() +
    labs( title = "y residuals vs. x residuals",
          x = "x residuals", y = "y residuals") +
    geom_hline(yintercept = 0, linetype = "dashed", col = "red") +
    geom_vline(xintercept = 0, linetype = "dashed", col = "red") +
    theme(plot.title = element_text(hjust = 0.5))

  ddf4 <- data.frame( x = york.output$fitted.y, y = york.output$x.residuals)
  plot.4 <- ggplot(aes(x = x, y = y), data = ddf4) +
    geom_point() +
    labs( title = "x residuals vs. fitted y",
          x = "fitted y", y = "x residuals") +
    geom_hline(yintercept = 0, linetype = "dashed", col = "red") +
    theme(plot.title = element_text(hjust = 0.5))

  ddf5 <- data.frame(x = york.output$fitted.y, y = york.output$y.residuals)
  plot.5 <- ggplot(aes(x = x, y = y), data = ddf5) +
    geom_point() +
    labs(title = "y residuals vs. fitted y",
         x = "fitted y", y = "y residuals") +
    geom_hline(yintercept = 0, linetype = "dashed", col = "red") +
    theme(plot.title = element_text(hjust = 0.5))

  if (york.output$york.arguments$approx.solution == F) {
    ddf6 <- data.frame(x = 1:york.output$number.of.iterations, y = york.output$slope.after.each.iteration[,1])
    plot.6 <- ggplot(aes(x = x, y = y), data = ddf6)+
      geom_line() +
      geom_point() +
      labs(title = "Trace plot",
           x = "Number of iterations", y = "slope coefficient") +
      geom_hline(yintercept = york.output$coefficients.york[2,1], col = "darkblue") +
      theme(plot.title = element_text(hjust = 0.5))
  } else {
    plot.6 <- NULL
  }
  return(list(plot.1, plot.2, plot.3, plot.4, plot.5, plot.6))
}

