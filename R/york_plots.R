load("R/original_data.RData")

library(ggplot2)
york.plots <- function(x, y, tolerance = 1e-10, weights.x = NULL, weights.y = NULL,
                       r.xy = NULL, sd.x = NULL, sd.y = NULL, mult.samples = F) {
  york.output <- york(x, y, tolerance, weights.x, weights.y,
                      r.xy, sd.x, sd.y)
  ddf <- data.frame(x=x,y=y)
  plot.1 <- ggplot(data=ddf, aes(x=york.output$original.x.values,
                                 y=york.output$original.y.values)) +
    geom_abline(aes(slope = york.output$coefficients.york[2,1],
                    intercept = york.output$coefficients.york[1,1],colour="York"),
                key_glyph = draw_key_rect) +
    geom_abline(aes(slope = york.output$coefficients.ols[2,1],
                    intercept = york.output$coefficients.ols[1,1],colour="OLS")) +
    geom_abline(aes(slope = york.output$coefficients.mayor[2,1], intercept =
                      york.output$coefficients.mayor[1,1],colour="Orthogonal")) +
    labs(colour="") + scale_colour_manual(values=c("blue","green","red")) +
    geom_point() +
    labs(title="York' best fit straight line compared to OLS and orthogonal ",
         x ="x data", y = "y data") +
    geom_vline(xintercept = york.output$mean.x, linetype="dashed",
               color = "red", size=0.4) +
    geom_hline(yintercept = york.output$mean.y, linetype="dashed",
               color = "red", size=0.4) +
    geom_vline(xintercept = mean(york.output$original.x.values),
               linetype="dashed",color = "blue", size = 0.4) +
    geom_hline(yintercept = mean(york.output$original.y.values),
               linetype="dashed",color = "blue", size = 0.4) +
    geom_smooth(method ="lm") +
    geom_point(aes(x = mean(york.output$original.x.values),
                   y= mean(york.output$original.y.values)), col = "blue") +
    geom_point(aes(x = york.output$mean.x,
                   y= york.output$mean.y), col = "red") +
    theme(plot.title =element_text(hjust = 0.5))
  ddf2 <- data.frame( x = york.output$x.residuals,  y = york.output$y.residuals)
  plot.2 <- ggplot(aes(x = x, y = y), data = ddf2) +
    geom_point() +
    labs( title = "y residuals vs. x residuals",
          x = "x residuals", y = "y residuals") +
    geom_hline(yintercept = 0, linetype = "dashed", col = "red") +
    theme(plot.title =element_text(hjust = 0.5))

  ddf3 <- data.frame( x = york.output$fitted.y, y = york.output$x.residuals)
  plot.3 <- ggplot(aes(x = x, y = y), data = ddf3) +
    geom_point() +
    labs( title = "x residuals vs. fitted y",
          x = "fitted y", y = "x residuals") +
    geom_hline(yintercept = 0, linetype = "dashed", col = "red") +
    theme(plot.title = element_text(hjust = 0.5))
  ddf4 <- data.frame(x = york.output$fitted.y, y = york.output$y.residuals)
  plot.4 <- ggplot(aes(x=x, y=y), data = ddf4) +
    geom_point() +
    labs(title = "y residuals vs. fitted y",
         x = "fitted y", y= "y residuals") +
    geom_hline(yintercept = 0, linetype = "dashed", col = "red") +
    theme(plot.title =element_text(hjust = 0.5))
  ddf5 <- data.frame(x = 1:york.output$number.of.iterations, y = york.output$slope.after.each.iteration[,1])
  plot.5 <- ggplot(aes(x=x, y=y), data = ddf5)+
    geom_line() +
    geom_point() +
    labs(title = "Trace plot",
         x = "Number of iterations", y = "slope coefficient") +
    geom_hline(yintercept = york.output$coefficients.york[2,1], col = "darkblue") +
    theme(plot.title = element_text(hjust = 0.5))
  return(list(plot.1, plot.2, plot.3, plot.4, plot.5))
}
york.plots(x=x,y=y, weights.x =weights.x, weights.y=weights.y,r.xy=0.1)
