york.plots <- function(york.output){
  plot(york.output$original.x.values, york.output$original.y.values,
       xlab = "x data",
       ylab = "y data",
       main = "Pearson's data with York's weights: Best fit straight line",
       col = "black",
       pch = 16,
       las = 1,
       ylim = c(1.5,6))

  lines(york.output$original.x.values, york.output$fitted.y , col = "red",lwd = 2)
  lines(york.output$original.x.values, york.output$fitted.ols, col = "blue", lty = "dashed",lwd = 2)
  legend("topright",legend = c("OLS","York"), fill = c("blue","red"))
  #the York regression line and the OLS regression line both go to the "center of gravity"
  abline(v = york.output$mean.x, h = york.output$mean.y, lty = "dashed", col = "red")
  points(x = york.output$mean.x, y = york.output$mean.y, col = "red", pch = 16)
  abline(v = mean(york.output$original.x.values), h = mean(york.output$original.y.values) ,lty = "dashed", col = "blue")
  points(x = mean(york.output$original.x.values), y = mean(york.output$original.y.values), col = "blue", pch = 16)
  compare.OLS.York <- recordPlot()
  dev.off()

  plot(york.output$x.residuals, york.output$y.residuals,
       xlab = "x residuals",
       ylab = "y residuals",
       main = "x residuals vs. y residuals",
       col = "black",
       pch = 16,
       las = 1)
  abline(h = 0, lty = "dashed", col = "red")
  x.residuals.vs.y.residuals <- recordPlot()
  dev.off()

  plot(york.output$fitted.y, york.output$x.residuals,
       xlab = "fitted y",
       ylab = "x residuals",
       main = "x residuals vs. fitted plot",
       col = "black",
       pch = 16,
       las = 1)
  abline(h = 0, lty = "dashed", col = "red")

  x.residuals.vs.fitted.y <- recordPlot()
  dev.off()

  plot(york.output$fitted.y, york.output$y.residuals,
       xlab = "fitted y",
       ylab = "y residuals",
       main = "y residuals vs. fitted plot",
       col = "black",
       pch = 16,
       las = 1)

  y.residuals.vs.fitted.y <- recordPlot()
  dev.off()

  return(list("compare.OLS.York.with.center.of.gravity" = compare.OLS.York, "residuals" = residuals, "yresiduals.vs.fitted" = yresiduals.vs.fitted.))
}

my.plots <- york.plots(york.output)

my.plots$compare.OLS.York.with.center.of.gravity
my.plots$residuals

##end of relevant script

load("original_data.RData")

x <- matrix(x,ncol=1)
x
(P <- x%*%solve((t(x)%*%x))%*%t(x))

diag(P)






