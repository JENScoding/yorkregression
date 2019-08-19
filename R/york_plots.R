york.plots <- function(york.output){
  plot(york.output$original.x.values, york.output$original.y.values, pch = 16)
  lines(york.output$original.x.values, york.output$fitted_y , col = "red",lwd = 2)
  lines(york.output$original.x.values, york.output$fitted_ols, col = "blue", lty = "dashed",lwd = 2)
  legend("topright",legend = c("OLS","York"), fill = c("blue","red"))
  #the York regression line and the OLS regression line both go to the "center of gravity"
  abline(v = york.output$mean.x, h = york.output$mean.y, lty = "dashed", col = "red")
  points(x = york.output$mean.x, y = york.output$mean.y, col = "red", pch = 16)
  abline(v = mean(york.output$original.x.values), h = mean(york.output$original.y.values) ,lty = "dashed", col = "blue")
  points(x = mean(york.output$original.x.values), y = mean(york.output$original.y.values), col = "blue", pch = 16)
  compare_OLS_York <- recordPlot()
  dev.off()
  
  plot(york.output$x.residuals, york.output$y.residuals)
  residuals <- recordPlot()
  dev.off()
  
  return(list("compare_OLS_York_with_center_of_gravity" = compare_OLS_York, "residuals" = residuals))
}

my_plots <- york.plots(york.output)

my_plots$compare_OLS_York_with_center_of_gravity
