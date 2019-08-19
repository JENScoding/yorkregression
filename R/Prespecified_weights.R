### Prespecify weights

define.weights <- function(x, y, factor.x = 0.02, factor.y = 0.003){
  weights.x <- 1 / (x * factor.x)^2
  weights.y <- 1 / (y * factor.y)^2
  return(data.frame(weights.x, weights.y))
}