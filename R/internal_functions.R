#'
#'

calc.var <- function(x, mean.x) {
  sum((x - mean.x)^2) / (ncol(x) - 1)
}
calc.corr <- function(x, y, mean.x, mean.y) {
  sum((x - mean.x) * (y - mean.y)) /
    sqrt(calc.var(x, mean.x) * calc.var(y, mean.y))
}
