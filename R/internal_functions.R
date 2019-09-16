#'
#'

calc.var <- function(x) {
  sum((x - apply(x, 1, mean))^2) / (length(x) - 1)
}
calc.corr <- function(x, y, mean.x, mean.y) {
  sum((x - apply(x, 1, mean)) * (y - apply(y, 1, mean))) /
    (sqrt(calc.var(x) * calc.var(y)) * (length(x) - 1))
}

