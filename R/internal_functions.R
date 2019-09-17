#'
#'

rewrite <- function(x, y, weights.x = NULL, weights.y = NULL,
                    sd.x = NULL, sd.y = NULL, r.xy = NULL) {

  # if input was only 1 value repeat it to adjust to sample size
  if (length(weights.x) == 1) {
    weights.x = rep(weights.x, length(x))
  }
  if (length(weights.y) == 1) {
    weights.y = rep(weights.y, length(y))
  }
  if (length(sd.x) == 1) {
    sd.x = rep(sd.x, length(x))
  }
  if (length(sd.y) == 1) {
    sd.y = rep(sd.y, length(y))
  }
  if (length(r.xy) == 1) {
    r.xy = rep(r.xy, length(x))
  }

  # specify weights and standard errors
  if (all(sapply(list(sd.x, sd.y, weights.x, weights.y),
                 function(x) !is.null(x)))) {
    stop("You can't specify weights and standard errors at the same time!")
  }
  if (is.null(weights.x) & is.null(weights.y)) {
    weights.x <- 1 / sd.x^2
    weights.y <- 1 / sd.y^2
  }
  if (is.null(sd.x) & is.null(sd.y)) {
    sd.x <- 1 / sqrt(weights.x)
    sd.y <- 1 / sqrt(weights.y)
  }

  # delete rows with NA values
  omit_na  <- c(which(is.na(x)), which(is.na(y)),
                  which(is.na(weights.x)), which(is.na(weights.y)),
                  which(is.na(sd.x)), which(is.na(sd.y)))
  if (length(omit_na) > 0){
    y <- y[-omit_na]
    x <- x[-omit_na]
    weights.x <- weights.x[-omit_na]
    weights.y <- weights.y[-omit_na]
    sd.x <- sd.x[-omit_na]
    sd.y <- sd.y[-omit_na]
    r.xy <- r.xy[-omit_na]
    omitted.share <- length(omit_na) / length(x)
    if (omitted.share > 0.1) {
      warning(omitted.share * 100,
              "% of the data were removed due to missing values!")
    }
  }

  input <- list("x" = x,
                "y" = y,
                "weights.x" = weights.x,
                "weights.y" = weights.y,
                "sd.x" = sd.x,
                "sd.y" = sd.y,
                "r.xy" = r.xy)
  return(input)
}

calc.var <- function(x) {
  sum((x - apply(x, 1, mean))^2) / (length(x) - 1)
}
calc.corr <- function(x, y, mean.x, mean.y) {
  sum((x - apply(x, 1, mean)) * (y - apply(y, 1, mean))) /
    (sqrt(calc.var(x) * calc.var(y)) * (length(x) - 1))
}

