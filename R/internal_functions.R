#'
#'

# Make data suitable for Algorithm
f_rewrite <- function(x, y, weights.x = NULL, weights.y = NULL,
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

# functions for the variance and correlation of the row
f_var_row <- function(x) {
  sum((x - apply(x, 1, mean))^2) / (length(x) - 1)
}
f_corr_row <- function(x, y) {
  sum((x - apply(x, 1, mean)) * (y - apply(y, 1, mean))) /
    (sqrt(f_var_row(x) * f_var_row(y)) * (length(x) - 1))
}

# simple linear regression to get ols
f_ols_reg <- function(x, y) {

  # find ols slope
  x_input <- matrix(c(rep(1, length(x)), x), ncol =2)
  lm <- solve(t(x_input) %*% x_input) %*% t(x_input) %*% y
  slope <- lm[2]
  intercept <- lm[1]
  if (any(is.na(c(slope,intercept)))){
    stop("Cannot fit a line through these data!")
  }

  # additional ols model information
  fitted_y <- x_input %*% lm
  residuals <- y - fitted_y
  RSS <- sum(residuals^2)
  sigma_squared_hat <- (1 / (length(x) - 2)) * RSS
  se_of_reg <- sqrt(sigma_squared_hat)
  mean_x <- mean(x)
  mean_y <- mean(y)
  centered_x <- x - mean_x
  centered_y <- y - mean_y
  SS_x <- sum(centered_x^2)
  SS_y <- sum(centered_y^2)
  S_x <- sum(x^2)
  SS_xy <- sum((centered_x) * (centered_y))
  se_intercept <- sqrt(sigma_squared_hat * (S_x / (length(x) * SS_x)))
  se_slope <- sqrt(sigma_squared_hat / SS_x)
  r_squared <- 1 - RSS / SS_y
  r_squared_adjusted <- r_squared - (1 - r_squared)*
    (1 / (length(x)-2))
  f_statistic <- (r_squared / (1- r_squared)) * ((length(x)-2))

  coef <- matrix(c(intercept, slope, se_intercept,
                       se_slope),
                     nrow = 2)
  rownames(coef) <- c("intercept", "slope")
  colnames(coef) <- c("Estimate", "Std_Error")

  ols_output <- list("slope" = slope,
                     "se_slope" = se_slope,
                     "coefficients_ols" = coef,
                     "fitted_y_ols" = fitted_y,
                     "residuals_ols" = residuals,
                     "residual_sum_of_squares_ols" = RSS,
                     "total_sum_of_squares_ols" = SS_y,
                     "se_of_reg_ols" = se_of_reg,
                     "r_squared_ols" = r_squared,
                     "r_squared_adjusted_ols" = r_squared_adjusted,
                     "f_statistic_ols" = f_statistic)
}
