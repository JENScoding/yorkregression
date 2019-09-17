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

  # matrix for coefficients and their se
  coef <- matrix(c(intercept, slope, se_intercept,
                       se_slope),
                     nrow = 2)
  rownames(coef) <- c("intercept", "slope")
  colnames(coef) <- c("Estimate", "Std_Error")

  # define output
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

  return(ols_output)
}

# Solve cubic root problem (approximate solution for slope from York (1966)
f_cubic_root <- function(x, y, weights.x, weights.y, r.xy, slope_ols, se_slope_ols) {

  if (any(r.xy != 0)) {
    stop(paste("There is no approximate solution in case of correlation",
               "between x and y errors!", sep = " "))
  }
  # use ols slope as intitial slope value and estimate weight and
  # centered data
  alpha <- sqrt(weights.x * weights.y)
  Weight <- alpha^2 / (slope_ols^2 * weights.y + weights.x)
  Weight_sum <- sum(Weight)
  x_bar <- sum(Weight * x) / Weight_sum
  y_bar <- sum(Weight * y) / Weight_sum
  x_centered <- x - x_bar
  y_centered <- y - y_bar

  # calculate alpha_cubic, beta_cubic and gamma_cubic. See York (1966) p. 1084
  xy <- x_centered * y_centered
  xW_w <- x_centered^2 * Weight^2 / weights.x
  yW_w <- y_centered^2 * Weight^2 / weights.x
  alpha_cubic <- 2 * sum(xy * Weight^2 / weights.x) /
    (3 * sum(xW_w))
  beta_cubic <- (sum(yW_w) - sum(Weight * x_centered^2)) /
    (3 * sum(xW_w))
  gamma_cubic <- - sum(xy * Weight) / (sum(x_centered^2 *
                                             Weight^2 / weights.x))

  # determine phi and solve cubic equation
  phi <- acos((alpha_cubic^3 - 3 /2 * alpha_cubic * beta_cubic + 0.5 *
                 gamma_cubic) /
                (alpha_cubic^2 - beta_cubic)^(3 / 2))

  sol_cubic <- alpha_cubic + 2 * (alpha_cubic^2 - beta_cubic)^0.5 *
    cos( 1 / 3 *(phi + 2 * pi * c(0:2)))

  # pick the root that is closest to ols
  ols_range <- c(slope_ols - 4 * se_slope_ols,
                 slope_ols + 4 * se_slope_ols)
  pick_right_root <- which(sol_cubic >= ols_range[1] &
                             sol_cubic <= ols_range[2])
  slope <- sol_cubic[pick_right_root]
  if (length(slope) == 0 | length(slope) > 1) {
    stop("An approximate solution does not exist!")
  }

  # later needed
  beta <- Weight * (x_centered / weights.y) + (slope * y_centered / weights.x)

  # define output
  cubic_root <- list("slope" = slope,
                     "Weight" = Weight,
                     "Weight_sum" = Weight_sum,
                     "x_bar" = x_bar,
                     "y_bar" = y_bar,
                     "beta" = beta)

  return(cubic_root)
}
