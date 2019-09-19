#' @title Internal Functions I
#'
#' @description Functions used in york function. For Internal Functions I
#' the output is used for proceeding calculations or displaying. All the
#' Internal Functions I functions begin with f_ and then the name of the
#' function is given.
#'
#' @keywords internal
#'

# Make data suitable for Algorithm
f_rewrite <- function(x, y, weights_x = NULL, weights_y = NULL,
                    sd_x = NULL, sd_y = NULL, r_xy_errors = NULL) {

  # if input was only 1 value repeat it to adjust to sample size
  if (length(weights_x) == 1) {
    weights_x = rep(weights_x, length(x))
  }
  if (length(weights_y) == 1) {
    weights_y = rep(weights_y, length(y))
  }
  if (length(sd_x) == 1) {
    sd_x = rep(sd_x, length(x))
  }
  if (length(sd_y) == 1) {
    sd_y = rep(sd_y, length(y))
  }
  if (length(r_xy_errors) == 1) {
    r_xy_errors = rep(r_xy_errors, length(x))
  }

  # specify weights and standard errors
  if (all(sapply(list(sd_x, sd_y, weights_x, weights_y),
                 function(x) !is.null(x)))) {
    stop("You can't specify weights and standard errors at the same time!")
  }
  if (is.null(weights_x) & is.null(weights_y)) {
    weights_x <- 1 / sd_x^2
    weights_y <- 1 / sd_y^2
  }
  if (is.null(sd_x) & is.null(sd_y)) {
    sd_x <- 1 / sqrt(weights_x)
    sd_y <- 1 / sqrt(weights_y)
  }

  # delete rows with NA values
  omit_na  <- c(which(is.na(x)), which(is.na(y)),
                  which(is.na(weights_x)), which(is.na(weights_y)),
                  which(is.na(sd_x)), which(is.na(sd_y)))
  if (length(omit_na) > 0){
    x_original <- x
    y <- y[-omit_na]
    x <- x[-omit_na]
    weights_x <- weights_x[-omit_na]
    weights_y <- weights_y[-omit_na]
    sd_x <- sd_x[-omit_na]
    sd_y <- sd_y[-omit_na]
    r_xy_errors <- r_xy_errors[-omit_na]
    omitted.share <- length(omit_na) / length(x_original)
    if (omitted.share > 0.1) {
      warning(omitted.share * 100,
              "% of the data were removed due to missing values!")
    }
  }

  input <- list("x" = x,
                "y" = y,
                "weights_x" = weights_x,
                "weights_y" = weights_y,
                "sd_x" = sd_x,
                "sd_y" = sd_y,
                "r_xy_errors" = r_xy_errors)
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
f_cubic_root <- function(x, y, weights_x, weights_y, r_xy, slope_ols, se_slope_ols) {

  if (any(r_xy != 0)) {
    stop(paste("There is no approximate solution in case of correlation",
               "between x and y errors!", sep = " "))
  }
  # use ols slope as intitial slope value and estimate weight and
  # centered data
  alpha <- sqrt(weights_x * weights_y)
  Weight <- alpha^2 / (slope_ols^2 * weights_y + weights_x)
  Weight_sum <- sum(Weight)
  x_bar <- sum(Weight * x) / Weight_sum
  y_bar <- sum(Weight * y) / Weight_sum
  x_centered <- x - x_bar
  y_centered <- y - y_bar

  # calculate alpha_cubic, beta_cubic and gamma_cubic. See York (1966) p. 1084
  xy <- x_centered * y_centered
  xW_w <- x_centered^2 * Weight^2 / weights_x
  yW_w <- y_centered^2 * Weight^2 / weights_x
  alpha_cubic <- 2 * sum(xy * Weight^2 / weights_x) /
    (3 * sum(xW_w))
  beta_cubic <- (sum(yW_w) - sum(Weight * x_centered^2)) /
    (3 * sum(xW_w))
  gamma_cubic <- - sum(xy * Weight) / (sum(x_centered^2 *
                                             Weight^2 / weights_x))

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
  beta <- Weight * (x_centered / weights_y) + (slope * y_centered / weights_x)

  # define output
  cubic_root <- list("slope" = slope,
                     "Weight" = Weight,
                     "Weight_sum" = Weight_sum,
                     "x_bar" = x_bar,
                     "y_bar" = y_bar,
                     "alpha" = alpha,
                     "beta" = beta)

  return(cubic_root)
}


f_define_output <- function(intercept, slope, sigma_intercept, sigma_slope,
                            weights_x, weights_y, mult_samples, x, y, sd_x,
                            sd_y, r_xy_errors, x_original, y_original,
                            x_errors, y_errors, mean_x_i, mean_y_i,
                            slope_per_iteration, tolerance, max_iterations,
                            approx_solution, S, chisq_df, p_value,
                            test_result) {
  york_reg <- matrix(c(intercept, slope, sigma_intercept, sigma_slope),
                     nrow = 2)
  rownames(york_reg) <- c("intercept", "slope")
  colnames(york_reg) <- c("Estimate", "Std_Error")

  weights_matrix <- matrix(c(weights_x, weights_y), ncol = 2)
  colnames(weights_matrix) <- c("weights of x", "weights of y")

  slope_per_iteration <- data.frame("slope" = slope_per_iteration)
  york_arguments <- list("tolerance" = tolerance, "max_iterations" = max_iterations,
                         "mult_samples" = mult_samples, "approx_solution" =
                           approx_solution)
  chisq_test_results <- list("chisq_statistic" = S, "chisq_df" = chisq_df,
                             "p_value" = p_value, "test_result" = test_result)

  if (mult_samples == FALSE) {
    data <- matrix(c(x, y, sd_x, sd_y, r_xy_errors), ncol = 5)
    colnames(data) <- c("x", "y", "sd_x", "sd_y", "r_xy_errors")
  } else {
    data <- list("x" = x_original, "y" = y_original, "sd_x" = sd_x,
                 "sd_y" = sd_y, "r_xy_errors" = r_xy_errors, "x_errors" = x_errors,
                 "y_errors" = y_errors, "mean_x_i" = mean_x_i,
                 "mean_y_i" = mean_y_i)
  }

  output <- list("york_reg" = york_reg,
                 "weights_matrix" = weights_matrix,
                 "slope_per_iteration" = slope_per_iteration,
                 "york_arguments" = york_arguments,
                 "chisq_test_results" = chisq_test_results,
                 "data" = data)

  return(output)
}


