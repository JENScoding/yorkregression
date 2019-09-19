### Test internal functions of the internal_functions.R file

test_that("Test f_rewrite function", {
  ## Input from Table I and Table II in York 1966
  x <- c(0,0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4)
  y <- c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5)
  x2 <- x
  x2[4] <- NA
  y2 <- y
  y2[2] <- NA

  weights_x <- NULL
  weights_y <- NULL
  sd_x <- 0.4
  sd_y <- 0.8
  r_xy_errors <- 0

  ## Test
  expect_error(f_rewrite(x, y, weights_x, weights_y,
                         sd_x, sd_y, r_xy_errors), NA)
  expect_warning(f_rewrite(x2, y2, weights_x, weights_y,
                         sd_x, sd_y, r_xy_errors))
  suppressWarnings(
    expect_equal(
      length(f_rewrite(x2, y2, weights_x, weights_y,
                       sd_x, sd_y, r_xy_errors)$x), 8))
  suppressWarnings(
    expect_equal(
      length(f_rewrite(x2, y2, weights_x, weights_y,
                       sd_x, sd_y, r_xy_errors)$x),
      length(f_rewrite(x2, y2, weights_x, weights_y,
                       sd_x, sd_y, r_xy_errors)$y)))
})


test_that("Test f_rewrite_mult function", {
  ## Input
  x <- data.frame(c(2,3,5), c(4,2,1))
  y <- data.frame(c(2,2,5), c(4,5,3))
  x2 <- x
  x2[1, 1] <- NA
  y2 <- y
  y2[3, 2] <- NA

  ## Test
  expect_error(f_rewrite_mult(x, y), NA)
  expect_warning(f_rewrite_mult(x2, y2))
  suppressWarnings(
    expect_equal(
      length(f_rewrite_mult(x2, y2)$x),
      length(f_rewrite_mult(x2, y2)$y)))
})


test_that("Test variance and correlation functions", {
  ## Input from Table I and Table II in York 1966
  x <- c(0,0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4)
  y <- c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5)
  x2 <- data.frame(t(x))
  y2 <- data.frame(t(y))

  x3 <- c(0,0.9, 1.8, 2.6, 5.3, 14.4, 5.2, 16.1, 6.5, 7.4)
  y3 <- c(5.9, 5.4, 4.4, 9.6, 3.5, 3.7, 12.8, 2.8, 2.4, 1.5)
  x4 <- data.frame(rbind(t(x), t(x3)))
  y4 <- data.frame(rbind(t(y), t(y3)))

  ## Test
  expect_error(f_var_row(x2), NA)
  library(stats)
  expect_equal(f_var_row(x2),
               var(x))
  expect_error(f_corr_row(x2, y2), NA)
  expect_equal(f_corr_row(x2, y2),
               cor(x, y))
  expect_equal(f_corr_row(x4, y4),
               c(cor(x, y), cor(x3, y3)))
})


test_that("Test OLS implementation", {
  ## Input from Table I and Table II in York 1966
  x <- c(0,0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4)
  y <- c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5)

  weights_y <- c(1, 1.8, 4, 8, 20, 20, 70, 70, 1e+2, 5e+2)
  weights_x <- c(1e+3, 1e+3, 5e+2, 8e+2, 2e+2, 8e+1, 6e+1, 2e+1, 1.8, 1)

  ## Test
  first <- york(x, y, weights_x = weights_x, weights_y = weights_y,
                r_xy_errors = 0)
  expect_equal(first$ols_summary$coefficients_ols[1,2],
               summary(lm(y ~ x))$coefficients[1, 2])
  expect_equal(first$ols_summary$coefficients_ols[2,2],
               summary(lm(y ~ x))$coefficients[2, 2])
  expect_equal(first$ols_summary$coefficients_ols[2,2],
               summary(lm(y ~ x))$coefficients[2, 2])
  expect_true(first$ols_summary$total_sum_of_squares >
                first$ols_summary$residual_sum_of_squares)
  expect_equal(first$ols_summary$r_squared_ols,
               summary(lm(y ~ x))$r.squared)
  expect_equal(first$ols_summary$r_squared_ols,cor(x,y, method = "pearson")^2)
  expect_equal(sum(first$ols_summary$residuals_ols),
               0)
  expect_equal(mean(first$ols_summary$fitted_y_ols),
               mean(y))
  expect_equal(sum(x * first$ols_summary$residuals_ols), 0)
  expect_equal(sum(first$ols_summary$residuals_ols *
                     first$ols_summary$fitted_y_ols), 0)
  expect_equal(first$ols_summary$coefficients_ols[1,1] /
                 first$ols_summary$coefficients_ols[1,2],
               summary(lm(y ~ x))$coefficients[1, 3])
  expect_equal(first$ols_summary$coefficients_ols[2,1] /
                 first$ols_summary$coefficients_ols[2,2],
               summary(lm(y ~ x))$coefficients[2, 3])
  expect_equal(first$ols_summary$f_statistic_ols,
               as.numeric(summary(lm(y ~ x))$fstatistic[1]))
  expect_equal(first$ols_summary$r_squared_adjusted_ols,
               summary(lm(y ~ x))$adj.r.squared)
})


test_that("Test f_cubic_root function", {
  ## Input
  x <- c(0,0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4)
  y <- c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5)
  weights_x <- 0.4
  weights_y <- 0.8
  no_corr <- 0
  corr <- 0.1
  slope <- -0.5
  se <- 0.03

  ## Test
  expect_error(f_cubic_root(x, y, weights_x, weights_y,
                            no_corr, slope, se), NA)
  expect_error(f_cubic_root(x, y, weights_x, weights_y,
                              corr, slope, se))
  expect_equal(length(f_cubic_root(x, y, weights_x, weights_y,
                                   no_corr, slope, se)), 7)
})
