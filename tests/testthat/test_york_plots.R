test_that("Test york.plots function", {
  ## Input from Table I and Table II in York 1966
  x <- c(0,0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4)
  y <- c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5)

  weights_y = c(1, 1.8, 4, 8, 20, 20, 70, 70, 1e+2, 5e+2)
  weights_x = c(1e+3, 1e+3, 5e+2, 8e+2, 2e+2, 8e+1, 6e+1, 2e+1, 1.8, 1)

  ## Test
  first <- york(x, y, weights_x = weights_x, weights_y = weights_y, r_xy_errors = 0,
                mult_samples = FALSE, approx_solution = TRUE)
  expect_error(plot(first), NA)
  expect_error(plot.york(first), NA)
  second <- data.frame("x" = rnorm(100), "y" = rnorm(100))
  expect_error(plot.york(second))
})




