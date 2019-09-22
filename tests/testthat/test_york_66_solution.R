test_that("Test implementation of exact solution", {
  ## Input from Table I and Table II in York 1966
  x <- c(0,0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4)
  y <- c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5)

  weights_y = c(1, 1.8, 4, 8, 20, 20, 70, 70, 1e+2, 5e+2)
  weights_x = c(1e+3, 1e+3, 5e+2, 8e+2, 2e+2, 8e+1, 6e+1, 2e+1, 1.8, 1)

  ## Test
  expect_error(york(x, y, weights_x = weights_x, weights_y = weights_y,
                    r_xy_errors = 0, approx_solution = TRUE), NA)

  # if weights for x and y are 1 we have orthogonal regression
  first <- york(x, y, weights_x = 1, weights_y = 1, r_xy_errors = 0)
  expect_equal(round(first$coefficients[2,1], 3),
               -0.546)

  # if weights x goes to inf and weights y are 1 we have ols
  second <- york(x, y, weights_x = 1e10, weights_y = 1, r_xy_errors = 0)
  expect_equal(round(second$coefficients[2,1], 3),
               round(second$ols_summary$coefficients_ols[2,1], 3))

  # compare approx value with the one in York (1966)
  third <- york(x, y, weights_x = weights_x, weights_y = weights_y,
                r_xy_errors = 0, approx_solution = TRUE)
  expect_type(third$coefficients[2, 1], "double")
  expect_true(third$coefficients[2, 1] < -0.477 && third$coefficients[2, 1] >
                -0.478)
})





