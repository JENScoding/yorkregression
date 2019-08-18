test_that("Test that our york implementation is correct", {
  x <- c(0,0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4)
  y <- c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5)

  weights_y = c(1, 1.8, 4, 8, 20, 20, 70, 70, 1e+2, 5e+2)
  weights_x = c(1e+3, 1e+3, 5e+2, 8e+2, 2e+2, 8e+1, 6e+1, 2e+1, 1.8, 1)

  expect_error(york(y, x, weights.y = weights_y, weights.x = weights_x), NA)

  first <- york(y, x, weights.y = weights_y, weights.x = weights_x)
  expect_type(first$coefficients[2, 1], "double")
#  expect_equal(first$coefficients[2, 1], -0.481)
  expect_true(first$coefficients[2, 1] < -0.47 && first$coefficients[2, 1] > -0.55)
})


## Input from Table I and Table II in York 1966


