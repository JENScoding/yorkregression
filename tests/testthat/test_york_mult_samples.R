load("R/multiple_samples.RData")

test_that("Test that our york implementation is correct", {

  expect_error(york(x, y, mult.samples = T), NA)

  first <- york(x, y, mult.samples = T)
  expect_type(first$coefficients.york[2, 1], "double")
  expect_true(first$coefficients.york[2, 1] < -0.47 && first$coefficients.york[2, 1] > -0.55)
})

