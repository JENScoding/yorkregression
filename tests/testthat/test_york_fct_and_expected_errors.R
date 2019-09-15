test_that("Test york implementation and error messages", {
  ## Input from Table I and Table II in York 1966
  x <- c(0,0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4)
  y <- c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5)

  weights.y = c(1, 1.8, 4, 8, 20, 20, 70, 70, 1e+2, 5e+2)
  weights.x = c(1e+3, 1e+3, 5e+2, 8e+2, 2e+2, 8e+1, 6e+1, 2e+1, 1.8, 1)

  ## Test
  expect_error(york(x, y, weights.y = weights.y , weights.x = weights.x,
                    r.xy = 0), NA)

  expect_error(york(x,y,weights.x = NULL, weights.y = NULL, r.xy = 0.1,
                    sd.x = NULL, sd.y = NULL))
  expect_error(york(x, y, sd.x = NULL, sd.y = NULL, weights.x = NULL,
                    weights.y = NULL, r.xy = 0.2))
  expect_error(york(x = c(1:5), y = c(1:6), sd.x = 4, sd.y = 2,
                    weights.x = NULL, weights.y = NULL, r.xy = 0.2))
  expect_error(york(x , y , sd.x = 4, sd.y = 2,
                    weights.x = weights.x, weights.y = weights.y, r.xy = 0.2))
  x[2:3] <- NA
  expect_warning(york(x, y, weights.y = weights.y , weights.x = weights.x,
                      r.xy = 0))
})
