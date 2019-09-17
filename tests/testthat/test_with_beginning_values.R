test_that("Test single values compared to what they were at the beginning", {
  ## Input from Table I and Table II in York 1966
  x <- c(0,0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4)
  y <- c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5)

  weights.y = c(1, 1.8, 4, 8, 20, 20, 70, 70, 1e+2, 5e+2)
  weights.x = c(1e+3, 1e+3, 5e+2, 8e+2, 2e+2, 8e+1, 6e+1, 2e+1, 1.8, 1)

  ## Test
  first <- york(x, y, weights.x = weights.x, weights.y = weights.y, r.xy = 0,
                mult.samples = FALSE, approx.solution = FALSE)
  expect_true(all(round(first$coefficients, 4) == matrix(c(5.4799, -0.4805, 0.0794, 0.0156), nrow = 2)))
  expect_true(all(first$weights == cbind(weights.x, weights.y)))
  expect_true(all(round(first$x.residuals, 4) == c(-0.0002, -0.0003,  0.0008, -0.0018,  0.0185, -0.0380,
                                                    0.0800, -0.2338, -0.0841, 0.8747)))
  expect_true(all(round(first$fitted.y, 4) == c(5.4799, 5.0474, 4.6150, 4.2305, 3.8942, 3.3656, 2.9811,
                                                2.5487, 2.3564, 1.9240)))
  expect_true(round(first$weighted.mean.x, 4) == 4.9110)
  expect_true(round(first$weighted.mean.y, 4) == 3.1200)
  expect_true(all(round(first$reduced.chisq, 4) == 1.4833))
  expect_true(all(round(first$Overall.significance.of.fit$p.value, 4) == 0.1573))
  expect_true(first$number.of.iterations == 4)
  expect_true(round(first$ols.summary$r.squared.ols, 4) == 0.9535)
  expect_true(first$york.arguments[2] == 50)
  expect_true(all(round(first$data[,3], 4) == c(0.0316, 0.0316, 0.0447, 0.0354, 0.0707, 0.1118, 0.1291,
                                                0.2236, 0.7454, 1.0000)))
})
