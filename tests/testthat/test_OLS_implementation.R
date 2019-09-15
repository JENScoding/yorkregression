test_that("Test OLS implementation", {
  ## Input from Table I and Table II in York 1966
  x <- c(0,0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4)
  y <- c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5)

  weights.y = c(1, 1.8, 4, 8, 20, 20, 70, 70, 1e+2, 5e+2)
  weights.x = c(1e+3, 1e+3, 5e+2, 8e+2, 2e+2, 8e+1, 6e+1, 2e+1, 1.8, 1)

  ## Test
  first <- york(x, y, weights.x = weights.x, weights.y = weights.y,
                r.xy = 0)
  expect_equal(first$ols.summary$coefficients.ols[1,2],
               summary(lm(y ~ x))$coefficients[1, 2])
  expect_equal(first$ols.summary$coefficients.ols[2,2],
               summary(lm(y ~ x))$coefficients[2, 2])
  expect_equal(first$ols.summary$coefficients.ols[2,2],
               summary(lm(y ~ x))$coefficients[2, 2])
  expect_true(first$ols.summary$total.sum.of.squares > first$ols.summary$residual.sum.of.squares)
  expect_equal(first$ols.summary$r.squared.ols,
               summary(lm(y ~ x))$r.squared)
  expect_equal(first$ols.summary$r.squared.ols,cor(x,y, method = "pearson")^2)
  expect_equal(sum(first$ols.summary$residuals.ols),
               0)
  expect_equal(mean(first$ols.summary$fitted.y.ols),
               mean(y))
  expect_equal(sum(x*first$ols.summary$residuals.ols), 0)
  expect_equal(sum(first$ols.summary$residuals.ols*first$ols.summary$fitted.y.ols), 0)
  expect_equal(first$ols.summary$coefficients.ols[1,1] / first$ols.summary$coefficients.ols[1,2],
               summary(lm(y ~ x))$coefficients[1, 3])
  expect_equal(first$ols.summary$coefficients.ols[2,1] / first$ols.summary$coefficients.ols[2,2],
               summary(lm(y ~ x))$coefficients[2, 3])
  expect_equal(first$ols.summary$f.statistic.ols,
               as.numeric(summary(lm(y ~ x))$fstatistic[1]))
  expect_equal(first$ols.summary$r.squared.adjusted.ols,
               summary(lm(y ~ x))$adj.r.squared)
})

