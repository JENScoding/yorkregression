test_that("Test internal functions", {
  ## Input from Table I and Table II in York 1966
  x <- c(0,0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4)
  y <- c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5)
  x2 <- data.frame(t(x))
  y2 <- data.frame(t(y))

  ## Test
  expect_error(calc.var(x2), NA)
  library(stats)
  expect_equal(calc.var(x2),
               var(x))
  expect_error(calc.corr(x2, y2), NA)
  expect_equal(calc.corr(x2, y2),
               cor(x, y))

  expect_equal(first$ols.summary$r.squared.adjusted.ols,
               summary(lm(y ~ x))$adj.r.squared)
})

