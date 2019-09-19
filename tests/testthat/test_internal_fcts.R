test_that("Test internal functions", {
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

