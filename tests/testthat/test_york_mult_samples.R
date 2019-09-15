test_that("Test implementation in the multi sample case", {
  ## slightly vary x and y and build 5 different samples of each
  x.error <- list()
  y.error <- list()
  x <- list()
  y <- list()
  set.seed(42)
  for (i in 1:5){
    x.error[[i]] <-  rnorm(10,sd = 0.1)
    y.error[[i]] <-  rnorm(10,sd = 0.05)
    x[[i]] <- c(0.1, 0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4) + x.error[[i]]
    y[[i]] <- c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5) + x.error[[i]] + y.error[[i]]
  }

  y <- data.frame(y)
  colnames(y) <- 1:5
  x <- data.frame(x)
  colnames(x) <- 1:5

  ## test

  expect_error(york(x, y, mult.samples = T), NA)

  first <- york(x, y, mult.samples = T)
  expect_type(first$coefficients[2, 1], "double")
  expect_true(first$coefficients[2, 1] < -0.47 && first$coefficients[2, 1] > -0.55)

  expect_error(york(x, y, mult.samples = T, exact.solution = T))
})

