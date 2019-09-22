test_that("Test implementation in the multi sample case", {
  ## slightly vary x and y and build 5 different samples of each
  x_error <- list()
  y_error <- list()
  x <- list()
  y <- list()
  set.seed(42)
  for (i in 1:12) {
    x_error[[i]] <-  rnorm(10,sd = 0.1)
    y_error[[i]] <-  rnorm(10,sd = 0.05)
    x[[i]] <- c(0.1, 0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4) + x_error[[i]]
    y[[i]] <- c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5) +
      x_error[[i]] + y_error[[i]]
  }

  y <- data.frame(y)
  colnames(y) <- 1:12
  x <- data.frame(x)
  colnames(x) <- 1:12

  ## test

  expect_error(york(x, y, mult_samples = TRUE), NA)

  first <- york(x, y, mult_samples = TRUE)
  expect_type(first$coefficients[2, 1], "double")
  expect_true(first$coefficients[2, 1] < -0.4 &&
                first$coefficients[2, 1] > -0.6)

  expect_error(york(x, y, mult_samples = TRUE, exact_solution = TRUE))

  ## only 5 samples
  x <- x[,-c(1:7)]
  y <- y[,-c(1:7)]

  ## test
  expect_warning(york(x, y, mult_samples = TRUE))

  ## only 4 samples
  x <- x[,-1]
  y <- y[,-1]

  ## test
  expect_error(york(x, y, mult_samples = TRUE))
})

