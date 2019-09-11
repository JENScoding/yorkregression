test_that("Test implementation of york.predicts function", {
  ## Input from Table I and Table II in York 1966
  x <- c(0,0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4)
  y <- c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5)

  weights.y = c(1, 1.8, 4, 8, 20, 20, 70, 70, 1e+2, 5e+2)
  weights.x = c(1e+3, 1e+3, 5e+2, 8e+2, 2e+2, 8e+1, 6e+1, 2e+1, 1.8, 1)

  ## Test
  first <- york(x, y, weights.x = weights.x, weights.y = weights.y, r.xy = 0, mult.samples = F, exact.solution = T)
  new <- c(2,3,4,3.6)
  expect_error(york.predict(first, new), NA)

  second <- york.predict(first, new)
  expect_true(all(round(second$prediction[,2], 4) == c(4.5080, 4.0308, 3.5535, 3.7444)))
  expect_true(all(colnames(second$prediction) == c("x", "predicted.y")))

  new <- data.frame(c(2,3,4,3.6), fix.empty.names = F)
  second <- york.predict(first, new)
  expect_true(all(round(second$prediction[,2], 4) == c(4.5080, 4.0308, 3.5535, 3.7444)))
  expect_true(all(colnames(second$prediction) == c("x", "predicted.y")))
})




