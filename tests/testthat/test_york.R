library(testthat)
setwd("")
load("original_data.RData")
test_that("Test that our york implementation is correct", {
    expect_error(york(x, y, weights.y = weights.y , weights.x = weights.x,
                      rxy = 0), NA)
        first <- york(x, y, weights.x = weights.x, weights.y = weights.y,
                      rxy = 0)
    #expect_true(first$coefficients[2, 1] < 1 && first$coefficients[2, 1] > -1)
    expect_type(first$coefficients.york[1, 1], "double")
    expect_equal(first$coefficients.ols[1,2], summary(lm(y ~ x))$coefficients[1, 2])
    expect_equal(first$coefficients.ols[2,2], summary(lm(y ~ x))$coefficients[2, 2])
    expect_error(york(x,y,weights.x = NULL, weights.y = NULL, rxy = 0.1,
                      sd.x = NULL, sd.y = NULL))
    expect_error(york(x, y, sd.x = NULL, sd.y = NULL, weights.x = NULL,
                          weights.y = NULL, rxy = 0.2))
    expect_error(york(x = c(1:5), y = c(1:6), sd.x = 4, sd.y = 2,
                      weights.x = NULL, weights.y = NULL, rxy = 0.2))
    expect_error(york(x , y , sd.x = 4, sd.y = 2,
                      weights.x = weights.x, weights.y = weights.y, rxy = 0.2))

})


## Input from Table I and Table II in York 1966


