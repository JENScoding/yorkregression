### Examples

library(York)

load("R/original_data.RData") # York 66

# 1.
(york.output <- york(x, y)) # no weights are specified

## work on algo file:
#  -- under if(exact.solution == T) {
# stop("There is no exact solution in case of multiple samples!")
# } --
# if(is.data.frame(x) & is.data.frame(y) == F) {
#   stop("Inputs x and y must be of class data.frame!")
# }
# if(ncol(x) | ncol(y) == 1) {
#   stop("You need more than one sample if you specify mult.samples = T")
# }
# if(ncol(x) != ncol(y) & nrow(x) != nrow(y)) {
#   stop("x and y must have the same number of columns/ rows")
# }

# 2.
york.output <- york(x, y, weights.x = weights.x, weights.y = weights.y, r.xy = 0)
york.output$data
york.output$coefficients.york

york.output <- york(x, y, weights.x = weights.x, weights.y = weights.y, r.xy = 0, exact.solution = T)
york.output$data
york.output$coefficients.york

york.output <- york(x, y, weights.x = weights.x, weights.y = weights.y, r.xy = 0.9) # Weight in
york.output$data
york.output$coefficients.york

# 3.
york.output <- york(x, y, sd.x = 2, sd.y = 5, r.xy = 0.8)
york.output$data
york.output$coefficients.york

# 4.
(york.output <- york(x, y, weights.x = weights.x, weights.y = weights.y, r.xy = 0, mult.samples = F))

class(york.output)

york.plots(york.output = york.output)

(york.output <- york(x, y, weights.x = weights.x, weights.y = weights.y, r.xy = 0, mult.samples = F, exact.solution = T))

york.plots(york.output = york.output)


load("R/multiple_samples.RData") # multiple samples

(york.output <- york(x, y, mult.samples = T))

york.plots(york.output = york.output)


load("R/original_data_with_NA.RData") # Show NA case

View(x)
(york.output <- york(x, y, weights.x = weights.x, weights.y = weights.y, r.xy = 0))
