## slightly vary x and y
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

rm(x.error, y.error, i)
