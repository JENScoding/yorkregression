## slightly vary x and y
vary <- list()
x <- list()
y <- list()
set.seed(42)
for (i in 1:5){
  vary[[i]] <- sample(seq(0.1,2,0.01), 10, replace = T)
  x[[i]] <- c(0, 0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4) * vary[[i]]
  y[[i]] <- c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5) * vary[[i]]
}

y <- data.frame(y)
colnames(y) <- 1:5
x <- data.frame(x)
colnames(x) <- 1:5

rm(vary)

weights.y = c(1, 1.8, 4, 8, 20, 20, 70, 70, 1e+2, 5e+2)
weights.x = c(1e+3, 1e+3, 5e+2, 8e+2, 2e+2, 8e+1, 6e+1, 2e+1, 1.8, 1)