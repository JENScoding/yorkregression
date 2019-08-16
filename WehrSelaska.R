### Algorithm from Wehr&Saleska in The Long-Solved Problem Of The Best-Fit Straigth Line ###


# Input from Table I and Table II in York 1966

X <- c(0,0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4)
Y <- c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5)

w_X <- c(1e+3, 1e+3, 5e+2, 8e+2, 2e+2, 8e+1, 6e+1, 2e+1, 1.8, 1)
w_Y <- c(1, 1.8, 4,8, 20, 20, 70, 70, 1e+2, 5e+2)





#initial value of b and tolerance level
tolerance <- 1e-6
b_ols <- function(y, x) { 
  fit <- lm(y ~ x)
  return(as.numeric(fit$coefficients[2]))
}
b <- b_ols(Y, X)
b_diff <- 1


#function for single variance
#single_var <- function(x, element){
 # n <- length(x)
  #mu <- mean(x)
  #return(1 / (n -1) * (mu - x[element])**2)
#}

#single_var(X,3)

count <- 0

while (b_diff > tolerance) {
  b_old <- b
  X_bar <- 0
  Y_bar <- 0
  W_sum <- 0
  alpha <- rep(0,length(X))
  W <- rep(0,length(X))
  for (i in 1:length(X)) {
    # omega, Jonas hat aber gesagt: "Lass es weg"
    alpha[i] <- sqrt(w_X[i] * w_Y[i])
    W[i] <- alpha[i]^2 / (b^2 * w_Y[i] + w_X[i])
    X_bar <- X_bar + W[i] * X[i]
    Y_bar <- Y_bar + W[i] * Y[i]
    W_sum <- W_sum + W[i]
  }
  X_bar <- X_bar / W_sum 
  Y_bar <- Y_bar / W_sum 
  
  Q1 <- 0
  Q2 <- 0
  U <- rep(0,length(X))
  V <- rep(0,length(X))
  beta <- rep(0,length(X))
  for (i in 1:length(X)) {
    U[i] <- X[i] - X_bar
    V[i] <- Y[i] - Y_bar
    beta[i] <- W[i] * ((U[i] / w_Y[i]) + (b * V[i] / w_X[i]) - (b * U[i] + V[i]) * 0 / alpha[i])
    Q1 <- Q1 + W[i] * beta[i] * V[i]
    Q2 <- Q2 + W[i] * beta[i] * U[i]
  }
 #Q1 <- sum(W*beta*V)
  #Q2 <- sum(W*beta*U)
  b <- Q1 / Q2
  b_diff <- abs(b - b_old)
  count <- count + 1
  print(b)
  
}

a <- Y_bar - b * X_bar


x_mean <- 0



x <- rep(0,length(X))

for(i in 1:length(X)){
  x[i] <- X_bar + beta[i]
  x_mean <- x_mean + W[i] * beta[i]
}

x_mean <- x_mean / (W_sum * (length(X) - 2))

sigma_b <- 0
chisq_w <- 0
u <- rep(0,length(X))

for(i in 1:length(X)){
  u[i] <- x[i] - x_mean
  sigma_b <- sigma_b + W[i]*u[i]^2
  chisq_w <- chisq_w + W[i]*(Y[i]-b*X[i]-a)^2
}

sigma_b <- sqrt(1 / sigma_b)
sigma_a <- sqrt(x_mean^2 * sigma_b^2 + 1 / W_sum)
chisq_w <- chisq_w / (length(X) - 2)
sigma_x <- sqrt(2 / (length(X) - 2))




## end of script


