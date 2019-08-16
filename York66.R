
### Least- Squares Fitting Of A Straight Line; Derek York, 1966 ###

# Input from Table I and Table II

X <- c(0,0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4)

Y <- c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5)

fit_lm <- lm(Y ~ X)

w_X <- c(1e+3, 1e+3, 5e+2, 8e+2, 2e+2, 8e+1, 6e+1, 2e+1, 1.8, 1)
w_Y <- c(1, 1.8, 4,8, 20, 20, 70, 70, 1e+2, 5e+2)


## Apply formula and use lm estimate as intitial value for b 

b_init <- fit_lm$coefficients[2]
b <- b_init
W <- w_X * w_Y / (b**2 * w_Y + w_X) # see formula 14
  # see formula 19 and 20 for the following:
X_bar <- sum(W * X) / sum(W) 
Y_bar <- sum(W * Y) / sum(W)
U <- X - X_bar 
V <- Y - Y_bar

alpha <- 2 * sum(U * V * W**2 / w_X) / (3 * sum(U**2 * W**2 / w_X)) # see page 1084
beta <- (sum(V**2 * W**2 / w_X) - sum(W * U**2)) / (3 * sum(U**2 * W**2 / w_X)) # see page 1084
gamma <- - sum(U * V * W) / (sum(U**2 * W**2 / w_X)) # see page 1084
 

# try with polyroot function to find root of cubic equation
b_cubic <- polyroot(c(-gamma, 3 * beta, -3 * alpha ,1))

b_solution <- as.numeric(b_cubic[1])
b_solution

b1 <- b_solution
b0 <-Y_bar - b1 * X_bar # see equtation 17
fitted_Y <- b0 + b1*X

plot(X,Y)
lines(X, fitted_Y, col = "red")
lines(X, fit_lm$fitted.values, col = "blue", lty = "dashed")
legend("topright",legend = c("OLS","York"), fill = c("blue","red"))


# use formula of York to find b, given on page 1084
phi <- acos((alpha**3 - 3 /2 * alpha * beta + 0.5 * gamma) / (alpha**2 - beta)**(3 / 2))

b_cubic2 <- alpha + 2 * (alpha**2 - beta)**0.5 * cos( 1 / 3 *(phi + 2 * pi * c(0:2)))
b_cubic2
b_york <- b_cubic2[3]
b_york # same as in paper

b1 <- b_york
b0 <-Y_bar - b1 * X_bar # see equtation 17
fitted_Y <- b0 + b1*X

plot(X,Y)
lines(X, fitted_Y, col = "red")
lines(X, fit_lm$fitted.values, col = "blue", lty = "dashed")
legend("topright",legend = c("OLS","York"), fill = c("blue","red"))

#error variance of the slope
(var_b1 <- (1 / (length(X)-2)) * (sum(W * (b_york * U - V)**2)) / (sum(W * U**2)))
sd_b1 <- sqrt(var_b1)

#error variance of the intercept
(var_b0 <- var_b1 * (sum(W * X**2)) / (sum(W)))
sd_b0 <- sqrt(var_b0)


## end of script

