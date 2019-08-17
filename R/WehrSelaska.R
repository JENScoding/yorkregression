### Algorithm from Wehr&Saleska in The Long-Solved Problem Of The Best-Fit Straigth Line ###

#' Function for york algo
#' @name york
#' @param y y vector
#' @param x x vector
#' @export
york <- function(y, x, tolerance = 1e-10, weights){
  #initial value of b is OLS
  b <- as.numeric(lm(Y~X)[[1]][2])

  b_diff <- 10
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
      alpha[i] <- sqrt(weights[i, 1] * weights[i, 2])
      W[i] <- weights[i, 1] * weights[i, 2] / (b^2 * weights[i, 1] + weights[i, 2])
      X_bar <- X_bar + W[i] * x[i]
      Y_bar <- Y_bar + W[i] * y[i]
      W_sum <- W_sum + W[i]
    }
    X_bar <- X_bar / W_sum
    Y_bar <- Y_bar / W_sum

    Q1 <- 0
    Q2 <- 0
    U <- rep(0,length(x))
    V <- rep(0,length(x))
    beta <- rep(0,length(x))
    for (i in 1:length(x)) {
      U[i] <- x[i] - X_bar
      V[i] <- y[i] - Y_bar
      beta[i] <- W[i] * ((U[i] / weights[i, 1]) + (b * V[i] / weights[i, 2]) - (b * U[i] + V[i]) * 0 / alpha[i])
      Q1 <- Q1 + W[i] * beta[i] * V[i]
      Q2 <- Q2 + W[i] * beta[i] * U[i]
    }
   #Q1 <- sum(W*beta*V)
    #Q2 <- sum(W*beta*U)
    b <- round(Q1 / Q2, 100)
    b_diff <- abs(b - b_old)
    count <- count + 1

  }

  a <- Y_bar - b * X_bar

  x_mean <- 0
  x_adj <- rep(0,length(X))

  for(i in 1:length(X)){
    x_adj[i] <- X_bar + beta[i]
    x_mean <- x_mean + W[i] * beta[i]
  }

  x_mean <- x_mean / (W_sum * (length(X) - 2))

  sigma_b <- 0
  chisq_w <- 0
  u <- rep(0,length(X))

  for(i in 1:length(X)){
    u[i] <- x_adj[i] - x_mean
    sigma_b <- sigma_b + W[i]*u[i]^2
    chisq_w <- chisq_w + W[i]*(Y[i]-b*X[i]-a)^2
  }

  sigma_b <- sqrt(1 / sigma_b)
  sigma_a <- sqrt(x_mean^2 * sigma_b^2 + 1 / W_sum)
  chisq_w <- chisq_w / (length(X) - 2)
  sigma_x <- sqrt(2 / (length(X) - 2))

  residuals <- y - (a + b * x)

  mt <- matrix(c(a, b, sigma_a, sigma_b), nrow = 2)
  rownames(mt) <- c("intercept", "slope")
  colnames(mt) <- c("Estimate", "Std.Error")
  est <- list("coefficients" = mt, "weighting.vector" = W, "residuals" = residuals, "mean.x" = X_bar, "mean.y" = Y_bar ,"Std.Error.chi" = sigma_x, "iterations" = count, U, V, x, y, x_mean, u)
  return(est)
}


## end of (relevant) script
