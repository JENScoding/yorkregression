### Algorithm from Wehr&Saleska in The Long-Solved Problem Of The Best-Fit Straigth Line ###


## Input from Table I and Table II in York 1966

X <- c(0,0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4)
Y <- c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5)

w_X <- c(1e+3, 1e+3, 5e+2, 8e+2, 2e+2, 8e+1, 6e+1, 2e+1, 1.8, 1)
w_Y <- c(1, 1.8, 4,8, 20, 20, 70, 70, 1e+2, 5e+2)

my_w <- data.frame(
  "Y" = c(1, 1.8, 4,8, 20, 20, 70, 70, 1e+2, 5e+2),
  "X" = c(1e+3, 1e+3, 5e+2, 8e+2, 2e+2, 8e+1, 6e+1, 2e+1, 1.8, 1))

## function for algo

york <- function(y, x, tolerance = 1e-6){
  #initial value of b and tolerance level
  b_ols <- function(y, x) { 
    fit <- lm(y ~ x)
    return(as.numeric(fit$coefficients[2]))
  }
  b <- b_ols(Y, X)
  
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
   # print(b)
    
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
  
  coef <- c(a,b)
  coefn <- attr(coef, "coefficients")
  est <- list("coefficients" = coefn, W, X_bar, Y_bar, sigma_a, sigma_b, sigma_x, count, U, V, x, x_mean, u)
  return(est)
}

first <- york(Y, g)

first$coefficients[1]

## end of (relevant) script


## Testing 

sec <- lm(Y~X)
lm
function (formula, data, subset, weights, na.action, method = "qr", 
          model = TRUE, x = FALSE, y = FALSE, qr = TRUE, singular.ok = TRUE, 
          contrasts = NULL, offset, ...) 
{
  ret.x <- x
  ret.y <- y
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", 
               "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  if (method == "model.frame") 
    return(mf)
  else if (method != "qr") 
    warning(gettextf("method = '%s' is not supported. Using 'qr'", 
                     method), domain = NA)
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  w <- as.vector(model.weights(mf))
  if (!is.null(w) && !is.numeric(w)) 
    stop("'weights' must be a numeric vector")
  offset <- as.vector(model.offset(mf))
  if (!is.null(offset)) {
    if (length(offset) != NROW(y)) 
      stop(gettextf("number of offsets is %d, should equal %d (number of observations)", 
                    length(offset), NROW(y)), domain = NA)
  }
  if (is.empty.model(mt)) {
    x <- NULL
    z <- list(coefficients = if (is.matrix(y)) matrix(, 0, 
                                                      3) else numeric(), residuals = y, fitted.values = 0 * 
                y, weights = w, rank = 0L, df.residual = if (!is.null(w)) sum(w != 
                                                                                0) else if (is.matrix(y)) nrow(y) else length(y))
    if (!is.null(offset)) {
      z$fitted.values <- offset
      z$residuals <- y - offset
    }
  }
  else {
    x <- model.matrix(mt, mf, contrasts)
    z <- if (is.null(w)) 
      lm.fit(x, y, offset = offset, singular.ok = singular.ok, 
             ...)
    else lm.wfit(x, y, w, offset = offset, singular.ok = singular.ok, 
                 ...)
  }
  class(z) <- c(if (is.matrix(y)) "mlm", "lm")
  z$na.action <- attr(mf, "na.action")
  z$offset <- offset
  z$contrasts <- attr(x, "contrasts")
  z$xlevels <- .getXlevels(mt, mf)
  z$call <- cl
  z$terms <- mt
  if (model) 
    z$model <- mf
  if (ret.x) 
    z$x <- x
  if (ret.y) 
    z$y <- y
  if (!qr) 
    z$qr <- NULL
  z
}
<bytecode: 0x563be91068e0>
  <environment: namespace:stats>
