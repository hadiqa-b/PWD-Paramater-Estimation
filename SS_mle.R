

library(maxLik)
library(ggplot2)

x = c( 0.097, 0.014, 0.03, 0.134, 0.240, 0.084, 0.146, 0.024, 0.045, 0.004, 0.099, 0.277, 0.472, 0.094, 0.023, 0.146, 0.030, 0.031, 0.104, 0.105, 0.036, 0.065, 0.022, 0.098, 0.178, 0.059, 0.014, 0.007, 0.007, 0.286)
y = c( 0.084, 0.236, 0.315, 0.199, 0.252, 0.103, 0.455, 0.135, 0.348, 0.321, 0.166, 0.04, 0.027, 0.519, 0.017, 0.821, 0.942, 0.27, 0.008, 0.03, 0.177, 0.268, 0.18, 0.796, 0.245, 0.703, 0.045, 0.314, 0.281, 0.652)



n = length(x)
m = length(y)

loglik <- function(params, x, y) {
  a1 <- params[1]
  b1 <- params[2]
  v1 <- params[3]
  a2 <- params[4]
  b2 <- params[5]
  v2 <- params[6]
  
  # Compute the terms of the PDF
  term1 <-  (b1 - 1) * log(x / a1)
  term2 <- -2 * (x / a1)^b1 + v1 * (-1 + exp(-(x / a1)^b1))
  term3 <- log(v1 + exp((x / a1)^b1))
  
  term4 <- (b2 - 1) * log(y / a2)
  term5 <- -2 * (y / a2)^b2 + v2 * (-1 + exp(-(y / a2)^b2))
  term6 <- log(v2 + exp((y / a2)^b2))
  
  
  # Calculate the log likelihood
  ll <- (n*log(b1 / a1)) +  (m*log(b2 / a2)) + sum(term1 + term2 + term3 + term4 + term5 + term6)
  
  # Return the negative log likelihood since maxLik minimizes
  -ll
}


initial_values <- c(a1 = 1.7, b1 = 1.5, v1 = 0.5, a2 = 1.7, b2 = 1.5, v2 = 0.5)

lower_bounds <- c(0.01, 1, 0.01, 0.01, 1, 0.01)
upper_bounds <- c(Inf, Inf, Inf, Inf, Inf, Inf)

fit <- optim(initial_values, loglik, method = "L-BFGS-B", lower = lower_bounds, upper = upper_bounds, x = x, y=y, hessian = TRUE)
estimated_params <- fit$par
hessian_matrix <- fit$hessian
cov_matrix <- solve(hessian_matrix)
cov_matrix

standard_errors <- sqrt(diag(cov_matrix))
standard_errors

estimated_params

