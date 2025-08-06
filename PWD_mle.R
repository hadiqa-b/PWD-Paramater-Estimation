library(maxLik)
install.packages("KScorrect")
library(KScorrect)

x <- c(0.77, 1.74, 0.81, 1.20, 1.95, 1.20, 0.47, 1.43, 3.37, 2.20, 3, 3.09, 1.51, 2.10, 0.52, 1.62, 1.31, 0.32, 0.59, 0.81, 2.81, 1.87, 1.18, 1.35, 4.75, 2.48, 0.96, 1.89, 0.90, 2.05)
n = length(x)

loglik <- function(params, x) {
  a <- params[1]
  b <- params[2]
  v <- params[3]
  
  
  # Compute the terms of the PDF
  term1 <- (b - 1) * log(x / a)
  term2 <- -2 * (x / a)^b + v * (-1 + exp(-(x / a)^b))
  term3 <- log(v + exp((x / a)^b))
  
  # Calculate the log likelihood
  ll <- (n*log(b / a)) + sum(term1 + term2 + term3)
  
  # Return the negative log likelihood since maxLik minimizes
  -ll
}


pdf <- function(params, x) {
  
  a <- params[1]
  b <- params[2]
  v <- params[3]
  
  
  # Calculate the terms of the PDF
  term1 <- (b / a)
  term2 <- (x / a)^(b - 1)
  term3 <- exp(-2 * (x / a)^b + v * (-1 + exp(-(x / a)^b)))
  term4 <- v + exp((x / a)^b)
  
  # Return the product of all terms
  term1 * term2 * term3 * term4
}

cdf <- function(params, x) {
  
  a <- params[1]
  b <- params[2]
  v <- params[3]
  
  
  # Compute the inner term (x/a)^b
  inner_term <- (x / a)^b
  
  # Calculate the full expression inside the exp function
  exp_term <- -(inner_term) + v * (exp(-inner_term) - 1)
  
  1 - exp(exp_term)
  
}

initial_values <- c(a = 1.7, b = 1.5, v = 0.5)
lower_bounds <- c(0, 1, 0)
upper_bounds <- c(Inf, Inf, Inf)

fit <- optim(initial_values, loglik, method = "L-BFGS-B", lower = lower_bounds, upper = upper_bounds, x = x, hessian = TRUE)
estimated_params <- fit$par
hessian_matrix <- fit$hessian
cov_matrix <- solve(hessian_matrix)
cov_matrix

standard_errors <- sqrt(diag(cov_matrix))
standard_errors

estimated_params
logLik_value <- -fit$value


x_jittered <- jitter(x, amount = 1e-5)

cdf_with_params <- function(q) {
  cdf(fit$par, q)
}

ks_test <- ks.test(x_jittered, cdf_with_params)

print(ks_test)

k <- length(estimated_params)  
aic <- -2 * logLik_value + 2 * k
bic <- -2 * logLik_value + log(n) * k
hqic <- -2 * logLik_value + 2 * k * log(log(n))


list(estimated_params = estimated_params, standard_errors = standard_errors, ks_test = ks_test, AIC = aic, BIC = bic, HQIC = hqic)
