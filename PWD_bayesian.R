# Load required libraries
library(rstan)

# Define the data
x <- c(0.77, 1.74, 0.81, 1.20, 1.95, 1.20, 0.47, 1.43, 3.37, 2.20, 3, 3.09, 1.51, 2.10, 0.52, 1.62, 1.31, 0.32, 0.59, 0.81, 2.81, 1.87, 1.18, 1.35, 4.75, 2.48, 0.96, 1.89, 0.90, 2.05)

# Specify prior parameters
shape_a <- 1
rate_a <- 1
shape_b <- 1.5
rate_b <- 1
shape_v <- 1
rate_v <- 1

# Define the Stan model code

model_code <- "
data {
  int<lower=0> N;
  real x[N];
  real<lower=0> shape_a;
  real<lower=0> rate_a;
  real<lower=0> shape_b;
  real<lower=0> rate_b;
  real<lower=0> shape_v;
  real<lower=0> rate_v;
}
parameters {
  real<lower=0> a;
  real<lower=1> b;
  real<lower=0> v;
}
model {
  a ~ gamma(shape_a, rate_a);
  b ~ gamma(shape_b, rate_b);
  v ~ gamma(shape_v, rate_v);
  
  for (i in 1:N) {
    target += (b - 1) * log(x[i] / a)
              - 2 * pow(x[i] / a, b) + v * (-1 + exp(-pow(x[i] / a, b)))
              + log(v + exp(-pow(x[i] / a, b)));
  }
}
"


fit <- stan(
  model_code = model_code,
  data = list(
    N = length(x),
    x = x,
    shape_a = shape_a,
    rate_a = rate_a,
    shape_b = shape_b,
    rate_b = rate_b,
    shape_v = shape_v,
    rate_v = rate_v
  ),
  chains = 4,
  iter = 2000,
  warmup = 1000,
  control = list(adapt_delta = 0.95)
)

# Summarize the results
print(fit)