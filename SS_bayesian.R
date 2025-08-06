library(rstan)
library(bayesplot)
library(ggplot2)
library(dplyr)
library(zoo)

x = c( 0.097, 0.014, 0.03, 0.134, 0.240, 0.084, 0.146, 0.024, 0.045, 0.004, 0.099, 0.277, 0.472, 0.094, 0.023, 0.146, 0.030, 0.031, 0.104, 0.105, 0.036, 0.065, 0.022, 0.098, 0.178, 0.059, 0.014, 0.007, 0.007, 0.286)
y = c( 0.084, 0.236, 0.315, 0.199, 0.252, 0.103, 0.455, 0.135, 0.348, 0.321, 0.166, 0.04, 0.027, 0.519, 0.017, 0.821, 0.942, 0.27, 0.008, 0.03, 0.177, 0.268, 0.18, 0.796, 0.245, 0.703, 0.045, 0.314, 0.281, 0.652)


model_code <- "
data {
  int<lower=0> N;
  real x[N];
  int<lower=0> M;
  real y[M];
  real<lower=0> shape_a1;
  real<lower=0> rate_a1;
  real<lower=0> shape_b1;
  real<lower=0> rate_b1;
  real<lower=0> shape_v1;
  real<lower=0> rate_v1;
  real<lower=0> shape_a2;
  real<lower=0> rate_a2;
  real<lower=0> shape_b2;
  real<lower=0> rate_b2;
  real<lower=0> shape_v2;
  real<lower=0> rate_v2;
}
parameters {
  real<lower=0> a1;
  real<lower=1> b1;
  real<lower=0> v1;
  real<lower=0> a2;
  real<lower=1> b2;
  real<lower=0> v2;
}
model {
  a1 ~ gamma(shape_a1, rate_a1);
  b1 ~ gamma(shape_b1, rate_b1);
  v1 ~ gamma(shape_v1, rate_v1);
  a2 ~ gamma(shape_a2, rate_a2);
  b2 ~ gamma(shape_b2, rate_b2);
  v2 ~ gamma(shape_v2, rate_v2);
  
  for (i in 1:N) {
    // Model for x[i]
    target += log(b1/a1) + (b1 - 1) * log(x[i] / a1) - 2 * pow(x[i] / a1, b1)
              + v1 * (-1 + exp(-pow(x[i] / a1, b1))) + log(v1 + exp(pow(x[i] / a1, b1)));
    // Model for y[i]
    target += log(b2/a2) + (b2 - 1) * log(y[i] / a2) - 2 * pow(y[i] / a2, b2)
              + v2 * (-1 + exp(-pow(y[i] / a2, b2))) + log(v2 + exp(pow(y[i] / a2, b2)));
  }
}
"

shape_a1 <- 9
rate_a1 <- 5
shape_b1 <- 2
rate_b1 <- 1
shape_v1 <- 1
rate_v1 <- 5
shape_a2 <- 9
rate_a2 <- 5
shape_b2 <- 2
rate_b2 <- 1
shape_v2 <- 1
rate_v2 <- 3


fit <- stan(
  model_code = model_code,
  data = list(
    N = length(x),
    M = length(y),
    x = x,
    y=y, 
    shape_a1 = shape_a1,
    rate_a1 = rate_a1,
    shape_b1 = shape_b1,
    rate_b1 = rate_b1,
    shape_v1 = shape_v1,
    rate_v1 = rate_v1,
    shape_a2 = shape_a2,
    rate_a2 = rate_a2,
    shape_b2 = shape_b2,
    rate_b2 = rate_b2,
    shape_v2 = shape_v2,
    rate_v2 = rate_v2
  ),
  chains = 5,
  iter = 20000,
  warmup = 1000,
  control = list(adapt_delta = 0.95)
)

# Summarize the results
print(fit)

posterior <- as.matrix(fit)

posterior

plot_title <- ggtitle("Posterior distributions of parameters",
                      "with medians and 80%intervals")

mcmc_areas(posterior,
           pars = c("a1", "b1", "v1", "a2", "b2", "v2"),
           prob = 0.8) 



########################
posterior2 <- extract(fit, inc_warmup = TRUE, permuted = FALSE)

color_scheme_set("teal")
p <- mcmc_trace(posterior2, pars = c("a1", "b1", "v1", "a2", "b2", "v2"), n_warmup = 1000,
                facet_args = list(nrow = 2, ncol = 3, labeller = label_parsed))
p + facet_text(size = 16)
