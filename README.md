
# BINtools

<!-- badges: start -->
<!-- badges: end -->

The goal of BINtools is to implement a BIN model, a Bayesian approach to decomposing forecasting accuracy into three components: bias, partial information, and noise. 

## Installation

You can install the released version of BINtools from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("BINtools")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(BINtools)
# An example with two forecasting groups
# a) Simulate synthetic data:
synthetic_data = simulate_data(list(mu_star = -0.8,mu_0 = -0.5,mu_1 = 0.2,gamma_0 = 0.1,
gamma_1 = 0.3, rho_0 = 0.05,delta_0 = 0.1, rho_1 = 0.2, delta_1 = 0.3,rho_01 = 0.05), 300,100,100)
# b) Estimate the BIN-model on the synthetic data:
full_bayesian_fit = estimate_BIN(synthetic_data$Outcomes,synthetic_data$Control,synthetic_data$Treatment,warmup = 1000, iter = 2000)
# c) Analyze the results:
complete_summary(full_bayesian_fit)

```

