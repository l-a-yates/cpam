# Simulate p-values using multivariate normal distribution

Simulate p-values using multivariate normal distribution

## Usage

``` r
simulate_p_mvn(score_table, nsim = 10000, reg = 0.001)
```

## Arguments

- score_table:

  Table of model scores

- nsim:

  Number of simulations (default: 1e4)

- reg:

  Regularization parameter for covariance matrix (default: 1e-3)

## Value

Probability that the null model is the best scoring model

## Details

Uses Monte Carlo simulation to compute the probability that the null
model is the best scoring model under multivariate normal assumptions.
