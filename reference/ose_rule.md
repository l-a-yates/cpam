# Apply the one-standard-error rule for model selection

Apply the one-standard-error rule for model selection

## Usage

``` r
ose_rule(tab, nse = 1)
```

## Arguments

- tab:

  A table containing model scores

- nse:

  Number of standard errors to use (default: 1)

## Value

The selected model that is most parsimonious within nse standard errors
of the best model

## Details

Implements the one-standard-error rule which selects the most
parsimonious model within a specified number of standard errors of the
best-scoring model.
