# Wrapper for estimating dispersions

Wrapper for estimating dispersions

## Usage

``` r
estimate_dispersions_wrapper(
  count_matrix_filtered,
  exp_design,
  bootstrap,
  catch
)
```

## Arguments

- count_matrix_filtered:

  Filtered count matrix

- exp_design:

  Experimental design

- bootstrap:

  Whether bootstrap was used

- catch:

  Overdispersion calculations from bootstrap

## Value

Estimated dispersions
