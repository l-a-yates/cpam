# Lancaster's method for p-value aggregation

Combines multiple p-values using Lancaster's weighted generalization of
Fisher's method, based on the chi-squared distribution.

## Usage

``` r
lancaster(pvalues, weights)
```

## Arguments

- pvalues:

  A numeric vector of p-values to combine.

- weights:

  A numeric vector of non-negative weights (same length as `pvalues`).

## Value

A single aggregated p-value.

## Details

NA p-values are removed, and entries with zero or NA weights are
excluded. If a single p-value remains after filtering, it is returned
as-is.

## References

Yi L, Pachter L (2018). aggregation: p-Value Aggregation Methods. R
package version 1.0.1, <https://CRAN.R-project.org/package=aggregation>.

Lancaster, H.O. (1961). The combination of probabilities: an application
of orthonormal functions. *Australian Journal of Statistics*, 3, 20-33.

## Examples

``` r
lancaster(c(0.01, 0.05, 0.1), c(1, 1, 1))
#> [1] 0.004259304
```
