# Aggregated Cauchy Association Test

A p-value combination method using the Cauchy distribution.

## Usage

``` r
acat(Pvals, weights = NULL, is.check = TRUE)
```

## Arguments

- Pvals:

  a numeric vector/matrix of p-values. When it is a matrix, each column
  of p-values is combined by ACAT.

- weights:

  a numeric vector/matrix of non-negative weights for the combined
  p-values. When it is NULL, the equal weights are used.

- is.check:

  logical. Should the validity of *Pvals* (and *weights*) be checked?
  When the size of *Pvals* is large and one knows *Pvals* is valid, then
  the checking part can be skipped to save memory.

## Value

The p-value(s) of acat

## References

Liu, Y., & Xie, J. (2019). Cauchy combination test: a powerful test with
analytic p-value calculation under arbitrary dependency structures.
*Journal of American Statistical Association*,115(529), 393-402.
([pub](https://amstat.tandfonline.com/doi/abs/10.1080/01621459.2018.1554485))

## Author

Yaowu Liu

## Examples

``` r
p.values<-c(2e-02,4e-04,0.2,0.1,0.8);acat(Pvals=p.values)
#> [1] 0.001953404
acat(matrix(runif(1000),ncol=10))
#>  [1] 0.69790512 0.97150685 0.62340647 0.01036855 0.38888582 0.56808756
#>  [7] 0.66476915 0.61309425 0.08586323 0.38686721
```
