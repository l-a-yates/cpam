# Compute p-values for each target ID

Compute p-values for each target ID

## Usage

``` r
compute_p_values(
  cpo,
  subset = NULL,
  p_adj_method = "BH",
  aggregation_method = "lancaster",
  gam_method = "REML",
  gam_optimizer = "efs",
  silent = TRUE
)
```

## Arguments

- cpo:

  a cpam object

- subset:

  a character vector of target_id names

- p_adj_method:

  method for p-value adjustment

- aggregation_method:

  method for aggregating target-level p-values to gene-level; either
  "lancaster" (default) or "acat"

- gam_method:

  fitting method for
  [`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html) (default is
  "REML")

- gam_optimizer:

  optimization method for
  [`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html) (default is
  "efs")

- silent:

  logical; silences warnings from model fitting (default is TRUE)

## Value

an updated cpam object with raw, adjusted, and possibly aggregated
p-values stored in the new slot "p_table"

## Details

This function computes p-values for each target_id in the supplied cpam
object. The p-values are computed from a negative binomial GAM model
with a thin-plate spline basis function(s) for time using the `mgcv`
package.

The p-values are stored in the new slot `p_table` in the cpam object. If
`aggregate_to_gene` is set to `TRUE` (default), the target p-values are
aggregated to the gene level using the `lancaster` method. The columns
`p_val_target` and `p_val_gene` store the raw p-values for target- and
gene-level, respectively. The function also computes adjusted p-values
using the `p_adj_method`. The default method is "BH"
(Benjamini-Hochberg), but any methods supported by the function
`p.adjust` can be used. The adjusted p-values are stored in the columns
`q_val_target` and `q_val_gene`.

## References

Wood, S.N. (2013a) On p-values for smooth components of an extended
generalized additive model. Biometrika 100:221-228
doi:10.1093/biomet/ass048

Yi L, Pachter L (2018). aggregation: p-Value Aggregation Methods. R
package version 1.0.1, https://CRAN.R-project.org/package=aggregation.

## Examples

``` r
library(cpam)

# load gene-only example cpam object
load(system.file("extdata", "cpo_example.rda", package = "cpam"))

# run on a small subset of the example data
cpo <- compute_p_values(cpo_example, subset = paste0("g00",1:9))
cpo$p_table
#> # A tibble: 9 × 4
#>   target_id counts_mean p_val_target q_val_target
#>   <chr>           <dbl>        <dbl>        <dbl>
#> 1 g001           92078.   4.65 e-  1    8.37e-  1
#> 2 g002            8589.   9.00 e-  1    9.09e-  1
#> 3 g003            2147.   1    e-319    9   e-319
#> 4 g004           25248.   8.93 e-  1    9.09e-  1
#> 5 g005          116963.   9.09 e-  1    9.09e-  1
#> 6 g006            2122.   6.62 e-  1    9.09e-  1
#> 7 g007            3337.   3.74 e-  1    8.37e-  1
#> 8 g008            4646.   2.51 e-  1    7.54e-  1
#> 9 g009           37448.   2.37 e-  1    7.54e-  1
```
