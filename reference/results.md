# Create a results table from a cpam object

Create a results table from a cpam object

## Usage

``` r
results(
  cpo,
  p_threshold = 0.05,
  p_type = c("p_gam", "p_mvn"),
  min_lfc = 0,
  min_count = 0,
  aggregate_to_gene = cpo$aggregate_to_gene,
  add_lfc = TRUE,
  add_counts = TRUE,
  cp_type = c("cp_1se", "cp_min"),
  shape_type = c("shape1", "shape2"),
  summarise_to_gene = FALSE,
  remove_null_targets = TRUE
)
```

## Arguments

- cpo:

  a cpam object

- p_threshold:

  numerical; threshold for adjusted p-values; default is 0.05

- p_type:

  character; choose the type of p-value. Options are "p_gam" (default)
  or "p_mvn" (see
  [`compute_p_values()`](https://l-a-yates.github.io/cpam/reference/compute_p_values.md)
  for details).

- min_lfc:

  numerical; maximum absolute log (base 2) fold change must exceed this
  minimum value; default is 0

- min_count:

  numerical; maximum of the modelled counts evaluated at the set of
  observed time points must exceed this minimum value for

- aggregate_to_gene:

  logical; filter by gene-aggregated p-values

- add_lfc:

  logical; add log (base 2) fold changes for each time point

- add_counts:

  logical; add modelled counts for each time point

- cp_type:

  character; model-selection rule used to select the changepoint

- shape_type:

  character; "shape1" to include unconstrained or otherwise "shape2"

- summarise_to_gene:

  logical; return gene-level results only

- remove_null_targets:

  logical; remove targets with null shapes (default is T). If F, targets
  with null shapes will be included if the aggregated p-value for the
  corresponding gene passes the specified filtering thresholds.

## Value

a tibble

## Details

This function is usually called after
[`compute_p_values()`](https://l-a-yates.github.io/cpam/reference/compute_p_values.md),
[`estimate_changepoint`](https://l-a-yates.github.io/cpam/reference/estimate_changepoint.md),
and `select_shape` have been run. The function has several useful
filters such as adjusted p-value thresholds, minimum log-fold changes,
and minimum counts.

## Examples

``` r
library(cpam)

# load gene-only example cpam object
load(system.file("extdata", "cpo_example.rda", package = "cpam"))

results(cpo_example)
#> # A tibble: 104 × 16
#>    target_id         p    cp shape lfc.1  lfc.2  lfc.3  lfc.4  lfc.5  lfc.6
#>    <chr>         <dbl> <dbl> <chr> <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>
#>  1 g003      9.26e-319     1 mdcx      0 -0.349 -0.634 -0.856 -1.02  -1.11 
#>  2 g013      9.26e-319     3 ilin      0  0      0      0.327  0.655  0.982
#>  3 g055      9.26e-319     1 ilin      0  0.195  0.389  0.584  0.778  0.973
#>  4 g063      9.26e-319     1 cv        0  0.595  0.937  1.03   0.861  0.443
#>  5 g069      9.26e-319     2 ilin      0  0      0.225  0.451  0.676  0.902
#>  6 g090      9.26e-319     2 ilin      0  0      0.239  0.477  0.716  0.954
#>  7 g106      9.26e-319     1 mdcx      0 -0.470 -0.778 -0.892 -0.902 -0.902
#>  8 g126      9.26e-319     1 ilin      0  0.225  0.450  0.675  0.900  1.12 
#>  9 g128      9.26e-319     1 ilin      0  0.209  0.418  0.628  0.837  1.05 
#> 10 g129      9.26e-319     1 cv        0  0.626  0.979  1.11   1.08   0.888
#> # ℹ 94 more rows
#> # ℹ 6 more variables: counts.1 <dbl>, counts.2 <dbl>, counts.3 <dbl>,
#> #   counts.4 <dbl>, counts.5 <dbl>, counts.6 <dbl>

# Add filters
results(cpo_example, p_threshold = 0.01, min_lfc = 1)
#> # A tibble: 22 × 16
#>    target_id         p    cp shape lfc.1  lfc.2  lfc.3  lfc.4 lfc.5 lfc.6
#>    <chr>         <dbl> <dbl> <chr> <dbl>  <dbl>  <dbl>  <dbl> <dbl> <dbl>
#>  1 g063      9.26e-319     1 cv        0  0.595  0.937  1.03  0.861 0.443
#>  2 g126      9.26e-319     1 ilin      0  0.225  0.450  0.675 0.900 1.12 
#>  3 g128      9.26e-319     1 ilin      0  0.209  0.418  0.628 0.837 1.05 
#>  4 g129      9.26e-319     1 cv        0  0.626  0.979  1.11  1.08  0.888
#>  5 g171      9.26e-319     1 cv        0  0.593  0.922  1.01  0.901 0.616
#>  6 g186      9.26e-319     1 cx        0 -0.233 -0.324 -0.209 0.219 1.05 
#>  7 g331      9.26e-319     3 ilin      0  0      0      0.339 0.678 1.02 
#>  8 g334      9.26e-319     1 micv      0  0.402  0.692  0.886 1.00  1.06 
#>  9 g341      9.26e-319     1 cx        0 -0.465 -0.593 -0.383 0.164 1.05 
#> 10 g393      9.26e-319     1 cv        0  0.432  0.747  0.945 1.02  0.987
#> # ℹ 12 more rows
#> # ℹ 6 more variables: counts.1 <dbl>, counts.2 <dbl>, counts.3 <dbl>,
#> #   counts.4 <dbl>, counts.5 <dbl>, counts.6 <dbl>
```
