# Use model selection to select a shape for each target

Use model selection to select a shape for each target

## Usage

``` r
select_shape(
  cpo,
  subset = NULL,
  sp = NULL,
  bss = c("micv", "mdcx", "cv", "cx", "lin", "tp", "null"),
  family = c("nb", "gaussian"),
  score = "gcv",
  cp_type = c("cp_1se", "cp_min")
)
```

## Arguments

- cpo:

  a cpam object

- subset:

  character vector; names of targets or genes (if `cpo$gene_level = T`)
  for which changepoints will be estimated

- sp:

  numerical \>= 0; supply a fixed smoothing parameter. If NULL
  (default), the smoothing parameter is estimated. Note, this fixed
  value is in any case applied only to shape constrained bases (i.e.,
  not `bs = 'tp'`).

- bss:

  character vector; names of candidate spline bases (i.e., candidate
  shape types).

- family:

  character; negative binomial ("nb", default) or Gaussian ("gaussian")

- score:

  character; model selection score, either Generalised Cross Validation
  ("gcv") or Akaike Information Criterion ("aic")

- cp_type:

  character; if changepoints have been estimated using
  [`estimate_changepoint()`](https://l-a-yates.github.io/cpam/reference/estimate_changepoint.md),
  which selection rule should be used. See
  [`estimate_changepoint()`](https://l-a-yates.github.io/cpam/reference/estimate_changepoint.md)
  for details.

## Value

a cpam object with the selected shapes added to the slot "shapes"

## Details

The function selects the best shape from a list of candidate shapes for
each target. It is typically the last step in the analysis, called after
p-values have been estimated using
[`compute_p_values()`](https://l-a-yates.github.io/cpam/reference/compute_p_values.md)
and changepoints have been estimated using
[`estimate_changepoint()`](https://l-a-yates.github.io/cpam/reference/estimate_changepoint.md).

Two shape selections are generated. The first selecting among linear,
convex and concave shape classes and their monotonic variants (or the
shape set given by bss), and the second selecting among the first
options plus an 'unconstrained' smooth. The inclusion of the
'unconstrained' type provides the flexibility to detect targets beyond
simpler trends. For computational reasons, as per the changepoint
estimation, shapes are selected only for those genes, or their isoforms,
identified as significant at the chosen FDR threshold. This is
overridden by providing a subset of target names to the `subset`
argument, provided these targets have estimated changepoints.

## Examples

``` r
library(cpam)

# load example data
load(system.file("extdata", "exp_design_example.rda", package = "cpam"))
load(system.file("extdata", "count_matrix_example.rda", package = "cpam"))

# Using a small subset of the example data
cpo <- prepare_cpam(exp_design = exp_design_example,
                    count_matrix = count_matrix_example[1:20,],
                    gene_level = TRUE,
                    num_cores = 1)
#> ℹ Processing count matrix
#> ✔ Processing count matrix [23ms]
#> 
#> ℹ Filtering low count genes
#> ℹ Estimating dispersions using edgeR
#> ✔ Estimating dispersions using edgeR [47ms]
#> 
#> ℹ Filtering low count genes
#> ✔ Filtering low count genes [83ms]
#> 
cpo <- compute_p_values(cpo)
cpo <- estimate_changepoint(cpo)
#> Estimating changepoints for 3 targets
#> Candidate changepoints are t = 1, 2, 3, 4, 5, and 6
cpo <- select_shape(cpo)
#> Estimating shapes for 3 targets
#> Candidate changepoints are bs = "micv", "mdcx", "cv", "cx", "lin", "tp", and
#> "null"
cpo$shapes
#> # A tibble: 3 × 4
#>   target_id    cp shape1 shape2
#>   <chr>     <dbl> <chr>  <chr> 
#> 1 g003          1 mdcx   mdcx  
#> 2 g013          3 ilin   ilin  
#> 3 g014          5 null   null  

```
