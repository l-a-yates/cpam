# Use model selection to estimate changepoints

Use model selection to estimate changepoints

## Usage

``` r
estimate_changepoint(
  cpo,
  cps = NULL,
  degs_only = TRUE,
  deg_threshold = 0.05,
  subset = NULL,
  sp = NULL,
  bss = "tp",
  family = c("nb", "gaussian"),
  score = "aic",
  compute_mvn = TRUE
)
```

## Arguments

- cpo:

  a cpam object

- cps:

  vector of candidate changepoints. Defaults to the set of observed
  timepoints

- degs_only:

  logical; should changepoints only be estimated for differentially
  expressed genes

- deg_threshold:

  logical; the threshold value for DEGs (ignored of `degs_only = F`)

- subset:

  character vector; names of targets or genes (if `cpo$gene_level = T`)
  for which changepoints will be estimated

- sp:

  numerical \>= 0; supply a fixed smoothing parameter. This can decrease
  the fitting time but it is not recommended as changepoints estimation
  is sensitive to smoothness.

- bss:

  character vector; names of candidate spline bases (i.e., candidate
  shape types). Default is thin plate ("tp") splines.

- family:

  character; negative binomial ("nb", default) or Gaussian ("gaussian",
  not currently supported)

- score:

  character; model selection score, either Generalised Cross Validation
  ("gcv") or Akaike Information Criterion ("aic")

- compute_mvn:

  Use simulation to compute p-value under multivariate normal model of
  the model scores

## Value

a cpam object with the estimated changepoint table added to the slot
"changepoints"

## Details

This function estimates changepoints for each target_id. The assumed
trajectory type for this modelling stage is initially constant followed
by a changepoint into thin-plate smoothing spline.

By default, candidate time points are limited to the discrete observed
values in the series, since, despite the use of smoothing constraints,
there is generally insufficient information to infer the timing of
changepoints beyond the temporal resolution of the data. In any case,
the candidate points can be set manually using the `cps` argument.

To estimate changepoints, a model is fit for each candidate changepoint
and generalised cross-validation (GCV, default) or the Akaike
Information Criterion (AIC) are used to select among them.
Model-selection uncertainty is dealt with by computing the
one-standard-error rule, which identifies the least complex model within
one standard error of the best scoring model.

Both the minimum and the one-standard-error (default) models are stored
in the returned slot "changepoints" so that either can be used. In
addition to these, this function also computes the probability (denoted
`p_mvn`) that the null model is the best scoring model, using a
simulation based approach based on the multivariate normal model of the
pointwise model scores.

Given the computational cost of fitting a separate model for each
candidate changepoint, cpam only estimates changepoints for targets
associated with 'significant' genes at the chosen threshold
`deg_threshold`.

## References

Yates, L. A., S. A. Richards, and B. W. Brook. 2021. Parsimonious model
selection using information theory: a modified selection rule. Ecology
102(10):e03475. 10.1002/ecy.3475

## Examples

``` r
library(cpam)
library(dplyr)
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union

# load example data
load(system.file("extdata", "exp_design_example.rda", package = "cpam"))
load(system.file("extdata", "count_matrix_example.rda", package = "cpam"))

cpo <- prepare_cpam(exp_design = exp_design_example,
                    count_matrix = count_matrix_example[1:20,],
                    gene_level = TRUE,
                    num_cores = 1)
#> ℹ Processing count matrix
#> ✔ Processing count matrix [17ms]
#> 
#> ℹ Filtering low count genes
#> ℹ Estimating dispersions using edgeR
#> ✔ Estimating dispersions using edgeR [48ms]
#> 
#> ℹ Filtering low count genes
#> ✔ Filtering low count genes [90ms]
#> 
cpo <- compute_p_values(cpo)
cpo <- estimate_changepoint(cpo)
#> Estimating changepoints for 3 targets
#> Candidate changepoints are t = 1, 2, 3, 4, 5, and 6
cpo$changepoints
#> # A tibble: 3 × 7
#>   target_id cp_min cp_1se  p_mvn bs    family score_table      
#>   <chr>      <dbl>  <dbl>  <dbl> <chr> <chr>  <list>           
#> 1 g003           1      1 0      tp    nb     <tibble [30 × 6]>
#> 2 g013           3      3 0      tp    nb     <tibble [30 × 6]>
#> 3 g014           1      5 0.0767 tp    nb     <tibble [30 × 6]>
```
