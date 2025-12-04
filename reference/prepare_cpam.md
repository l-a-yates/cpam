# Prepare a cpam object

Prepare a cpam object

## Usage

``` r
prepare_cpam(
  exp_design,
  count_matrix = NULL,
  t2g = NULL,
  import_type = NULL,
  model_type = c("case-only", "case-control"),
  bootstrap = TRUE,
  filter_fun = "ts_filter",
  filter_fun_args = list(min_reads = 5, min_prop = 3/5),
  regularize = TRUE,
  gene_level = FALSE,
  aggregate_to_gene = !gene_level,
  condition_var = "condition",
  case_value = "treatment",
  num_cores = 1,
  normalize = TRUE,
  fixed_effects = NULL,
  intercept_cc = c("1", condition_var)
)
```

## Arguments

- exp_design:

  a dataframe or tibble with the experimental design, containing at
  least a 'time' and a 'sample' column

- count_matrix:

  a matrix of counts. Column names must be in 'sample' column of
  `exp_design`,

- t2g:

  a transcript to gene dataframe or tibble with columns target_id and
  gene_id

- import_type:

  software used for quantification, one of "kallisto", "salmon" ,....
  Ignored if `count_matrix` is supplied.

- model_type:

  "case-only" (default) or "case-control"

- bootstrap:

  logical; load bootstrap samples, also called inferential replicates,
  if available, and rescale counts.

- filter_fun:

  filter function to remove lowly expressed genes (default is
  `filter_fun()`)

- filter_fun_args:

  arguments for filter function

- regularize:

  logical; use empirical Bayes regularization of dispersions (default is
  TRUE)

- gene_level:

  logical; aggregate counts to gene level before data preparation and
  modelling (default is FALSE)

- aggregate_to_gene:

  logical; aggregate p values from transcript- to gene-level

- condition_var:

  string; column name in `exp_design` for the condition variable (for
  `model_type` = "case_control" only)

- case_value:

  value of `condition_var` that indicates the "case". All other values
  are deemed to be control

- num_cores:

  integer; number of cores to use for parallel computation

- normalize:

  logical; use model offsets based on sampling depth and gene length

- fixed_effects:

  a model formula of the form `~ effect1 + effect2`

- intercept_cc:

  string; intercept for case-control model: "1" (default) for common
  intercept or "condition"

## Value

an object of class
[`cpam-class`](https://l-a-yates.github.io/cpam/reference/cpam-class.md).
The returned object has methods `print` and `summary` for displaying
information. See
[`cpam-class`](https://l-a-yates.github.io/cpam/reference/cpam-class.md)
for details on the structure of the returned object.

## Details

This function prepares a cpam object for analysis. The function loads
count data from files or a matrix, filters lowly expressed genes,
computes normalisation factors, and estimates dispersions. Many of these
steps can be customised or turned off.

When bootstrap samples (inferential replicates) are available, it loads
and summarises these using means, standard errors, and estimated
overdispersions. The latter are a measure of quantification uncertainty
and they are used to rescale the counts which accounts for this
uncertainty during the modelling steps.

## References

Pedro L Baldoni, Yunshun Chen, Soroor Hediyeh-zadeh, Yang Liao, Xueyi
Dong, Matthew E Ritchie, Wei Shi, Gordon K Smyth, Dividing out
quantification uncertainty allows efficient assessment of differential
transcript expression with edgeR, Nucleic Acids Research, Volume 52,
Issue 3, 9 February 2024, Page e13, https://doi.org/10.1093/nar/gkad1167

Yunshun Chen, Lizhong Chen, Aaron T L Lun, Pedro L Baldoni, Gordon K
Smyth, edgeR v4: powerful differential analysis of sequencing data with
expanded functionality and improved support for small counts and larger
datasets, Nucleic Acids Research, Volume 53, Issue 2, 27 January 2025,
https://doi.org/10.1093/nar/gkaf018

## Examples

``` r
library(cpam)

# load gene-only example data
load(system.file("extdata", "exp_design_example.rda", package = "cpam"))
load(system.file("extdata", "count_matrix_example.rda", package = "cpam"))

cpo <- prepare_cpam(exp_design = exp_design_example,
                    count_matrix = count_matrix_example,
                    gene_level = TRUE)
#> ℹ Processing count matrix
#> ✔ Processing count matrix [20ms]
#> 
#> ℹ Filtering low count genes
#> ℹ Estimating dispersions using edgeR
#> ✔ Estimating dispersions using edgeR [571ms]
#> 
#> ℹ Filtering low count genes
#> ✔ Filtering low count genes [730ms]
#> 
cpo
#> 
#> ── cpam object ─────────────────────────────────────────────────────────────────
#> • case-only time series
#> • 30 samples
#> • 6 time points
#> • Counts aggregated for gene-level inference
```
