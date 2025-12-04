# Launches a Shiny app to visualise the data and fitted models of a cpam object

Launches a Shiny app to visualise the data and fitted models of a cpam
object

## Usage

``` r
visualise(
  cpo,
  subset = NULL,
  degs_only = TRUE,
  deg_threshold = 0.05,
  p_type = c("p_gam", "p_mvn"),
  shape_type = c("shape1", "shape2")
)
```

## Arguments

- cpo:

  A cpam object containing count data, model fits, and optional
  changepoint/shape estimates

- subset:

  Character vector; names of targets or genes (if
  `cpo$gene_level = TRUE`) to load into the Shiny app. If NULL, all
  genes/targets are included based on `degs_only`.

- degs_only:

  Logical; if TRUE, display only differentially expressed genes/targets
  with adjusted p-value below `deg_threshold`. Default is TRUE.

- deg_threshold:

  Numeric; significance threshold for differentially expressed
  genes/targets. Only used when `degs_only = TRUE`. Default is 0.05.

- p_type:

  character; choose the type of p-value. Options are "p_gam" (default)
  or "p_mvn" (see
  [`compute_p_values()`](https://l-a-yates.github.io/cpam/reference/compute_p_values.md)
  for details).

- shape_type:

  character; "shape1" to include unconstrained or otherwise "shape2".
  Default is "shape1". In some instances, all of the transcripts for a
  gene may be "null" shaped, but the p-value for the gene may still be
  significant. This is due to the different methods of determining
  significance for the changepoints and the gene-level p-values. Here,
  conservatively, we remove these null-shaped genes from the DEG list.

## Value

None (launches Shiny app in browser)

## Examples

``` r
if (interactive()){

# A simple gene-level example
cpo <- cpo_example

# Launch visualization with all genes
visualise(cpo, degs_only = FALSE)

# Launch with only significant genes
visualise(cpo, deg_threshold = 0.05)

# Launch with specific genes
visualise(cpo, subset = c("g001","g002","g003"))
}
```
