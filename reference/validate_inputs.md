# Validate input parameters for prepare_cpam

Validate input parameters for prepare_cpam

## Usage

``` r
validate_inputs(
  exp_design,
  count_matrix,
  t2g,
  import_type,
  model_type,
  condition_var,
  case_value,
  fixed_effects,
  aggregate_to_gene,
  gene_level,
  import
)
```

## Arguments

- exp_design:

  Experimental design data frame

- count_matrix:

  Count matrix

- t2g:

  Transcript to gene mapping

- import_type:

  Import type

- model_type:

  Model type

- condition_var:

  Condition variable

- case_value:

  Case value

- fixed_effects:

  Fixed effects formula

- aggregate_to_gene:

  Whether to aggregate to gene level

- gene_level:

  Whether data is at gene level

- import:

  Whether to import data

## Value

NULL, throws error if validation fails
