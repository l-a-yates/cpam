# Changelog

## cpam (development version)

- Added [`acat()`](https://l-a-yates.github.io/cpam/reference/acat.md)
  function for p-value aggregation
- Added `aggregation_method` parameter to
  [`compute_p_values()`](https://l-a-yates.github.io/cpam/reference/compute_p_values.md)
  to allow users to choose between “lancaster” (default) and “acat”
  methods

## cpam 0.1.3

CRAN release: 2025-03-13

- First release of cpam on CRAN.
- Implements changepoint detection, smooth and shape-constrained trends
  for time series omics data.
- Supports transcript- and gene-level analysis, case-only and
  case-control studies, and uncertainty quantification.
- Includes an interactive ‘shiny’ interface.
