# Changelog

## cpam 0.2.1

- Fixed bug where
  [`lancaster()`](https://l-a-yates.github.io/cpam/reference/lancaster.md)
  was not exported, causing `match.fun("lancaster")` to fail for
  installed package users
- Replaced fragile
  [`match.fun()`](https://rdrr.io/r/base/match.fun.html) lookup with
  internal resolution in
  [`compute_p_values()`](https://l-a-yates.github.io/cpam/reference/compute_p_values.md)
- Added `aggregation_method` parameter to
  [`estimate_changepoint()`](https://l-a-yates.github.io/cpam/reference/estimate_changepoint.md)
  for consistency with
  [`compute_p_values()`](https://l-a-yates.github.io/cpam/reference/compute_p_values.md)
- Exported
  [`lancaster()`](https://l-a-yates.github.io/cpam/reference/lancaster.md)
  with documentation
- Added tests for p-value aggregation methods

## cpam 0.2.0

CRAN release: 2026-03-04

- Added [`acat()`](https://l-a-yates.github.io/cpam/reference/acat.md)
  function for p-value aggregation
- Added `aggregation_method` parameter to
  [`compute_p_values()`](https://l-a-yates.github.io/cpam/reference/compute_p_values.md)
  to allow users to choose between “lancaster” (default) and “acat”
  methods
- Improved performance by reusing GAM design matrices in
  [`compute_p_values()`](https://l-a-yates.github.io/cpam/reference/compute_p_values.md)
  and
  [`estimate_changepoint()`](https://l-a-yates.github.io/cpam/reference/estimate_changepoint.md)
- Changed `compute_mvn` default to `FALSE` in
  [`estimate_changepoint()`](https://l-a-yates.github.io/cpam/reference/estimate_changepoint.md)
  to skip costly multivariate normal simulation by default
- Fixed tidyselect deprecation warnings

## cpam 0.1.3

CRAN release: 2025-03-13

- First release of cpam on CRAN.
- Implements changepoint detection, smooth and shape-constrained trends
  for time series omics data.
- Supports transcript- and gene-level analysis, case-only and
  case-control studies, and uncertainty quantification.
- Includes an interactive ‘shiny’ interface.
