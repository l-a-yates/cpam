# cpam 0.2.1

* Fixed bug where `lancaster()` was not exported, causing `match.fun("lancaster")` to fail for installed package users
* Replaced fragile `match.fun()` lookup with internal resolution in `compute_p_values()`
* Added `aggregation_method` parameter to `estimate_changepoint()` for consistency with `compute_p_values()`
* Exported `lancaster()` with documentation
* Added tests for p-value aggregation methods

# cpam 0.2.0

* Added `acat()` function for p-value aggregation
* Added `aggregation_method` parameter to `compute_p_values()` to allow users to choose between "lancaster" (default) and "acat" methods
* Improved performance by reusing GAM design matrices in `compute_p_values()` and `estimate_changepoint()`
* Changed `compute_mvn` default to `FALSE` in `estimate_changepoint()` to skip costly multivariate normal simulation by default
* Fixed tidyselect deprecation warnings

# cpam 0.1.3

* First release of cpam on CRAN.
* Implements changepoint detection, smooth and shape-constrained trends for time series omics data.
* Supports transcript- and gene-level analysis, case-only and case-control studies, and uncertainty quantification.
* Includes an interactive 'shiny' interface.
