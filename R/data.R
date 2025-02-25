#' Example Experimental Design
#'
#' A tibble representing an example experimental design for cpam analysis.
#'
#' @format A tibble with columns:
#' \describe{
#'   \item{sample}{Character vector of sample names.}
#'   \item{time}{Numeric vector of time points.}
#'   \item{path}{Character vector of file paths.}
#'   ...
#' }
#' @source Simulated example.
#' @concept dataset
"exp_design_path"

#' Transcript-to-Gene Mapping
#'
#' A example tibble representing a mapping of transcript IDs to gene IDs.
#'
#' @format A tibble with columns:
#' \describe{
#'   \item{transcript}{Character vector of transcript IDs.}
#'   \item{gene}{Character vector of corresponding gene IDs.}
#' }
#' @source The TAIR database.
#' @concept dataset
"t2g_arabidopsis"


#' A small case-only count matrix to run light-weight examples
#'
#' See `cpo_example` for fitted CPAM object.
#' See `exp_design_example` for experimental design.
#' 1000 genes x 6 time points x 5 reps
#'
#' @format A matrix with 1000 rows and 30 columns:
#' \describe{
#'   \item{rownames}{Gene ID}
#'   \item{colnames}{Samples}
#' }
#' @source Simulated data, see data-raw directory for details
#' @concept dataset
"count_matrix_example"

#' An experimental design to run light-weight examples
#'
#' Six time points x 3 reps
#' #' See `cpo_example` for fitted CPAM object.
#'
#'#'
#' @format A tibble with columns:
#' \describe{
#'   \item{sample}{Character vector of sample names.}
#'   \item{time}{Numeric vector of time points.}
#'   ...
#' }
#' @source Simulated data, see data-raw directory for details
#' @concept dataset
"exp_design_example"


#' A fitted `cpam` object to run light-weight examples
#'
#' Six time points x 3 reps
#' See `cpo_example` for fitted CPAM object.
#' See `exp_design_example` for experimental design.
#'#'
#' @format A cpam object with p-values, changpoints, and shapes estimated
#'
#' @source Simulated data, see data-raw directory for details
#' @concept dataset
"cpo_example"
