#' Example Experimental Design
#'
#' A tibble representing an example experimental design for CPAM analysis.
#'
#' @format A tibble with columns:
#' \describe{
#'   \item{sample}{Character vector of sample names.}
#'   \item{time}{Numeric vector of time points.}
#'   \item{path}{Character vector of file paths.}
#'   ...
#' }
#' @source Simulated example.
"exp_design"

#' Transcript-to-Gene Mapping
#'
#' A tibble representing a mapping of transcript IDs to gene IDs.
#'
#' @format A tibble with columns:
#' \describe{
#'   \item{transcript}{Character vector of transcript IDs.}
#'   \item{gene}{Character vector of corresponding gene IDs.}
#' }
#' @source The TAIR database.
"t2g"
