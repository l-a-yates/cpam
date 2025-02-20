#' Validate fixed effects formula against experimental design
#'
#' @param fixed_effects A model formula (e.g., ~ effect1 + effect2) or NULL
#' @param exp_design A data frame containing the experimental design
#' @return The validated and processed fixed effects as character, or NULL
#' @keywords internal
validate_fixed_effects <- function(fixed_effects, exp_design) {
  if (is.null(fixed_effects)) {
    return(NULL)
  }

  if (!inherits(fixed_effects, "formula")) {
    stop("fixed_effects must be a formula (e.g., ~ effect1 + effect2)")
  }

  fe_terms <- attr(stats::terms(fixed_effects), "term.labels")

  if (length(fe_terms) == 0) {
    return(NULL)    }

  missing_terms <- fe_terms[!fe_terms %in% names(exp_design)]
  if (length(missing_terms) > 0) {
    stop(paste0("The following fixed effects variables are missing from exp_design: ",
                paste(missing_terms, collapse = ", ")))
  }

  collinear <- check_collinearity(exp_design, fixed_effects)
  if (collinear) {
    stop("The model matrix is not full rank. This is likely due to collinear fixed effects.",
            call. = FALSE)
  }

  check_perfect_separation(exp_design, fe_terms)
  paste(fe_terms, collapse = " + ")
}


#' Check for collinearity between predictors and with time variable
#'
#' @param formula Model formula
#' @param exp_design Experimental design data frame
#' @return Logical indicating whether collinearity was detected
check_collinearity <- function(exp_design, formula){
  design_matrix <- stats::model.matrix(formula, exp_design)
  design_matrix_time <- cbind(design_matrix, exp_design[["time"]])
  qr(design_matrix_time)$rank < ncol(design_matrix_time)
}


#' Check for perfect separation in categorical fixed effects
#'
#' @param exp_design Experimental design data frame
#' @param fe_terms Fixed effects terms to check
#' @return NULL, issues warning if perfect separation detected
#' @keywords internal
check_perfect_separation <- function(exp_design, fe_terms) {
  categorical_terms <- fe_terms[sapply(exp_design[fe_terms], function(x) {
    is.factor(x) || is.character(x)
  })]

  ps <- c()
  if (length(categorical_terms) > 0) {
    for (term in categorical_terms) {
      levels_count <- table(exp_design[[term]])
      if (any(levels_count == 1)) {
        ps <- c(ps, term)
      }
    }
    if (length(ps) > 0) {
      warning(paste0("Perfect separation detected in the following terms: ",
                     paste(ps, collapse = ", ")),
              call. = FALSE)
    }
  }
}

