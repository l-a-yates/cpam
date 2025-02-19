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

  # Ensure formula is properly formatted
  if (!inherits(fixed_effects, "formula")) {
    stop("fixed_effects must be a formula (e.g., ~ effect1 + effect2)")
  }

  # Validate that formula terms exist in exp_design
  fe_terms <- attr(stats::terms(fixed_effects), "term.labels")

  if (length(fe_terms) == 0) {
    return(NULL)    }

  # Check all terms exist in exp_design
  missing_terms <- fe_terms[!fe_terms %in% names(exp_design)]
  if (length(missing_terms) > 0) {
    stop(paste0("The following fixed effects variables are missing from exp_design: ",
                paste(missing_terms, collapse = ", ")))
  }

  # Check for perfect collinearity with time
  time_collinear <- check_collinearity_with_time(exp_design, fe_terms)
  if (time_collinear) {
    warning("Some fixed effects may be collinear with time, which could affect model estimation",
            call. = FALSE)
  }

  # Check for perfect separation
  check_perfect_separation(exp_design, fe_terms)

  # Return terms as character string
  paste(fe_terms, collapse = " + ")
}

#' Check for collinearity between fixed effects and time variable
#'
#' @param exp_design Experimental design data frame
#' @param fe_terms Fixed effects terms to check
#' @return Logical indicating whether collinearity was detected
#' @keywords internal
check_collinearity_with_time <- function(exp_design, fe_terms) {
  # For simpler cases, check direct correlation
  if (length(fe_terms) == 1) {
    # For categorical variables
    if (is.factor(exp_design[[fe_terms]]) || is.character(exp_design[[fe_terms]])) {
      # Check if time perfectly determines the fixed effect
      return(length(unique(paste0(exp_design$time, "_", exp_design[[fe_terms]]))) ==
               length(unique(exp_design$time)))
    } else {
      # For numeric variables, check correlation
      cor_val <- stats::cor(exp_design$time, exp_design[[fe_terms]])
      return(abs(cor_val) > 0.95) # High correlation threshold
    }
  }
  return(FALSE)
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

  if (length(categorical_terms) > 0) {
    for (term in categorical_terms) {
      levels_count <- table(exp_design[[term]])
      if (any(levels_count == 1)) {
        warning(paste0("Fixed effect '", term,
                       "' has levels with only one observation, which may cause estimation issues"),
                call. = FALSE)
      }
    }
  }
}
