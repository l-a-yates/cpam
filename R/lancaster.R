# This function is taken from the aggregation package to save a package dependency.
# Source: <https://cran.r-project.org/package=aggregation>
# Publication: <https://doi.org/10.1101%2F190199>
#
# The original code is licensed under GPL-3
# This modified version retains the original license and attribution.
#
# Copyright (C) 2018 Lynn Yi and Lior Pachter

lancaster <- function (pvalues, weights)
{
  if (length(weights) != length(pvalues)) {
    stop("Length of weights not equal to length of pvalues")
  }
  weights <- weights[!is.na(pvalues)]
  pvalues <- pvalues[!is.na(pvalues)]
  pvalues <- pvalues[weights > 0 | is.na(weights)]
  weights <- weights[weights > 0 | is.na(weights)]
  if (length(pvalues) == 0) {
    return(NA)
  }
  if (length(pvalues) == 1) {
    return(pvalues)
  }
  if (any(pvalues < 9.99988867182683e-320)) {
    warning("Extreme p-values around and below 10e-320 will produce a p-value of 0. Replace extreme p-values with 10e-320 to obtain an upper bound for the aggregated p-value.")
  }
  t <- sapply(1:length(pvalues), function(i) lts(pvalues[i],
                                                 weights[i]))
  t <- sum(t)
  p <- stats::pchisq(t, sum(weights), lower.tail = FALSE)
  p
}

lts <- function (pvalue, weight)
  {
    stats::qgamma(pvalue, shape = weight/2, scale = 2, lower.tail = FALSE)
  }

