# This function is take from the DESeq2 package to save a package dependency.
# Source: <https://code.bioconductor.org/browse/DESeq2/RELEASE_3_20/>
#
# The original code is licensed under LGPL-3
# This modified version retains the original license and attribution.
#
# Copyright (C) Michael Love, Simon Anders, Wolfgang Huber
#
estimateSizeFactorsForMatrix <- function (counts){

  loggeomeans <- rowMeans(log(counts))

  if (all(is.infinite(loggeomeans))) {
    stop("every gene contains at least one zero, cannot compute log geometric means")
  }
  sf <- apply(counts, 2, function(cnts) {
      exp(stats::median((log(cnts) - loggeomeans)[is.finite(loggeomeans) & cnts > 0]))
    })

  sf
}
