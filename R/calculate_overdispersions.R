# This function implements the method presented in https://doi.org/10.1093/nar/gkad1167
# to calculate overdispersions due to quantification uncertainty using inferential replicates.
#
# It is adapted from the code provided in the supplementary material of the paper at
# https://plbaldoni.rbind.io/TranscriptDE-code/casestudy.html#Differential_gene_expression_analysis
# to allow it to be applied directly to list objects returned from tximport.
#
# The original code is licensed under the MIT License. The license text is included below.
#
# The MIT License (MIT)
# Copyright (c) 2023 Pedro L. Baldoni
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the “Software”), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
# of the Software, and to permit persons to whom the Software is furnished to
# do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
calculate_overdispersions <- function(txi){
  cli::cli_progress_step("Calculating overdispersions from inferential replicates")
  NTx <- txi$counts %>% nrow
  NBoot <- txi$infReps[[1]] %>% ncol
  NSamples <- txi$counts %>% ncol
  DF <- rep_len(0L, NTx)
  OverDisp <- rep_len(0, NTx)

  for(j in 1L:NSamples){
    Boot <- txi$infReps[[j]]
    M <- rowMeans(Boot)
    i <- (M > 0)
    OverDisp[i] <- OverDisp[i] + rowSums((Boot[i, ] -
                                            M[i])^2)/M[i]
    DF[i] <- DF[i] + NBoot - 1L
  }

  i <- (DF > 0L)
  if (sum(i) > 0L) {
    OverDisp[i] <- OverDisp[i]/DF[i]
    DFMedian <- stats::median(DF[i])
    DFPrior <- 3
    OverDispPrior <- stats::median(OverDisp[i])/stats::qf(0.5, df1 = DFMedian,
                                            df2 = DFPrior)
    if (OverDispPrior < 1)
      OverDispPrior <- 1
    OverDisp[i] <- (DFPrior * OverDispPrior + DF[i] * OverDisp[i])/(DFPrior +
                                                                      DF[i])
    OverDisp <- pmax(OverDisp, 1)
    OverDisp[!i] <- OverDispPrior
  }  else {
    OverDisp[] <- NA_real_
    OverDispPrior <- NA_real_
  }
  names(OverDisp) <- rownames(txi$counts)

  list(overdispersion = OverDisp,
       overdispersion.prior = OverDispPrior)

}
