# this needs citing etc.
calculate_overdispersions <- function(txi){
  message("Calculating overdispersions from inferential replicates")
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
