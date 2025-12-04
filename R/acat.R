# This function is taken from the ACAT package to save a package dependency.
# Source: <https://github.com/yaowuliu/ACAT/blob/master/R/ACAT.R>
# Publication: <https://doi.org/10.1016/j.ajhg.2019.01.002>
#
# The original code is licensed under GPL-3
# This modified version retains the original license and attribution.

#'
#' Aggregated Cauchy Association Test
#'
#' A p-value combination method using the Cauchy distribution.
#'
#'
#'
#' @param weights a numeric vector/matrix of non-negative weights for the combined p-values. When it is NULL, the equal weights are used.
#' @param Pvals a numeric vector/matrix of p-values. When it is a matrix, each column of p-values is combined by ACAT.
#' @param is.check logical. Should the validity of \emph{Pvals} (and \emph{weights}) be checked? When the size of \emph{Pvals} is large and one knows \emph{Pvals} is valid, then the checking part can be skipped to save memory.
#' @return The p-value(s) of acat
#' @author Yaowu Liu
#' @examples p.values<-c(2e-02,4e-04,0.2,0.1,0.8);acat(Pvals=p.values)
#' @examples acat(matrix(runif(1000),ncol=10))
#' @references Liu, Y., & Xie, J. (2019). Cauchy combination test: a powerful test with analytic p-value calculation
#' under arbitrary dependency structures. \emph{Journal of American Statistical Association},115(529), 393-402. (\href{https://amstat.tandfonline.com/doi/abs/10.1080/01621459.2018.1554485}{pub})
#' @importFrom stats pcauchy
#' @export
acat<-function(Pvals,weights=NULL,is.check=TRUE){
  Pvals<-as.matrix(Pvals)
  if (is.check){
    #### check if there is NA
    if (sum(is.na(Pvals))>0){
      stop("Cannot have NAs in the p-values!")
    }
    #### check if Pvals are between 0 and 1
    if ((sum(Pvals<0)+sum(Pvals>1))>0){
      stop("P-values must be between 0 and 1!")
    }
    #### check if there are pvals that are either exactly 0 or 1.
    is.zero<-(colSums(Pvals==0)>=1)
    is.one<-(colSums(Pvals==1)>=1)
    if (sum((is.zero+is.one)==2)>0){
      stop("Cannot have both 0 and 1 p-values in the same column!")
    }

    if (sum(is.zero)>0){
      warning("There are p-values that are exactly 0!")
    }
    if (sum(is.one)>0){
      warning("There are p-values that are exactly 1!")
    }

  }
  #### Default: equal weights. If not, check the validity of the user supplied weights and standadize them.
  if (is.null(weights)){
    is.weights.null<-TRUE
  }else{
    is.weights.null<-FALSE
    weights<-as.matrix(weights)
    if (sum(dim(weights)!=dim(Pvals))>0){
      stop("The dimensions of weights and Pvals must be the same!")
    }else if (is.check & (sum(weights<0)>0)){
      stop("All the weights must be nonnegative!")
    }else{
      w.sum<-colSums(weights)
      if (sum(w.sum<=0)>0){
        stop("At least one weight should be positive in each column!")
      }else{
        for (j in 1:ncol(weights)){
          weights[,j]<-weights[,j]/w.sum[j]
        }
      }
    }

  }

  #### check if there are very small non-zero p values and calcuate the cauchy statistics
  is.small<-(Pvals<1e-15)
  if (is.weights.null){
    Pvals[!is.small]<-tan((0.5-Pvals[!is.small])*pi)
    Pvals[is.small]<-1/Pvals[is.small]/pi
    cct.stat<-colMeans(Pvals)
  }else{
    Pvals[!is.small]<-weights[!is.small]*tan((0.5-Pvals[!is.small])*pi)
    Pvals[is.small]<-(weights[is.small]/Pvals[is.small])/pi
    cct.stat<-colSums(Pvals)
  }
  #### return the ACAT p value(s).
  pval<-pcauchy(cct.stat,lower.tail = F)
  return(pval)
}
