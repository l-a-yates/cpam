#' @export
print.cpam <- function(x, ...){
  cat("cpam object\n-----------\n")
  cat(paste0(x$model_type, " time series\n"))
  cat(paste0(length(x$exp_design$sample), " samples\n"))
  cat(paste0(length(unique(x$exp_design$time)), " time points\n"))
  if(x$gene_level) cat("Counts aggregated for gene-level inference\n")
  if(x$bootstrap){
    cat(paste0("Overdispersion estimated using ",x$nboot, " inferential replicates\n"))
    cat(paste0("Counts rescaled by estimated overdispersion"))
  }
}
