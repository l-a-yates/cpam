#' @export
print.cpam <- function(x, ...){
  cli::cli_h1("cpam object")
  cli::cli_ul()
  cli::cli_li(cli::style_dim("{x$model_type} time series"))
  cli::cli_li(cli::style_dim("{length(x$exp_design$sample)} samples"))
  cli::cli_li(cli::style_dim("{length(unique(x$exp_design$time))} time points"))
  if(x$gene_level) cli::cli_li(cli::style_dim("Counts aggregated for gene-level inference"))
  if(x$bootstrap){
    cli::cli_li(cli::style_dim("Overdispersion estimated using {x$nboot} inferential replicates"))
    cli::cli_li(cli::style_dim("Counts rescaled by estimated overdispersion"))
  }
}
