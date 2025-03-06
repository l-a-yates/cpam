#' Plot clustered targets
#'
#' @param cpo a cpam object
#' @param res a tibble, output from [results()] containing columns target_id, cp, and shape
#' @param changepoints numerical or character; one or more changepoints (these should be the same as the ones used in [estimate_changepoint()]
#' @param shapes character; one or more shapes (these should be the same as the ones used in [select_shape()]
#' @param alpha numeric between 0 and 1; controls line transparency in plot (default: 0.1)

#'
#' @details
#' Plots the fitted trends for a set of targets whose estimated changepoints and
#' shapes are given by the arguments `changepoints` and `shapes`, respectively.
#'
#' Creates a combined plot showing fitted expression trends for all targets that
#' share specified changepoint times and shape patterns. Each line represents one
#' target's fitted trajectory, with transparency controlled by alpha.
#'
#' @return A ggplot object showing overlaid fitted trends, or NULL if no matching
#'   targets are found
#' @export
#'
#' @seealso [`results()`], [`plot_cpam()`]
#' @examples
#' library(cpam)
#'
#' # load gene-only example cpam object
#' load(system.file("extdata", "cpo_example.rda", package = "cpam"))
#'
#' # Generate results table
#' res_example <- results(cpo_example)
#'
#' # plot all targets with changepoint at timepoint 0 and shape "ilin" (increasing linear)
#' plot_cluster(cpo_example, res_example, changepoints = 2, shapes = "ilin")
#'
plot_cluster <- function(cpo, res, changepoints, shapes, alpha = 0.1){

  if (alpha < 0 || alpha > 1) {
    stop("alpha must be between 0 and 1")
  }

  txs <-
    res %>%
    dplyr::filter(.data$cp %in% changepoints,
                  .data$shape %in% shapes) %>%
    dplyr::pull("target_id")

  if(length(txs) == 0) {
    warning("No targets found for the selected shapes and timepoints")
    return(NULL)
  }

  message(paste0("Plotting ",length(txs)," targets"))

  plot_data <-
    txs %>%
    purrr::set_names() %>%
    purrr::map(~ plot_cpam(cpo,target_id = .x,return_fits_only = TRUE) %>%
                     {if(rlang::is_na(.)) NULL else predict_lfc(.)}) %>%
    purrr::compact() %>%
    purrr::list_rbind(names_to = "target_id")

  plot_data %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$time, y = .data$lfc, group = .data$target_id)) +
    ggplot2::geom_line(alpha = alpha) +
    ggplot2::theme_classic()

}



predict_lfc <- function(fit, length.out = 200) {

  newdata <- dplyr::tibble(time = seq(min(fit$data$time), max(fit$data$time),
                                      length.out = length.out),
                           td = pmax(fit$cp, .data$time))

  dplyr::tibble(time = newdata$time,
                lfc = fit %>%
                  stats::predict(newdata = newdata) %>%
                  {.-.[1]} %>% {.*log(exp(1), base = 2)} %>% as.numeric,
  )
}



#' Plot fitted changepoint additive models
#'
#' @param cpo A cpam object containing count data, model fits, and optional changepoint/shape estimates
#' @param gene_id character; gene_id (mutually exclusive with target_id)
#' @param target_id character; target_id (mutually exclusive with gene_id)
#' @param cp_type One of "cp_1se" or "cp_min"; rule for selecting changepoint from fitted models.
#' See [estimate_changepoint()] for details.
#' @param shape_type One of "shape1" or "shape2"; which set of fitted shape patterns to use.
#' See [select_shape()] for details.
#' @param bs Shape pattern to fit ("null", "lin", "ilin", "dlin", or from cpo$bss).
#' Use "auto" (default) to use estimated shapes as per `shape_type`.
#' @param cp_fix Numeric; fixed changepoint time. Set to -1 (default) to use estimated changepoints
#' @param facet Logical; for multiple transcripts, plot in separate facets?
#' @param sp numerical; set the smooth parameter. NULL (default) for automatic selection
#' @param show_fit logical; show the fitted trend?
#' @param show_data logical; show (possibly normalized and scaled) data points?
#' @param show_fit_ci logical; show credible interval for the fitted trend?
#' @param show_data_ci logical; show bootstrapped quantile for data points?
#' @param ci_prob  "se" for standard error bands (see [mgcv::predict.gam()]), or numeric for simulation-based intervals.
#' If numerical, sets the probability for the simulation-based estimates of credible interval.
#' @param remove_null logical; only plot differentially expressed transcripts (not applicable for gene-only analyses)
#' @param null_threshold numeric; P value threshold for filtering out NULL transcripts
#' @param null_threshold_adj logical; use adjusted (default) or non-adjusted p-values for filtering targets
#' @param k_mult numerical; multiplier for the number of knots in the spline. Not recommended to change this value.
#' @param return_fits_only logical; return the model fits. Does not plot the function
#' @param family character; negative binomial ("nb", default) or Gaussian ("gaussian")
#' @param common_y_scale logical; for faceted plots of multiple transcripts, should the scale of the y-axis
#' be common or free.
#' @param scaled logical; scaled data by overdispersions (for bootstrapped data only)
#' @param base_size numeric; base font size for the plot
#'
#' @details
#' Plots the fitted trend and data points for a given gene or target. If a gene ID
#' is supplied, the function will plot all transcripts for that gene.
#' The function can also be used to return the model fit(s) only, which are
#' `gamObject` objects from the `mgcv` package.
#'
#'
#' @return a ggplot object
#' @export
#'
#' @examples
#' library(cpam)
#'
#' # load gene-only example cpam object
#' load(system.file("extdata", "cpo_example.rda", package = "cpam"))
#'
#' # example gene
#' plot_cpam(cpo_example, gene_id = "g003")
#'
#' # gene with estimated changepoint at timepoint 3
#' plot_cpam(cpo_example, gene_id = "g013")
#'
#' # manually set the changepoint
#' plot_cpam(cpo_example, gene_id = "g013", cp_fix = 2)
#'
plot_cpam <- function(cpo,
                         gene_id = NULL,
                         target_id = NULL,
                         cp_type = c("cp_1se","cp_min"),
                         shape_type = "shape1",
                         bs = "auto",
                         cp_fix = -1,
                         facet = F,
                         sp = NULL,
                         show_fit = TRUE,
                         show_data = TRUE,
                         show_fit_ci = TRUE,
                         show_data_ci = TRUE,
                         ci_prob = "se",
                         remove_null = FALSE,
                         null_threshold =  0.05,
                         null_threshold_adj = TRUE,
                         k_mult = 1.2,
                         #logged = F,
                         #gene_level_plot = F,
                         return_fits_only = FALSE,
                         family = "nb",
                         common_y_scale = TRUE,
                         scaled = FALSE,
                         base_size = 12){

  if(family != "nb"){
    warning("Plotting is currently only suppported for negative binomial models. Setting 'family ='nb''")
    family <- "nb"
  }

  cp_type <- match.arg(cp_type)
  shape_type <- match.arg(shape_type, c("shape1","shape2"))
  bs <- match.arg(bs, c("auto","null","lin","ilin","dlin",cpo$bss))
  if(!is.numeric(cp_fix)) stop("The fixed changepoint must be numeric")
  if(!cpo$bootstrap) show_data_ci <- FALSE
  if(!identical(ci_prob,"se")){
    if(is.numeric(ci_prob)){
      if(ci_prob < 0 | ci_prob > 1) stop("ci_prob must be between 0 and 1")
    } else {
      stop("ci_prob must be 'se' or a numeric value")
    }
  }
  if(cpo$gene_level) remove_null <- FALSE

  if(is.null(gene_id)){
    if(is.null(target_id)) stop("gene_id and target_id cannot both be null")

    data =
      cpo$data_long %>%
      dplyr::filter(target_id == {{target_id}}) %>%
      tidyr::nest(.by = "target_id")

    if(nrow(data) == 0) stop("Invalid target_id")

  } else{

    data =
      cpo$data_long %>%
      dplyr::filter(gene_id == {{gene_id}}) %>%
      tidyr::nest(.by = "target_id")

    if(nrow(data) == 0) stop("Invalid gene_id")
  }

  if(!is.null(cpo$changepoints)){
    data <- data %>%
      dplyr::left_join(cpo$changepoints %>% dplyr::select("target_id", cp = dplyr::all_of(cp_type)), by = "target_id")
    #cp_estimated <- any(is.na(data$cp))
  } else {
    data <- data %>% dplyr::mutate(cp = 0)
    #cp_estimated <- F
  }

  if(!is.null(cpo$shapes)){
    data <- data %>%
      dplyr::left_join(cpo$shapes %>% dplyr::select("target_id", shape = dplyr::all_of(shape_type)), by = "target_id")
    #shape_estimated <- any(is.na(data$shape))
  } else {
    data <- data %>% dplyr::mutate(shape = "tp")
    #shape_estimated <- F
  }

  if (remove_null) {
    if (!is.null(cpo$p_table)) {
      pval <- "q_val_target"
      if(!null_threshold_adj) pval <- "p_val_target"

      data <- data %>%
        dplyr::left_join(cpo$p_table %>% dplyr::select(dplyr::all_of(c("target_id",pval))), by = "target_id") %>%
        dplyr::filter(.data[[pval]] <= null_threshold |
                      .data$shape != "null")
    }
  }

  n_target = nrow(data)
  names_target = data %>% dplyr::pull("target_id")

  if(!is.null(target_id)){
    if(!target_id %in% names_target){
      target_id <- NULL
      stop("Invalid target_id supplied.")
    } else {
      data <- data %>% dplyr::filter(.data$target_id == {{target_id}})
      gene_id <- target_id
    }
  }

  if(bs != "auto"){
    data <- data %>% dplyr::mutate(shape = bs)
  } else {
    data <- data %>% dplyr::mutate(shape = tidyr::replace_na(.data$shape, "tp"))
  }

  if(cp_fix >= min(cpo$times) & cp_fix <= max(cpo$times)){
    data <- data %>% dplyr::mutate(cp = cp_fix)
  } else {
    data <- data %>% dplyr::mutate(cp = tidyr::replace_na(.data$cp,0))
  }


  if(n_target == 1) facet <- FALSE
  if(cpo$model_type == "case-control" & n_target != 1) facet <- TRUE

  if(scaled){
    obs <- data %>% dplyr::select(dplyr::all_of(c("target_id","data"))) %>% tidyr::unnest(cols = "data") %>%
      dplyr::mutate(dplyr::across(dplyr::starts_with("q_"), ~ .x/(.data$overdispersions*.data$norm_factor)),
                    counts = .data$counts/.data$norm_factor)
  } else{
    obs <- data %>% dplyr::select(dplyr::all_of(c("target_id","data"))) %>% tidyr::unnest(cols = "data") %>%
      dplyr::mutate(dplyr::across(dplyr::starts_with("q_"), ~ .x/(.data$norm_factor)),
                    counts = .data$counts_raw/.data$norm_factor)
  }

  fits <-
    data %>%
    dplyr::rowwise() %>%
    dplyr::mutate(fit = list(cpgam(data = data,
                                   family = family,
                                   model_type = cpo$model_type,
                                   regularize = cpo$regularize,
                                   bs = .data$shape,
                                   cp = .data$cp,
                                   sp = sp,
                                   k_mult = k_mult))) %>%
    dplyr::filter(!is.logical(.data$fit)) %>%
    dplyr::mutate(pred = list(predict_cpgam(fit = .data$fit, ci_prob = ci_prob, scaled = scaled)))

  if(return_fits_only){
    return(switch(as.character(length(fits$fit)),
                  "0" = NA,
                  "1" = fits$fit[[1]],
                  fits$fit %>% purrr::set_names(fits[["target_id"]])))
  }

  preds <-
    fits %>%
    dplyr::select(dplyr::all_of(c("target_id","cp","pred"))) %>%
    tidyr::unnest(cols = c("pred"))

  # facet labels
  bs_labels <- paste0("(",data$cp,", ",data$shape,")")
  facet_labels <- paste0(data$target_id,"\n (",data$cp,", ",data$shape,")")
  names(facet_labels) <- data$target_id

  col_values <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8,"Dark2"))(n_target)
  names(col_values) <- names_target

  gg <- preds %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$time)) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::scale_color_manual(values = col_values, aesthetics = c("colour","fill")) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(strip.background = ggplot2::element_rect(linewidth = 0),
                   #strip.text = ggplot2::element_text(size = 12),
                   plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
                   ) +
    ggplot2::scale_x_continuous(breaks = cpo$times) +
    ggplot2::labs(title = gene_id, y = "counts")


  if(cpo$model_type == "case-control") gg <- gg + ggplot2::aes(group = factor(.data$case), lty = factor(.data$case)) + ggplot2::scale_linetype_manual(values = c("dashed","solid"))
  if(show_fit) gg <- gg + ggplot2::geom_line(ggplot2::aes(y = .data$counts, col = .data$target_id), linewidth = base_size/16)
  if(show_fit_ci) gg <- gg + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$q_lo, ymax = .data$q_hi, fill = target_id), alpha = 0.1)
  if(show_data) gg <- gg + ggplot2::geom_point(ggplot2::aes(y = .data$counts, col = .data$target_id), data = obs, size = base_size/7)
  if(show_data_ci) gg <-  gg + ggplot2::geom_linerange(ggplot2::aes(ymin = .data$q_lo, ymax = .data$q_hi, col = target_id),
                                                       data = obs,
                                                       alpha = 0.3)
  if(facet) gg <- gg + ggplot2::facet_wrap(~ target_id,
                                           scales = dplyr::if_else(common_y_scale, "fixed","free_y"),
                                           labeller = ggplot2::as_labeller(facet_labels))
  if(facet | n_target == 1) gg <- gg + ggplot2::theme(legend.position = "none")
  if(!facet) gg <- gg + ggplot2::labs(subtitle = paste0(bs_labels, collapse = ", "))

  gg
}



