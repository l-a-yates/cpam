#' Plot a fitted cpgam
#'
#' @param cpo a cpam object
#' @param gene_id character; gene_id
#' @param target_id character; target_id
#' @param cp_type character; if changepoints have been estimated using [estimate_changepoint()],
#' which selection rule should be used. See [estimate_changepoint()] for details.
#' @param shape_type character; if shapes have been estimated using [select_shape()],
#' which set of candidate shapes should be used. See [select_shape()] for details.
#' @param bs character; set the basis (i.e. shape)
#' @param cp_fix numerical; set the changepoint
#' @param facet logical; for genes with multiple transcripts, should the transcripts be plotted in separate facets
#' @param sp numerical; set the smooth parameter
#' @param show_fit logical; show the fitted trend
#' @param show_data logical; show (possibly normalized and scaled) data points
#' @param show_fit_ci logical; show credible interval for the fitted trend
#' @param show_data_ci logical; show bootstrapped quantile for data points
#' @param ci_prob numerical; set the probabilty for the credible interval of the fitted trend
#' @param remove_null logical; only plot differentially expressed transcripts
#' @param null_threshold numeric; P value threshold for filtering out NULL transcripts
#' @param gene_level_plot logical; plot gene-level data and fitted trend
#' @param family character; negative binomial ("nb", default) or Gaussian ("gaussian")
#'
#' @return a ggplot object
#' @export
#'
#' @examples 1+1
plot_gene_co <- function(cpo,
                         gene_id = NULL,
                         target_id = NULL,
                         cp_type = c("cp_1se","cp_min"),
                         shape_type = "shape1",
                         bs = "auto",
                         cp_fix = -1,
                         facet = F,
                         sp = NULL,
                         show_fit = T,
                         show_data = T,
                         show_fit_ci = T,
                         show_data_ci = T,
                         ci_prob = 2*pnorm(1)-1,
                         remove_null = F,
                         null_threshold =  0.05,
                         #logged = F,
                         gene_level_plot = F,
                         family = "nb"){

  if(family != "nb"){
    warning("Plotting is only suppported for negative binomial models. Setting 'family ='nb''")
    family <- "nb"
  }
  if(gene_level_plot) stop("Watch this space! Gene-level plotting is coming your way very soon")
  cp_type <- match.arg(cp_type)
  shape_type <- match.arg(shape_type, c("shape1","shape2"))
  bs <- match.arg(bs, c("auto","null",cpo$bss))
  if(!is.numeric(cp_fix)) stop("The fixed changepoint must be numeric")
  if(!cpo$bootstrap) show_data_ci <- F

  if(is.null(gene_id) & is.null(target_id)) stop("gene_id and target_id cannon both be null")

  data =
    cpo$data_long %>%
    dplyr::filter(gene_id == {{gene_id}}) %>%
    tidyr::nest(.by = "target_id")

  if(nrow(data) == 0) stop("Invalid gene_id")

  if(!is.null(cpo$changepoints)){
    data <- data %>%
      dplyr::left_join(cpo$changepoints %>% dplyr::select(.data$target_id, cp = dplyr::all_of(cp_type)), by = "target_id")
    cp_estimated <- any(is.na(data$cp))
  } else {
    data <- data %>% dplyr::mutate(cp = 0)
    cp_estimated <- F
  }

  if(!is.null(cpo$shapes)){
    data <- data %>%
      dplyr::left_join(cpo$shapes %>% dplyr::select(.data$target_id, shape = dplyr::all_of(shape_type)), by = "target_id")
    shape_estimated <- any(is.na(data$shape))
  } else {
    data <- data %>% dplyr::mutate(shape = "tp")
    shape_estimated <- F
  }

  if (remove_null & !gene_level_plot) {
    if (!is.null(cpo$p_table)) {
      data <- data %>%
        dplyr::left_join(cpo$p_table %>% dplyr::select(.data$target_id, .data$q_val_target), by = "target_id") %>%
        dplyr::filter(.data$q_val_target <= null_threshold)
    } #else{
    #data <- data %>% dplyr::filter(.data$cp != 240)
    #  }
  }

  n_target = nrow(data)
  names_target = data %>% dplyr::pull(.data$target_id)

  if(!is.null(target_id)){
    if(!target_id %in% names_target){
      target_id <- NULL
      stop("Invalid target_id supplied.")
    } else {
      data <- data %>% dplyr::filter(target_id == {{target_id}})
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


  if(n_target == 1) facet <- F

  obs <- data %>% dplyr::select(.data$target_id,.data$data) %>% tidyr::unnest(cols = "data") %>%
    dplyr::mutate(dplyr::across(dplyr::starts_with("q_"), ~ .x/(.data$overdispersions*.data$norm_factor)),
                  counts = .data$counts/.data$norm_factor)

  fits <-
    data %>%
    dplyr::rowwise() %>%
    dplyr::mutate(fit = list(cpgam(data = data,
                                   family = family,
                                   model_type = cpo$model_type,
                                   regularize = cpo$regularize,
                                   bs = .data$shape,
                                   cp = .data$cp,
                                   sp = sp))) %>%
    dplyr::filter(!is.logical(.data$fit)) %>%
    dplyr::mutate(pred = list(predict_cpgam(fit = .data$fit, ci_prob = ci_prob)))

  preds <-
    fits %>%
    dplyr::select(.data$target_id, .data$cp, .data$pred) %>%
    tidyr::unnest(cols = c("pred"))
  #dplyr::mutate(bcp = time<cp)

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
    ggplot2::theme_classic() +
    ggplot2::theme(strip.background = ggplot2::element_rect(linewidth = 0),
                   strip.text = ggplot2::element_text(size = 12),
                   plot.title = ggplot2::element_text(hjust = 0.5, size = 16, face = "bold")) +
    ggplot2::scale_x_continuous(breaks = cpo$times) +
    ggplot2::labs(title = gene_id, y = "counts")

  if(show_fit) gg <- gg + ggplot2::geom_line(ggplot2::aes(y = .data$counts, col = .data$target_id))
  if(show_fit_ci) gg <- gg + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$q_lo, ymax = .data$q_hi, fill = target_id), alpha = 0.1)
  if(show_data) gg <- gg + ggplot2::geom_point(ggplot2::aes(y = .data$counts, col = .data$target_id), data = obs)
  if(show_data_ci) gg <-  gg + ggplot2::geom_linerange(ggplot2::aes(ymin = .data$q_lo, ymax = .data$q_hi, col = target_id),
                                                       data = obs,
                                                       alpha = 0.3)
  if(facet) gg <- gg + ggplot2::facet_wrap(~ target_id, labeller = ggplot2::as_labeller(facet_labels))
  if(facet) gg <- gg + ggplot2::theme(legend.position = "none")
  if(!facet) gg <- gg + ggplot2::labs(subtitle = paste0(bs_labels, collapse = ", "))

  gg
}


