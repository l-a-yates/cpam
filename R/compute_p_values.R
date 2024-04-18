compute_p_values <- function(cpo,
                             subset = NULL,
                             p_adj_method = "BH",
                             fixed_effects = NULL,
                             gam_method = "REML",
                             gam_optimizer = "efs"){

  if(!is.null(subset)){
    if(!is.character(subset)){
      stop("subset must be character vector of target_id values")
    } else {
      if(!all(subset %in% rownames(cpo$count_matrix_raw))) warning("Subset contains invalid targets")
      if(!all(subset %in% cpo$target_to_keep)) warning("Subset contains targets that have been filtered out")
    }
  }

  regularize <- cpo$regularize

  if(is.null(fixed_effects)){
    fixed_effects <- ""
  } else {
    fixed_effects <- paste(fixed_effects,"+")
  }

  f_string <- paste0("counts ~",
                     fixed_effects,
                     " s(time, bs = 'tp', k = ",
                     length(unique(cpo$exp_design$time)),
                     ")")

  if(cpo$model_type == "case-control"){
    f_string <- paste0(f_string,
                       " + s(td, bs = 'tp', k = ",
                       length(unique(cpo$exp_design$time)),
                       ")" )
    test_var <- "s(td)"
  } else test_var <- "s(time)"

  data_nest <- cpo$data_long %>%
    {
      if(cpo$model_type == "case-control"){
        dplyr::mutate(.,td = .data$time*.data$case)
      } else{
        .
      }
    } %>%
    tidyr::nest(.by = .data$target_id, .key = "data") %>%
    {
      if (is.null(subset))
        .
      else
        dplyr::filter(.,.data$target_id %in% subset)
    }

  p_table =
    data_nest %>%
    dplyr::rowwise() %>%
    dplyr::transmute(.data$target_id, counts_mean = mean(.data$data$counts)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      p_val_target = data_nest$data %>%
        pbmcapply::pbmclapply(function(d) {
          p.val = try(mgcv::gam(
            stats::as.formula(f_string),
            data = d,
            method = gam_method,
            optimizer = gam_optimizer,
            offset = log(d$norm_factor),
            family = mgcv::nb(theta = if (regularize)
              as.numeric(1/d$disp[1])
              else
                NULL)
          ) %>%
            summary %>%
            {
              .$s.table[test_var, "p-value"]
            })
          if ("try-error" %in% class(p.val)) {
            return(NA)
          }
          p.val
        }, mc.cores = cpo$num_cores) %>% purrr::list_c()
    ) %>%
    dplyr::mutate(q_val_target = stats::p.adjust(.data$p_val_target, method = p_adj_method))

  if(cpo$aggregate_to_gene){
    p_table <-
      p_table %>%
      dplyr::left_join(cpo$t2g, by = c("target_id")) %>%
      dplyr::group_by(.data$gene_id) %>%
      dplyr::mutate(p_val_gene = lancaster(pmax(.data$p_val_target,10e-320),.data$counts_mean/sum(.data$counts_mean))) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(q_val_gene = stats::p.adjust(.data$p_val_gene, method = p_adj_method)) %>%
      dplyr::relocate(.data$target_id, .data$gene_id, .data$counts_mean)
  }

  cpo$p_table <- p_table

  cpo

}