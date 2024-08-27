#' Compute P values
#'
#' @param cpo a cpam object
#' @param subset a character vector of target_id names
#' @param p_adj_method method for P value adjustment
#' @param fixed_effects a model formula (RHS only) to provide or update
#' fixed affects. Terms must correspond to columns in `exp_design`
#' @param gam_method fitting method for `mgcv::gam`
#' @param gam_optimizer optimization method for `mgcv::gam`
#' @param silent logical; silences warnings from model fitting
#'
#' @return an updated cpam object with P values stored in the new slot "p_table"
#' @export
#'
#' @examples 1 + 1
compute_p_values <- function(cpo,
                             subset = NULL,
                             p_adj_method = "BH",
                             fixed_effects = NULL,
                             gam_method = "REML",
                             gam_optimizer = "efs",
                             silent = TRUE){

  if(!is.null(subset)){
    if(!is.character(subset)){
      stop("subset must be character vector of target_id values")
    } else {
      if(!all(subset %in% rownames(cpo$count_matrix_raw))) warning("Subset contains invalid targets")
      if(!all(subset %in% cpo$target_to_keep)) warning("Subset contains targets that have been filtered out")
    }
  }

  regularize <- cpo$regularize


  if (is.null(cpo$fixed_effects)) {
    if (is.null(fixed_effects)) {
      fe_string <- ""
    } else {
      message("Setting fixed effects")
      cpo$fixed_effects <- deparse(fixed_effects[[2]])
      fe_string <- paste(cpo$fixed_effects,"+")
    }
  } else {
    if(!is.null(fixed_effects)){
      message("Updating fixed effects")
      cpo$fixed_effects <- deparse(fixed_effects[[2]])
    }
    fe_string <- paste(cpo$fixed_effects,"+")
  }

  f_string <- paste0("counts ~",
                     fe_string,
                     " s(time, bs = 'tp', k = ",
                     length(unique(cpo$exp_design$time)),
                     ")")

  data_nest <- cpo$data_long %>%
    {
      if(cpo$model_type == "case-control"){
        dplyr::mutate(.,td = .data$time*(.data$condition=="treatment"))
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


  if(cpo$model_type == "case-only"){
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
              .$s.table["s(time)", "p-value"]
            }, silent = silent)
          if ("try-error" %in% class(p.val)) {
            return(NA)
          }
          p.val
        }, mc.cores = cpo$num_cores) %>% purrr::list_c()
    ) %>%
    dplyr::mutate(p_val_target = pmax(.data$p_val_target, 10e-320),
                  q_val_target = stats::p.adjust(.data$p_val_target, method = p_adj_method))
  } # end case-only


  if(cpo$model_type == "case-control"){

    f_string_cc <- paste0("counts ~",
                          fe_string,
                          " s(time, bs = 'tp', k = ",
                          length(unique(cpo$exp_design$time)),
                          ") + s(td, bs = 'tp', k = ",
                          length(unique(cpo$exp_design$time)),
                          ")")
    p_table =
      data_nest %>%
      dplyr::rowwise() %>%
      dplyr::transmute(.data$target_id, counts_mean = mean(.data$data$counts)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        p_val_target = data_nest$data %>%
          pbmcapply::pbmclapply(function(d) {
            p.val = try(mgcv::anova.gam(
              mgcv::gam(
                formula = stats::as.formula(f_string_cc),
                data = d,
                method = gam_method,
                optimizer = gam_optimizer,
                offset = log(d$norm_factor),
                family = mgcv::nb(theta = if (regularize)
                  as.numeric(1 / d$disp[1])
                  else
                    NULL)
              ),
              mgcv::gam(
                formula = stats::as.formula(f_string),
                data = d,
                method = gam_method,
                optimizer = gam_optimizer,
                offset = log(d$norm_factor),
                family = mgcv::nb(theta = if (regularize)
                  as.numeric(1 / d$disp[1])
                  else
                    NULL)
              ),
              test = "F"
            )$`Pr(>Chi)`[2]
            ,
            silent = silent)
            if ("try-error" %in% class(p.val)) {
              return(NA)
            }
            p.val
          }, mc.cores = cpo$num_cores) %>% purrr::list_c()
      ) %>%
      dplyr::mutate(p_val_target = pmax(.data$p_val_target, 10e-320),
                    q_val_target = stats::p.adjust(.data$p_val_target, method = p_adj_method))
  } # end case-control


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
