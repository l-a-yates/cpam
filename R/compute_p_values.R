#' Compute p-values for each target ID
#'
#' @param cpo a cpam object
#' @param subset a character vector of target_id names
#' @param p_adj_method method for p-value adjustment
#' @param aggregation_method method for aggregating target-level p-values to gene-level;
#'   either "lancaster" (default) or "acat"
#' @param gam_method fitting method for `mgcv::gam` (default is "REML")
#' @param gam_optimizer optimization method for `mgcv::gam` (default is "efs")
#' @param silent logical; silences warnings from model fitting (default is TRUE)
#'
#' @details
#' This function computes p-values for each target_id in the supplied cpam object.
#' The p-values are computed from a negative binomial GAM model
#' with a thin-plate spline basis function(s) for time
#' using the `mgcv` package.
#'
#' The p-values are stored in the new slot `p_table` in the cpam object.
#' If `aggregate_to_gene` is set to `TRUE` (default),
#' the target p-values are aggregated to the gene level using the `lancaster` method.
#' The columns `p_val_target` and `p_val_gene` store the raw p-values for target- and gene-level, respectively.
#' The function also computes adjusted p-values using the `p_adj_method`.
#' The default method is "BH" (Benjamini-Hochberg),
#' but any methods supported by the function `p.adjust` can be used.
#' The adjusted p-values are stored in the columns `q_val_target` and `q_val_gene`.
#'
#' @return an updated cpam object with raw, adjusted, and possibly aggregated p-values stored in the new slot "p_table"
#' @export
#'
#' @references
#' Wood, S.N. (2013a) On p-values for smooth components of an extended
#' generalized additive model. Biometrika 100:221-228 doi:10.1093/biomet/ass048
#'
#' Yi L, Pachter L (2018). aggregation: p-Value Aggregation Methods. R package version 1.0.1,
#' https://CRAN.R-project.org/package=aggregation.
#'
#' @examples
#' library(cpam)
#'
#' # load gene-only example cpam object
#' load(system.file("extdata", "cpo_example.rda", package = "cpam"))
#'
#' # run on a small subset of the example data
#' cpo <- compute_p_values(cpo_example, subset = paste0("g00",1:9))
#' cpo$p_table
#'

compute_p_values <- function(cpo,
                             subset = NULL,
                             p_adj_method = "BH",
                             aggregation_method = "lancaster",
                             gam_method = "REML",
                             gam_optimizer = "efs",
                             silent = TRUE){

  if(!is.null(subset)){
    if(!is.character(subset)){
      stop("subset must be character vector of target_id values")
    } else {
      if(!all(subset %in% rownames(cpo$count_matrix_raw))){
        stop("Subset contains invalid targets")}
      else if(!all(subset %in% cpo$target_to_keep)) {
        stop("Subset contains targets that have been filtered out")
      }
    }
  }

  aggregation_method <- match.arg(aggregation_method, choices = c("lancaster", "acat"))

  regularize <- cpo$regularize

  if(!is.null(cpo$fixed_effects)){
    fe_string <- paste(cpo$fixed_effects,"+")
  } else {
    fe_string <- ""
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
    tidyr::nest(.by = "target_id", .key = "data") %>%
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
                          cpo$intercept_cc,
                          " + ",
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

    aggregation_fn <- match.fun(aggregation_method)

    p_table <-
      p_table %>%
      dplyr::left_join(cpo$t2g, by = c("target_id")) %>%
      dplyr::group_by(.data$gene_id) %>%
      dplyr::mutate(p_val_gene = aggregation_fn(pmax(.data$p_val_target,10e-320),.data$counts_mean/sum(.data$counts_mean))) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(q_val_gene = stats::p.adjust(.data$p_val_gene, method = p_adj_method)) %>%
      dplyr::relocate("target_id", "gene_id", "counts_mean")
  }

  cpo$p_table <- p_table
  cpo$p_adj_method <- p_adj_method
  cpo

}
