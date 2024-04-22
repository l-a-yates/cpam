#' Use model selection to select a shape from a list of candidates
#'
#' @param cpo a cpam object
#' @param subset character vector; names of targets or genes (if `cpo$gene_level = T`)
#' for which changepoints will be estimated
#' @param sp numerical >= 0; supply a fixed smoothing parameter.
#' @param bss character vector; names of candidate spline bases (i.e., candidate shape types).
#' @param include_tp logical; should the non-shape-constrained basis "tp" be included as a candidate
#' @param family character; negative binomial ("nb", default) or Gaussian ("gaussian")
#' @param score character; model selection score, either Generalised Cross Validation ("gcv") or
#' Akaike Information Criterion ("aic")
#' @param cp_type character; if changepoints have been estimated using [estimate_changepoint()],
#' which selection rule should be used. See [estimate_changepoint()] for details.
#'
#' @return a cpam object with the selected shapes added to the slot "shapes"
#' @export
#'
#' @examples
select_shape <- function(cpo,
                         subset = NULL,
                         sp = 0.001,
                         bss = c("micv","mdcx","cv","cx"),
                         include_tp = T,
                         family = c("nb","gaussian"),
                         score = "aic",
                         cp_type = c("cp_min","cp_1se")) {


  cp_type <- match.arg(cp_type)
  family <- match.arg(family)

  if(is.null(cpo$changepoints)) stop("Changepoints must be estimated prior to shape selection")

  data_nest <-
    cpo$data_long %>%
    dplyr::right_join(cpo$changepoints %>%
                       dplyr::select(.data$target_id, cp = dplyr::all_of(cp_type)),
                     by = "target_id") %>%
    dplyr::mutate(cp = as.numeric(.data$cp)) %>%
    tidyr::nest(.by = .data$target_id, .key = "data") %>%
    {
      if (is.null(subset))
        .
      else
        dplyr::filter(.,.data$target_id %in% subset)
    }

  message(paste0("Estimating shapes for ", nrow(data_nest), " targets"))
  if(include_tp) bss <- c(bss,"tp")
  message(paste0("Candidate shapes are bs = ", paste0(bss, collapse = ", "),"."))

  regularize <- cpo$regularize
  fixed_effects <- cpo$fixed_effects
  model_type <- cpo$model_type

  if(family == "nb"){
    if(score!="aic") warning("The score has been set to `aic` for negative binomial models")
    score <- "aic_negbin"
  }

  cpo$shapes <-
    data_nest %>%
    dplyr::select(.data$target_id) %>%
    dplyr::mutate(shape =
                    data_nest$data %>%
                    pbmcapply::pbmclapply(function(d) {
                        bss %>%
                        purrr::set_names() %>%
                        purrr::map(~ cpgam(
                          data = d,
                          family = family,
                          regularize = regularize,
                          model_type = model_type,
                          cp = d$cp[1],
                          sp = sp,
                          bs = .x,
                          fixed_effects = fixed_effects
                        ))  %>%
                        {.[!is.na(.)]} %>%
                        purrr::map_dbl(~ do.call(score, args = list(fit = .x)) %>% sum) %>%
                        {c(shape1 = names(which.min(.)),
                           shape2 = names(which.min(.[names(.) != "tp"])),
                           family = family)}
                    }, mc.cores = cpo$num_cores) %>%
                    dplyr::bind_rows()) %>%
    tidyr::unnest(cols = "shape")
  cpo
}
