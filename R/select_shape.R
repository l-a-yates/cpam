#' Use model selection to select a shape from a list of candidates
#'
#' @param cpo a cpam object
#' @param subset character vector; names of targets or genes (if `cpo$gene_level = T`)
#' for which changepoints will be estimated
#' @param sp numerical >= 0; supply a fixed smoothing parameter.
#' Note, the fixd value is applied to shape constrained bases only (i.e., not `bs = 'tp'`).
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
                         sp = NULL,
                         bss = c("micv","mdcx","cv","cx","lin"),
                         include_tp = T,
                         family = c("nb","gaussian"),
                         score = "gcv",
                         cp_type = c("cp_1se","cp_min")) {

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
    } %>%
    dplyr::rowwise() %>%
    dplyr::mutate(cp = .data$data$cp[1],
           k = sum(unique(.data$data$time)>= .data$cp)) %>%
    dplyr::ungroup()

  message(paste0("Estimating shapes for ", nrow(data_nest), " targets"))
  if(include_tp) bss <- c(bss,"tp")
  message(paste0("Candidate shapes are bs = ", paste0(bss, collapse = ", "),"."))

  regularize <- cpo$regularize
  fixed_effects <- cpo$fixed_effects
  model_type <- cpo$model_type

  if(family == "nb" & score == "aic"){
    score <- "aic_negbin"
  }

  #return(data_nest)

  shapes_gt2 <-
    data_nest %>%
    dplyr::filter(.data$k >= 2) %>%
    dplyr::mutate(x =
                    .data$data %>%
                    pbmcapply::pbmclapply(function(d) {
                        bss %>%
                        purrr::set_names() %>%
                        purrr::map(~ cpgam(
                          data = d,
                          family = family,
                          regularize = regularize,
                          model_type = model_type,
                          cp = d$cp[1],
                          sp = if(.x == "tp") NULL else sp,
                          bs = .x,
                          fixed_effects = fixed_effects
                        ))  %>%
                        {.[!is.na(.)]} %>%
                        keep_converged() %>%
                        shape_selector(score)
                    }, mc.cores = cpo$num_cores)) %>%
    dplyr::select(.data$target_id, .data$x) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(shape1 = .data$x$shape1,
           shape2 = .data$x$shape2,
           lfc = .data$x$lfc) %>%
    dplyr::select(-.data$x) %>%
    dplyr::ungroup()

  null_tibble <- dplyr::tibble(time = cpo$times,lfc = 0)

  cpo$shapes <-
    data_nest %>%
    dplyr::select(-.data$data) %>%
    dplyr::left_join(dplyr::select(shapes_gt2, -.data$lfc), by = "target_id") %>%
    dplyr::mutate(shape1 = dplyr::if_else(.data$k==1,"null",.data$shape1),
                  shape2 = dplyr::if_else(.data$k==1,"null",.data$shape2),
                  k = NULL)


  cpo$lfc <-
    data_nest %>%
    dplyr::select(-.data$data) %>%
    dplyr::left_join(dplyr::select(shapes_gt2, .data$target_id,.data$lfc), by = "target_id") %>%
    dplyr::mutate(lfc = dplyr::if_else(.data$k==1, list(null_tibble), .data$lfc)) %>%
    dplyr::select(-.data$cp,-.data$k) %>%
    tidyr::unnest(.data$lfc) %>%
    tidyr::pivot_wider(id_cols = .data$target_id, names_from = .data$time, values_from = .data$lfc)

  cpo
}

shape_selector <- function(fits, score){
  fits <- fits[names(fits)==purrr::map_chr(fits,"bs")]

  if(length(fits)==0){
    shape1 <- shape2 <- NA
  } else if(length(fits)==1){
    shape1 <- shape2 <- names(fits)
  } else {
    score_table <- fits %>% purrr::map_dfc(~ do.call(score, args = list(fit = .x)))
    edf <- fits %>% purrr::map_dbl(~ .x[["edf"]] %>% sum)
    shape1 <- shape2 <- shape_ose(score_table,edf)

    if(shape1 == "tp"){
      shape2 <- shape_ose(score_table %>% dplyr::select(-.data$tp),edf)
    }
  }

  if(!is.na(shape1)){
    lfc = extract_lfc(fits[[shape1]])
  } else lfc <- NA

  if(shape2 == "lin"){
    if(fits[["lin"]]$coefficients["td"]>0){
      shape2 <- "ilin"
    } else shape2 <- "dlin"
  }
  if(shape1 != "tp") shape1 <- shape2


  dplyr::tibble(shape1 = shape1, shape2 = shape2, lfc = list(lfc))
}

shape_ose <- function(score_table, edf){
  m.min <- score_table %>% colSums() %>% which.min() %>% names
  st <-
    score_table %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ .x - .data[[m.min]])) %>%
    purrr::map(~ c(score_diff = sum(.x), se_diff = sd(.x)*sqrt(length(.x)))) %>%
    dplyr::bind_rows(.id = "model") %>%
    dplyr::mutate(dim = edf[.data$model]) %>%
    dplyr::filter(.data$se_diff >= .data$score_diff)

  if(all(c("mdcx","lin") %in% st$model)){
    if(edf["mdcx"]<=2.001) st <- st %>% dplyr::filter(.data$model!="mdcx")}
  if(all(c("micv","lin") %in% st$model)){
    if(edf["micv"]<=2.001) st <- st %>% dplyr::filter(.data$model!="mdcx")}

  st %>%
    dplyr::filter(.data$dim == min(.data$dim)) %>%
    dplyr::arrange(.data$score_diff) %>%
    dplyr::slice(1) %>%
    dplyr::pull(.data$model)
}

keep_converged <- function(fits){
  if(length(fits)==0){
    return(list())}
  else{
    fits[purrr::map_lgl(fits, ~ if(inherits(.x,"scam")) {.x$conv} else {.x$converged})]
  }
}

extract_lfc <- function(fit) {

  newdata = fit$data %>% dplyr::select(.data$time,.data$td) %>% dplyr::distinct()

  dplyr::tibble(time = newdata$time,
                 lfc = fit %>%
    stats::predict(newdata = newdata) %>%
    {.-.[1]} %>% {.*log(exp(1), base = 2)} %>% as.numeric,
   )
}
