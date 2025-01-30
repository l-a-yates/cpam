#' Use model selection to select a shape for each target
#'
#' @param cpo a cpam object
#' @param subset character vector; names of targets or genes (if `cpo$gene_level = T`)
#' for which changepoints will be estimated
#' @param sp numerical >= 0; supply a fixed smoothing parameter. If NULL (default), the smoothing
#' parameter is estimated. Note, this fixed value is in any case applied only to
#' shape constrained bases (i.e., not `bs = 'tp'`).
#' @param bss character vector; names of candidate spline bases (i.e., candidate shape types).
#' @param family character; negative binomial ("nb", default) or Gaussian ("gaussian")
#' @param score character; model selection score, either Generalised Cross Validation ("gcv") or
#' Akaike Information Criterion ("aic")
#' @param cp_type character; if changepoints have been estimated using [estimate_changepoint()],
#' which selection rule should be used. See [estimate_changepoint()] for details.
#'
#' @details
#' The function selects the best shape from a list of candidate shapes for each target.
#' It is typically the last step in the analysis, called after p-values have
#' been estimated using [compute_p_values()] and changepoints have
#' been estimated using [estimate_changepoint()].
#'
#' Two shape selections are generated. The first selecting among linear,
#' convex and concave shape classes and their monotonic variants (or the shape
#' set given by bss), and the second selecting among the first options plus an
#' 'unconstrained' smooth. The inclusion of the `unconstrained' type provides
#' the flexibility to detect targets beyond simpler trends.

#' For computational reasons, as per the changepoint estimation,
#' shapes are selected only for those genes, or their isoforms,
#' identified as significant at the chosen FDR threshold. This is
#' overridden by providing a subset of target names to the `subset` argument.
#'
#' @return a cpam object with the selected shapes added to the slot "shapes"
#' @export
#'
#' @examples
#' \dontrun{
#'
#' library(cpam)
#' library(dplyr)
#'
#' # Example Experimental Design
#' exp_design <- tibble(sample = paste0("s",1:50),
#'                      time = rep(c(0:4,10),
#'                      path = paste0("path/",sample,"/abundance.h5"))
#'
#' # Example Transcript-to-Gene Mapping
#' t2g <- readr::read_csv("path/to/t2g.csv")
#'
#' # Prepare a cpam object
#' cpo <- prepare_cpam(
#'  exp_design = exp_design,
#'  t2g = t2g,
#'  import_type = "kallisto",
#'  num_cores = 5)
#'
#'  # compute p-values
#'  cpo <- compute_p_values(cpo)
#'
#'  # estimate changepoints
#'  cpo <- estimate_changepoint(cpo)
#'
#'  # estimate shapes
#'  cpo <- select_shape(cpo)
#'
#'  # Inspect the shapes
#'  cpo$shapes
#'  }

select_shape <- function(cpo,
                         subset = NULL,
                         sp = NULL,
                         bss = c("micv","mdcx","cv","cx","lin","tp","null"),
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
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(.data$cp))

  if(nrow(data_nest) == 0) stop("No targets to estimate shapes for")
  message(paste0("Estimating shapes for ", nrow(data_nest), " targets"))
  message(paste0("Candidate shapes are bs = ", paste0(bss, collapse = ", "),"."))

  regularize <- cpo$regularize
  fixed_effects <- cpo$fixed_effects
  model_type <- cpo$model_type
  cp_max <- max(cpo$times)

  if(family == "nb" & score == "aic"){
    score <- "aic_negbin"
  }


  shapes <-
    data_nest %>%
    dplyr::mutate(x =
                    .data$data %>%
                    pbmcapply::pbmclapply(function(d) {
                        bss %>%
                        {
                          if(d$cp[1] == cp_max){
                            c("null")
                          } else {
                            .
                          }
                        } %>%
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
                        {
                          if(all(is.na(.)))
                            NA
                          else
                            keep_converged(.) %>%
                            shape_selector(score)
                        }

                    }, mc.cores = cpo$num_cores)) %>%
    dplyr::select(.data$target_id, .data$x) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(shape1 = .data$x$shape1,
           shape2 = .data$x$shape2,
           lfc = .data$x$lfc,
           pred = .data$x$pred) %>%
    dplyr::select(-.data$x) %>%
    dplyr::ungroup()

  cpo$shapes <-
    data_nest %>%
    dplyr::select(-.data$data) %>%
    dplyr::left_join(dplyr::select(shapes, -.data$lfc, -.data$pred), by = "target_id") %>%
    dplyr::mutate(k = NULL)

  cpo$changepoints <-
    cpo$changepoints %>%
    dplyr::left_join(cpo$shapes %>% dplyr::select(.data$target_id,.data$shape1), by = "target_id") %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(cp_type), ~ dplyr::if_else(.data$shape1 == "null", cp_max, .data[[cp_type]]))) %>%
    dplyr::select(-.data$shape1)

  pivot_names <- c("time")
  if(model_type == "case-control") pivot_names <- c("time","case")

  cpo$lfc <-
    data_nest %>%
    dplyr::select(-.data$data) %>%
    dplyr::left_join(dplyr::select(shapes, .data$target_id,.data$lfc), by = "target_id") %>%
    dplyr::select(-.data$cp,-.data$k) %>%
    tidyr::unnest("lfc") %>%
    {
      if(model_type == "case-control"){
        dplyr::mutate(.,case = c("ctrl","case")[.data$case + 1])
      } else
        .
    } %>%
    tidyr::pivot_wider(id_cols = "target_id", names_from = pivot_names, values_from = "lfc")

  cpo$pred <-
    data_nest %>%
    dplyr::select(-.data$data) %>%
    dplyr::left_join(dplyr::select(shapes, .data$target_id,.data$pred), by = "target_id") %>%
    dplyr::select(-.data$cp,-.data$k) %>%
    tidyr::unnest("pred") %>%
    {
      if(model_type == "case-control"){
        dplyr::mutate(.,case = c("ctrl","case")[.data$case + 1])
      } else
        .
    } %>%
    tidyr::pivot_wider(id_cols = "target_id", names_from = pivot_names, values_from = "pred")

  cpo
}

shape_selector <- function(fits, score){
  fits <- fits[names(fits)==purrr::map_chr(fits,"bs")]

  #remove fits with very large predictive standard errors
  fits <- fits[map_dbl(fits,"se_ratio")<5]

  # remove cv if it is an micv shape
  if(all(c("micv","cv") %in% names(fits))){
    d = extract_lfc(fits[["cv"]])
    if(which.max(dplyr::pull(d,"lfc")) == nrow(d)) fits[["cv"]] <- NULL
    }
  # remove cx if it is an mdcx shape
  if(all(c("mdcx","cx") %in% names(fits))){
    d = extract_lfc(fits[["cx"]])
    if(which.min(dplyr::pull(d,"lfc")) == nrow(d)) fits[["cx"]] <- NULL
  }

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
    pred = extract_pred(fits[[shape1]])
  } else lfc <- pred <- NA

  if(shape2 == "lin"){
    if(fits[["lin"]]$coefficients["td"]>0){
      shape2 <- "ilin"
    } else shape2 <- "dlin"
  }
  if(shape1 != "tp") shape1 <- shape2


  dplyr::tibble(shape1 = shape1, shape2 = shape2, lfc = list(lfc), pred = list(pred))
}

shape_ose <- function(score_table, edf, tol = 0.01){
  m.min <- score_table %>% colSums() %>% which.min() %>% names
  st <-
    score_table %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ .x - .data[[m.min]])) %>%
    purrr::map(~ c(score_diff = sum(.x), se_diff = sd(.x)*sqrt(length(.x)))) %>%
    dplyr::bind_rows(.id = "model") %>%
    dplyr::mutate(edf = edf[.data$model]) %>%
    dplyr::filter(.data$se_diff >= .data$score_diff)

  # remove shapes that are effectively of linear complexity if lin is present
  if("lin" %in% st$model){
    edf_lin <- st %>% dplyr::filter(.data$model == "lin") %>% dplyr::pull(.data$edf)
    st <- st %>% dplyr::filter(!(.data$model %in% c("mdcx","micv","cv","cx") &
                                   .data$edf - edf_lin[1] <= tol))
  }
  # remove shapes that are effectively of null complexity if null is present
  if("null" %in% st$model){
    edf_null <- st %>% dplyr::filter(.data$model == "null") %>% dplyr::pull(.data$edf)
    st <- st %>% dplyr::filter(!(.data$model %in% c("mdcx","micv","cv","cx") &
                                   .data$edf - edf_null[1] <= tol))
  }

  st %>%
    dplyr::filter(.data$edf == min(.data$edf)) %>%
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
  svars <- c("time","td")
  if(fit$model_type == "case-control") svars <- c("time","td","case")
  newdata = fit$data %>% dplyr::select(dplyr::all_of(svars)) %>% dplyr::distinct()

  newdata %>%
    dplyr::select(-.data$td) %>%
    dplyr::mutate(lfc = fit %>%
                    stats::predict(newdata = newdata) %>%
                    {. - .[1]} %>%
                    {. * log(exp(1), base = 2)} %>%
                    as.numeric)
}

extract_pred <- function(fit, scaled = F) {
  svars <- c("time","td")
  if(fit$model_type == "case-control") svars <- c("time","td","case")
  newdata = fit$data %>% dplyr::select(dplyr::all_of(svars)) %>% dplyr::distinct()

  od <- 1
  if(!is.null(fit$data$overdispersions) & !scaled){
    od <- as.numeric(fit$data$overdispersions[1])
  }

  newdata %>%
    dplyr::select(-.data$td) %>%
    dplyr::mutate(pred = fit %>%
                  stats::predict(newdata = newdata, type = "response") %>%
                  as.numeric) %>%
    dplyr::mutate(pred = .data$pred*od)
}


