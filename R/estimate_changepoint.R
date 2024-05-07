#' Use model selection to estimate changepoints
#'
#' @param cpo a cpam object
#' @param cps vector of candidate changepoints. Defaults to the set of observed timepoints
#' @param degs_only logical; should changepoints only be estimated for differentially expressed genes
#' @param deg_threshold logical; the threshold value for DEGs (ignored of `degs_only = F`)
#' @param subset character vector; names of targets or genes (if `cpo$gene_level = T`)
#' for which changepoints will be estimated
#' @param sp numerical >= 0; supply a fixed smoothing parameter.
#' @param bss character vector; names of candidate spline bases (i.e., candidate shape types).
#' Default is thin plate ("tp") splines.
#' @param family character; negative binomial ("nb", default) or Gaussian ("gaussian")
#' @param score character; model selection score, either Generalised Cross Validation ("gcv") or
#' Akaike Information Criterion ("aic")
#'
#' @return a cpam object with the estimated changepoint table added to the slot "changepoints"
#' @export
#'
#' @examples 1+1
estimate_changepoint <- function(cpo,
                        cps = NULL,
                        degs_only = T,
                        deg_threshold = 0.05,
                        subset = NULL,
                        sp = NULL,
                        bss = "tp",
                        family = c("nb","gaussian"),
                        score = "aic") {

  family <- match.arg(family)

  if(is.null(subset)){
    if(degs_only){
      if(cpo$gene_level){
        subset <- cpo$p_table %>%
          dplyr::filter(.data$q_val_target <= deg_threshold) %>%
          dplyr::pull(.data$target_id)
      } else {
        subset <- cpo$p_table %>%
          dplyr::filter(.data$q_val_gene <= deg_threshold) %>%
          dplyr::pull(.data$target_id)
      }
    }
  }

  data_nest <- cpo$data_long %>%
    tidyr::nest(.by = .data$target_id, .key = "data") %>%
    {
      if (is.null(subset))
        .
      else
        dplyr::filter(.,.data$target_id %in% subset)
    }

  message(paste0("Estimating changepoints for ", nrow(data_nest), " targets"))
  if(is.null(cps)) cps <- cpo$times
  message(paste0("Candidate changepoints are t = ", paste0(cps, collapse = ", "),"."))
  if(nrow(data_nest)>3000){
    message("Warning: model fitting may take several minutes for this many targets.")
    message("Fitting time can be reduced by setting a lower threshold for DEGs ",
    "(using `deg_threshold`) and/or reducing the ",
    "number of candidate time points (set using `cps`)")
  }

  regularize <- cpo$regularize
  model_type <- cpo$model_type

  cpo$changepoints <-
    data_nest %>%
    dplyr::select(.data$target_id) %>%
    dplyr::mutate(x =
                    data_nest$data %>%
                    pbmcapply::pbmclapply(function(d) {
                      score_tables =
                        bss %>%
                        purrr::map( ~ calc_score_table(
                          data = d,
                          cps = cps,
                          sp = sp,
                          bs = .x,
                          family = family,
                          score = score,
                          model_type = model_type,
                          regularize = regularize
                        ))

                      bs.min = score_tables %>% purrr::map( ~ colSums(.x) %>% min) %>% which.min()

                      dplyr::tibble(
                        cp_min = score_tables[[bs.min]] %>% colSums() %>% which.min() %>% names %>% as.numeric,
                        cp_1se = score_tables[[bs.min]] %>% ose_rule(nse = 1) %>% as.numeric,
                        p_mvn = simulate_p_mvn(score_tables[[bs.min]]),
                        bs = bss[bs.min],
                        family = family,
                        score_table = list(score_tables[[bs.min]])
                      )
                    }, mc.cores = cpo$num_cores)) %>%
   tidyr::unnest(cols = "x")

  cpo$p_mvn <-
    cpo$changepoints %>%
    dplyr::select(.data$target_id,p_mvn_target = .data$p_mvn) %>%
    dplyr::left_join(cpo$p_table %>% dplyr::select(.data$target_id, .data$gene_id, .data$counts_mean),
                     by = "target_id") %>%
    dplyr::relocate(.data$target_id, .data$gene_id, .data$counts_mean) %>%
    dplyr::mutate(q_mvn_target = stats::p.adjust(.data$p_mvn_target, method = "BH")) %>%
    dplyr::group_by(.data$gene_id) %>%
    dplyr::mutate(p_mvn_gene = aggregation::lancaster(pmax(.data$p_mvn_target,10e-320),.data$counts_mean/sum(.data$counts_mean))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(q_mvn_gene = stats::p.adjust(.data$p_mvn_gene, method = "BH"))

    cpo
}


ose_rule <- function(tab, nse=1){
  m.min <- tab %>% colSums %>% which.min %>% names

  if(nse == 0) return(m.min)

  ose_tab <-
    tab %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ .x - .data[[m.min]])) %>%
    purrr::map(~ c(score_diff = sum(.x), se_diff = sd(.x)*sqrt(length(.x)))) %>%
    dplyr::bind_rows(.id = "model") %>%
    dplyr::mutate(dim = rank(abs(as.numeric(.data$model))))

  ose_tab %>%
    dplyr::filter(.data$score_diff <= nse*.data$se_diff) %>%
    dplyr::filter(.data$dim == max(.data$dim)) %>%
    dplyr::filter(.data$score_diff == min(.data$score_diff)) %>%
    dplyr::pull(.data$model)

}


simulate_p_mvn <- function(score_table, nsim = 1e4, reg = 1e-3){
  mvnfast::rmvn(nsim,
                mu = score_table %>% colSums,
                sigma = stats::cov(score_table)*nrow(score_table) + diag(reg,ncol(score_table))) %>%
    `colnames<-`(colnames(score_table)) %>%
    apply(1,function(x) names(which.min(x))) %>%
    factor(levels = colnames(score_table)) %>%
    table %>%
    {.[length(.)]/nsim} %>%
    unname
}


calc_score_table <- function(data,
                             cps,
                             sp,
                             bs,
                             family = c("nb","gaussian"),
                             score = "aic",
                             model_type,
                             regularize){

  family <- match.arg(family)

  if(family == "nb"){
    if(score!="aic") warning("The score has been set to `aic` for negative binomial models")
    score <- "aic_negbin"
  }

  cps %>%
    purrr::set_names() %>%
    purrr::map(~ cpgam(data = data,
                       family = family,
                       model_type = model_type,
                       cp = .x,
                       regularize = regularize,
                       sp = sp,
                       bs = bs)) %>%
    {.[!is.na(.)]} %>%
    purrr::map(~ do.call(score, args = list(fit = .x))) %>%
    dplyr::bind_cols()
}

aic <- function(fit){
  if(fit$family$family == "gaussian"){
    if(fit$scale.estimated){
      sd = sqrt(mean(stats::residuals(fit, type = "response")^2))
    } else {
      sd = sqrt(fit$sig2)
    }
    ll <- stats::dnorm(stats::residuals(fit, type = "response"),
                0,
                sd,
                log = T)*(fit$weights)
  }

  if(fit$family$family %>% stringr::str_starts("Negative")){
    if("gam" %in% class(fit)) theta <- fit$family$getTheta(T)
    if("scam" %in% class(fit)) theta <- fit$family$getTheta()
    ll <- dnbl(fit$y,
               theta = theta,
               mu = fit$fitted.values)*fit$prior.weights
  }
  as.numeric(-2*ll + 2*attributes(stats::logLik(fit))$df/length(ll))
}

aic_negbin <- function(fit){
  y <- fit$y
  #if("gam" %in% class(fit)) Theta <- fit$family$getTheta(T)
  #if("scam" %in% class(fit)) Theta <- fit$family$getTheta()
  Theta <- fit$family$getTheta()
  mu <- fit$fitted.values
  wt <- fit$prior.weights

  term <- (y + Theta) * log(mu + Theta) - y * log(mu) +
    lgamma(y + 1) - Theta * log(Theta) + lgamma(Theta) - lgamma(Theta + y)

  2 * term * wt + 2*sum(fit$edf)/length(term)

}


# computes pointwise generalised cross validation (gcv) for a mgcv::gam fit
gcv.gam <- function(fit){
    gcvi = (length(fit$residuals)*((sqrt(fit$weights)*fit$residuals)^2))/(length(fit$residuals)-sum(fit$hat))^2
}

gcv <- function(fit){
  dev = stats::residuals(fit, "deviance")^2
  nobs = length(fit$y)

  if(inherits(fit,"gam")){
    gamma = 1
    trA = fit$hat %>% sum
  }

  if(inherits(fit,"scam")){
    gamma = fit$gamma
    trA = fit$trA
  }
  as.numeric(dev * nobs/(nobs - gamma * trA)^2)
}

# log probabilty density for negative binomial with continuous positive response values
dnbl <- function(x,theta,mu){
  lgamma(x + theta) -
    lgamma(theta) -
    lgamma(x + 1) +
    theta*log(theta/(theta + mu)) +
    x*log(mu/(mu+theta))
}





