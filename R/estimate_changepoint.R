#' Use model selection to estimate changepoints
#'
#' @param cpo a cpam object
#' @param cps vector of candidate changepoints. Defaults to the set of observed timepoints
#' @param degs_only logical; should changepoints only be estimated for differentially expressed genes
#' @param deg_threshold logical; the threshold value for DEGs (ignored of `degs_only = F`)
#' @param subset character vector; names of targets or genes (if `cpo$gene_level = T`)
#'  for which changepoints will be estimated
#' @param sp numerical >= 0; supply a fixed smoothing parameter.
#'  This can decrease the fitting time but it is not recommended as changepoints estimation
#'  is sensitive to smoothness.
#' @param bss character vector; names of candidate spline bases (i.e., candidate shape types).
#'  Default is thin plate ("tp") splines.
#' @param family character; negative binomial ("nb", default) or Gaussian ("gaussian", not currently supported)
#' @param score character; model selection score, either Generalised Cross Validation ("gcv") or
#'  Akaike Information Criterion ("aic")
#' @param compute_mvn Use simulation to compute p-value under multivariate normal model of the
#'  model scores
#'
#' @details
#' This function estimates changepoints for each target_id. The assumed
#' trajectory type for this modelling stage is initially constant followed by
#' a changepoint into thin-plate smoothing spline.
#'
#' By default, candidate time points are limited to the discrete observed values
#' in the series, since, despite the use of smoothing constraints,
#' there is generally insufficient information to infer the timing of
#' changepoints beyond the temporal resolution of the data. In any case, the
#' candidate points can be set manually using the `cps` argument.
#'
#' To estimate changepoints, a model is fit for each candidate changepoint and
#' generalised cross-validation (GCV, default) or the Akaike Information
#' Criterion (AIC) are used to select among them. Model-selection uncertainty
#' is dealt with by computing the one-standard-error rule, which identifies the
#' least complex model within one standard error of the best scoring model.
#'
#' Both the minimum and the one-standard-error (default) models are stored in the returned
#' slot "changepoints" so that either can be used. In addition to these, this function also computes the
#' probability (denoted `p_mvn`) that the null model is the best scoring model, using a simulation
#' based approach based on the multivariate normal model of the pointwise
#' model scores.
#'
#' Given the computational cost of fitting a separate model for each candidate
#' changepoint, cpam only estimates changepoints for targets associated with
#' 'significant' genes at the chosen threshold `deg_threshold`.
#'
#'
#' @return a cpam object with the estimated changepoint table added to the slot "changepoints"
#' @export
#'
#' @examples
#' library(cpam)
#' library(dplyr)
#'
#' # Using a small subset of the example data
#' cpo <- prepare_cpam(exp_design = exp_design_example,
#'                     count_matrix = count_matrix_example[1:20,],
#'                     gene_level = TRUE,
#'                     num_cores = 1)
#' cpo <- compute_p_values(cpo)
#' cpo <- estimate_changepoint(cpo)
#' cpo$changepoints
#'
#' \dontrun{
#'
#' # Example Experimental Design
#' exp_design <- tibble(sample = paste0("s",1:50),
#'                      time = rep(c(0:4), each = 10),
#'                      path = paste0("path/",sample,"/abundance.h5"))
#'
#' # Example Transcript-to-Gene Mapping
#' t2g <- readr::read_csv("path/to/t2g.csv")
#'
#' # Fit model
#' cpo <- prepare_cpam(
#'  exp_design = exp_design,
#'  t2g = t2g,
#'  import_type = "kallisto",
#'  num_cores = 5)
#' cpo <- compute_p_values(cpo)
#' cpo <- estimate_changepoint(cpo)
#' cpo$changepoints
#' }
#' @references
#' Yates, L. A., S. A. Richards, and B. W. Brook. 2021.
#' Parsimonious model selection using information theory:
#' a modified selection rule.
#' Ecology 102(10):e03475. 10.1002/ecy.3475
#'
estimate_changepoint <- function(cpo,
                        cps = NULL,
                        degs_only = T,
                        deg_threshold = 0.05,
                        subset = NULL,
                        sp = NULL,
                        bss = "tp",
                        family = c("nb","gaussian"),
                        score = "aic",
                        compute_mvn = T) {

  family <- match.arg(family)

  if(is.null(subset)){
    if(degs_only){
      if(cpo$gene_level){
        subset <- cpo$p_table %>%
          dplyr::filter(.data$q_val_target <= deg_threshold | is.na(.data$q_val_target)) %>%
          dplyr::pull(.data$target_id)
      } else {
        subset <- cpo$p_table %>%
          dplyr::filter(.data$q_val_gene <= deg_threshold | is.na(.data$q_val_target)) %>%
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

  if(nrow(data_nest) == 0) {
    cli::cli_abort("No targets selected. Please review subset and/or p-value thresholds.")
    return(cpo)
  }

  cli::cli_text("Estimating changepoints for {.val {nrow(data_nest)}} targets")
  if(is.null(cps)) cps <- cpo$times
  cli::cli_text("Candidate changepoints are t = {cps}")
  if(nrow(data_nest)>3000){
    cli::cli_alert_warning("Warning: model fitting may take several minutes for this many targets.")
    cli::cli_alert_info("Fitting time can be reduced by:")
    cli::cli_ul()
    cli::cli_li("setting a lower threshold for DEGs (using `deg_threshold`)")
    cli::cli_li("subsetting the list of targets (using `subset`)")
    cli::cli_li("reducing the number of candidate time points (using `cps`)")
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

                      if(all(is.na(score_tables))){
                        bs.min <- cp_min <- cp_1se <- p_mvn <- bs <- score_table <- NA
                      } else {
                        bs.min = score_tables %>% purrr::map( ~ colSums(.x) %>% min) %>% which.min()
                        cp_min = score_tables[[bs.min]] %>% colSums() %>% which.min() %>% names %>% as.numeric
                        cp_1se = score_tables[[bs.min]] %>% ose_rule(nse = 1) %>% as.numeric
                        if(compute_mvn){
                          p_mvn = simulate_p_mvn(score_tables[[bs.min]])
                        } else p_mvn <- NA
                        bs = bss[bs.min]
                        score_table = list(score_tables[[bs.min]])
                      }

                      dplyr::tibble(
                        cp_min = cp_min,
                        cp_1se = cp_1se,
                        p_mvn = p_mvn,
                        bs = bs,
                        family = family,
                        score_table = score_table
                      )
                    }, mc.cores = cpo$num_cores)) %>%
   tidyr::unnest(cols = "x")

  cpo[["problematic_targets"]] <-
    cpo$changepoints %>%
    dplyr::filter(is.na(.data$cp_min)) %>%
    dplyr::pull("target_id")

  cpo$changepoints <-
    cpo$changepoints %>%
    dplyr::filter(!is.na(.data$cp_min))

  if(compute_mvn){


    if(!is.null(cpo[["p_adj_method"]])){
      p_adj_method <- cpo[["p_adj_method"]]
    } else {
      p_adj_method <- "BH"
    }

      cpo$p_mvn <-
      cpo$changepoints %>%
      dplyr::select(.data$target_id,p_mvn_target = .data$p_mvn) %>%
      dplyr::left_join(cpo$p_table %>% dplyr::select(tidyr::any_of(c("target_id","gene_id","counts_mean"))),
                       by = "target_id") %>%
      dplyr::relocate(tidyr::any_of(c("target_id","gene_id","counts_mean"))) %>%
      dplyr::mutate(q_mvn_target = stats::p.adjust(.data$p_mvn_target, method = p_adj_method))

    if(cpo$aggregate_to_gene){
      cpo$p_mvn <-
        cpo$p_mvn %>%
        dplyr::group_by(.data$gene_id) %>%
        dplyr::mutate(p_mvn_gene = lancaster(pmax(.data$p_mvn_target,10e-320),.data$counts_mean/sum(.data$counts_mean))) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(q_mvn_gene = stats::p.adjust(.data$p_mvn_gene, method = p_adj_method))
    }
  }

    cpo
}

#' Apply the one-standard-error rule for model selection
#'
#' @param tab A table containing model scores
#' @param nse Number of standard errors to use (default: 1)
#' @return The selected model that is most parsimonious within nse standard errors of the best model
#' @details
#' Implements the one-standard-error rule which selects the most parsimonious model
#' within a specified number of standard errors of the best-scoring model.
#' @keywords internal
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

#' Simulate p-values using multivariate normal distribution
#'
#' @param score_table Table of model scores
#' @param nsim Number of simulations (default: 1e4)
#' @param reg Regularization parameter for covariance matrix (default: 1e-3)
#' @return Probability that the null model is the best scoring model
#' @details
#' Uses Monte Carlo simulation to compute the probability that the null model
#' is the best scoring model under multivariate normal assumptions.
#' @keywords internal
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

#' Calculate model selection scores for candidate changepoints
#'
#' @param data Input data frame
#' @param cps Vector of candidate changepoints
#' @param sp Smoothing parameter
#' @param bs Basis function type
#' @param family Model family ("nb" or "gaussian")
#' @param score Score type to compute
#' @param model_type Type of model to fit
#' @param regularize Regularization parameter
#' @return Table of scores for each candidate changepoint
#' @details
#' Fits models for each candidate changepoint and computes specified
#' model selection scores (AIC or GCV).
#' @keywords internal
#' @noRd
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
    {
      if (length(.) == 0)
        NA
      else
        purrr::map(., ~ do.call(score, args = list(fit = .x))) %>%
        dplyr::bind_cols() %>%
        dplyr::select(dplyr::where(~ !any(is.na(.x)))) %>%
        {
          if (nrow(.) == 0)
            NA
          else
            .
        }
    }

}

#' Calculate AIC for fitted models
#'
#' @param fit A fitted model object (GAM or SCAM)
#' @return AIC value
#' @details
#' Computes AIC for both Gaussian and negative binomial models,
#' handling different parameter estimation approaches for each family.
#' @keywords internal
#' @noRd
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

#' Calculate AIC specifically for negative binomial models
#'
#' @param fit A fitted negative binomial model
#' @return AIC value
#' @details
#' Implements AIC calculation for negative binomial distribution,
#' accounting for dispersion parameter estimation.
#' @keywords internal
#' @noRd
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

#' Calculate pointwise GCV for GAM models
#'
#' @param fit A fitted GAM object
#' @return Vector of pointwise GCV scores
#' @details
#' Implements the standard GCV calculation for GAM models on a pointwise basis.
# computes pointwise generalised cross validation (gcv) for a mgcv::gam fit
#' @keywords internal
#' @noRd
gcv.gam <- function(fit){
    gcvi = (length(fit$residuals)*((sqrt(fit$weights)*fit$residuals)^2))/(length(fit$residuals)-sum(fit$hat))^2
}

#' Calculate GCV score for fitted models
#'
#' @param fit A fitted model object (GAM or SCAM)
#' @return GCV score
#' @details
#' Computes GCV score handling both GAM and SCAM models with appropriate
#' effective degrees of freedom calculations.
#' @keywords internal
#' @noRd
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

#' Calculate log probability density for negative binomial distribution
#'
#' @param x Response values
#' @param theta Dispersion parameter
#' @param mu Mean parameter
#' @return Log probability density values
#' @details
#' Computes the log probability density for a negative binomial distribution
#' adapted for continuous positive response values.
#' @keywords internal
#' @noRd
dnbl <- function(x,theta,mu){
  lgamma(x + theta) -
    lgamma(theta) -
    lgamma(x + 1) +
    theta*log(theta/(theta + mu)) +
    x*log(mu/(mu+theta))
}





