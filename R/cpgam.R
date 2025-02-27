#' Fit a changepoint GAM model (internal)
#'
#' @description
#' Internal function to fit a generalized additive model (GAM) with changepoint at a specified time.
#' Supports both case-only and case-control designs with negative binomial or gaussian families.
#'
#' @param data A data frame containing the time series data
#' @param cp Numeric changepoint time
#' @param regularize Logical; whether to use regularization for negative binomial models
#' @param model_type Character; either "case-only" or "case-control"
#' @param bs Character; basis type. Can be GAM bases ("tp", "null", "lin", "ilin", "dlin")
#'           or SCAM bases ("micv", "mdcx", "cv", "cx", "micx", "mdcv")
#' @param family Character; either "nb" (negative binomial) or "gaussian"
#' @param gam_method Character; fitting method for GAM, defaults to "REML"
#' @param gam_optimizer Character; optimization method, defaults to "efs"
#' @param fixed_effects Character; additional fixed effects to include in model formula
#' @param sp Numeric; smoothing parameter(s)
#' @param k_mult Numeric; multiplier for basis dimension k
#' @param n_try Integer; number of fitting attempts
#' @param not_exp Logical; passed to scam::scam
#' @param silent Logical; whether to suppress warnings
#' @param debug Logical; whether to print debug messages
#'
#' @return A fitted GAM or SCAM model object with additional components:
#' \itemize{
#'   \item f: model formula
#'   \item cp: changepoint time
#'   \item bs: basis type used
#'   \item data: input data
#'   \item model_type: type of model fitted
#'   \item se_ratio: ratio of maximum standard error to maximum counts
#' }
#' Returns NA if model fitting fails.
#'
#' @keywords internal
#' @noRd
#'
cpgam <- function(data,
                  cp,
                  regularize,
                  model_type = c("case-only","case-control"),
                  bs = "tp",
                  family = c("nb","gaussian"),
                  gam_method = "REML",
                  gam_optimizer = "efs",
                  fixed_effects = NULL,
                  sp = NULL,
                  k_mult = 1.2,
                  n_try = 1,
                  not_exp = F,
                  silent = T,
                  debug = F){

  family <- match.arg(family)
  model_type <- match.arg(model_type)
  bs_scam <- c("micv","mdcx","cv","cx","micx","mdcv")
  bs_gam <- c("tp","null","lin","ilin","dlin")
  bs <- match.arg(bs,c(bs_gam,bs_scam))
  use_scam <-  !bs %in% bs_gam
  if(is.null(data[["time"]])) stop("data must contain a 'time' column")

  if (family == "nb") {
    if (!regularize) {
      gam_family <- mgcv::nb()
      if (bs != "tp") {
        warning(
          paste0(
            "Shape-constrained splines are unavilable for negative-binomial models when the dipsersion must be estimated. The basis has been to set
              to bs = 'tp' instead of the supplied value bs = ",bs))
        bs <- "tp"
        use_scam <- T
      }
    } else {
      gam_family <- mgcv::negbin(theta = as.numeric(1/data$disp[1]))
    }
    resp <- "counts"
  } else{
    gam_family <- stats::gaussian()
    resp <- "y"
    if(!is.null(sp)){
      gam_optimizer <- if(use_scam) "bfgs" else "outer"
    }
    gam_method <- "GCV.Cp"
  }

  k = sum(unique(data$time) >= cp)

  if (k == 1 | bs == "null"){
    f <- paste0(resp," ~ 1")
    sp <- NULL
    bs <- "null"
    use_scam <- F
    cp <- 0
    #gam_optimizer <- "outer"
  } else if (k == 2 | bs %in% c("lin","ilin","dlin")){
    f <- paste0(resp," ~ 1 + td")
    if(!bs %in% c("lin","ilin","dlin")) bs <- "lin"
    sp <- NULL
    use_scam <- F
    #gam_optimizer <- if(use_scam) "bfgs" else "outer"
  } else {
    if(bs != "tp") k <- k_mult*k
    f <- paste0(resp," ~ s(td, bs = '",bs,"', k =", k,")")
  }

  if(model_type == "case-control"){
    data <- data %>% dplyr::mutate(td = pmax(.data$time*.data$case,cp),
                                   y = log(.data$counts + 0.5))
    f <- paste(f,"+ s(time, k = ", length(unique(data$time)),")")
  } else{
    data <- data %>% dplyr::mutate(td = pmax(.data$time,cp),
                                   y = log(.data$counts + 0.5))
  }

  if(!is.null(fixed_effects)) f <- paste(f,"+",fixed_effects)

  f <- stats::as.formula(f)

  refit = F

  if(!use_scam){
    m = try(mgcv::gam(formula = f,
                      data = data,
                      family = gam_family,
                      sp = sp,
                      method = gam_method,
                      optimizer = gam_optimizer,
                      offset = log(data$norm_factor)) %>% suppressWarnings(),
            silent = silent)

    if(inherits(m, "try-error") & gam_optimizer != "outer" & n_try < 3){
      refit = T
      n_try = 3
      gam_optimizer = "outer"
      if(debug) print("refit-C")
    }
  } else {
    m = try(scam::scam(formula = f,
                       data = data,
                       family = gam_family,
                       sp = sp,
                       optimizer = gam_optimizer,
                       not.exp = not_exp,
                       offset = log(data$norm_factor)) %>% suppressWarnings(),
            silent = silent)

    if(inherits(m, "try-error")){
      if(stringr::str_detect(m,"Error in Rrank") & k_mult != 1 & n_try == 1){
        k_mult = 1
        n_try = 2
        refit = T
        if(debug) print("refit-A")
      } else if(n_try <= 2){
        n_try = 3
        not_exp = T
        refit = T
        if(debug) print("refit-B")
      } else refit = F
      }
  }

  if(refit){
    m <- cpgam(
      data = data,
      cp = cp,
      regularize = regularize,
      model_type = model_type,
      bs = bs,
      family = family,
      gam_method = gam_method,
      gam_optimizer = gam_optimizer,
      fixed_effects = fixed_effects,
      sp = sp,
      k_mult = k_mult,
      n_try = n_try,
      not_exp = not_exp,
      silent = silent
    )
  }

  if(inherits(m, "try-error")){
    #warning(paste0("The changepoint additive model for gene_id ",data$gene_id[1], " with basis '",bs,"' failed to converged at t0 =",cp))
    return(NA)
  }

  if(inherits(m,c("gam","scam"))){
    m$f <- f
    m$cp <- cp
    m$bs <- bs
    m$data <- data
    m$model_type <- model_type
    if(inherits(m,"scam")){
      m$se_ratio <- max(scam::predict.scam(m,se = T)$se.fit)/max(data$counts)
    } else {
      m$se_ratio <- max(mgcv::predict.gam(m,se = T)$se.fit)/max(data$counts)
    }

  }

  m
}
