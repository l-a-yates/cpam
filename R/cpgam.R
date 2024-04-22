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
                  silent = T){

  family <- match.arg(family)
  model_type <- match.arg(model_type)
  bs_scam <- c("micv","mdcx","cv","cx","micx","mdcv")
  bs <- match.arg(bs,c("tp","null",bs_scam))
  use_scam <-  !bs %in% c("tp","null")

  if (family == "nb") {
    if (!regularize) {
      gam_family <- mgcv::nb()
      if (bs != "tp") {
        warning(
          paste0(
            "Shape-constrained splines are unavilable for negative-binomial models when the dipsersion must be estimated. The basis has been to set
              to bs = 'tp' instead of the supplied value bs = ",bs))
        bs <- "tp"
        use_scam <- F
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
    gam_optimizer <- if(use_scam) "bfgs" else "outer"
  } else if (k == 2){
    f <- paste0(resp," ~ 1 + td")
    sp <- NULL
    gam_optimizer <- if(use_scam) "bfgs" else "outer"
  } else {
    if(strsplit(bs,"")[[1]][1] == "m") k <- 2*k
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

  if(!use_scam){
    m = try(mgcv::gam(formula = f,
                      data = data,
                      family = gam_family,
                      sp = sp,
                      method = gam_method,
                      optimizer = gam_optimizer,
                      offset = log(data$norm_factor)),
            silent = silent)
  } else if(bs %in% c("micv","mdcx","cv","cx","micx","mdcv")){
    m = try(scam::scam(formula = f,
                       data = data,
                       family = gam_family,
                       sp = sp,
                       optimizer = gam_optimizer,
                       offset = log(data$norm_factor)),
            silent = silent)
  }

  if(inherits(m, "try-error")){
    warning(paste0("The changepoint additive model for gene_id ",data$gene_id[1], " with basis '",bs,"' failed to converged at t0 =",cp))
    return(NA)
  }
  m$f <- f
  m$cp <- cp
  m$bs <- bs
  m$data <- data
  m$model_type <- model_type
  m
}

