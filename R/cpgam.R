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
                  silent = T){

  family <- match.arg(family)
  model_type <- match.arg(model_type)
  bs_scam <- c("micv","mdcx","cv","cx","micx","mdcv")
  bs_gam <- c("tp","null","lin","ilin","dlin")
  bs <- match.arg(bs,c(bs_gam,bs_scam))
  use_scam <-  !bs %in% bs_gam

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

  if(!use_scam){
    m = try(mgcv::gam(formula = f,
                      data = data,
                      family = gam_family,
                      sp = sp,
                      method = gam_method,
                      optimizer = gam_optimizer,
                      offset = log(data$norm_factor)) %>% suppressWarnings(),
            silent = silent)
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
        #print("refit-A")
      } else if(n_try <= 2){
        n_try = 3
        not_exp = T
        refit = T
        #print("refit-B")
      } else refit = F

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
      }
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
    m$se_ratio <- max(predict(m,se = T)$se.fit)/max(data$counts)
  }

  m
}
