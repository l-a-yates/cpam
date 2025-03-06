predict_cpgam <- function(fit,
                          ci_prob = "se",
                          length.out = 500,
                          lim_factor = 2,
                          nsim = 2e3,
                          scaled = FALSE,
                          logged = FALSE){


  newdata <- dplyr::tibble(time = seq(min(fit$data$time), max(fit$data$time),
                                      length.out = length.out),
                           td = pmax(fit$cp, .data$time))

  if(fit$model_type == "case-control"){
    newdata <- tidyr::expand_grid(case = 0:1, newdata) %>%
      dplyr::mutate(td = pmax(.data$time*.data$case,fit$cp))
  }

  od <- 1
  if("overdispersions" %in% names(fit$data) & !scaled){
    od <- fit$data$overdispersions[1] %>% as.numeric
  }

  if(ci_prob == "se"){

    if(logged){

      pred = stats::predict(fit, newdata = newdata, se.fit = TRUE, type = "link")

      newdata %>%
        dplyr::mutate(counts = pred$fit + log(od),
                      q_lo =  pred$fit + log(od) - pred$se.fit,
                      q_hi =  pred$fit + log(od) + pred$se.fit
        ) %>% return()
    } else{
    pred = stats::predict(fit, newdata = newdata, se.fit = TRUE, type = "response")

    newdata %>%
      dplyr::mutate(counts = pred$fit*od,
                    q_lo =  pred$fit*od - pred$se.fit*sqrt(od),
                    q_hi =  pred$fit*od + pred$se.fit*sqrt(od),
                    q_lo = pmax(.data$q_lo, 0)
                    ##q_hi = pmin(.data$q_hi, lim_factor*max(fit$data$counts)*od
      ) %>% return()
    }

  } else {

    if(!is.numeric(ci_prob)) stop("If not equal to 'se', ci_prob must be numeric.")

    if(inherits(fit,"scam")){
      beta <- fit$coefficients.t
      V <- fit$Vp.t
    } else{
      beta <- fit$coefficients
      V <- fit$Vp
    }

    Xp = stats::predict(fit, newdata = newdata, type = "lpmatrix")
    br = mgcv::rmvn(nsim, mu = as.vector(beta),V = V)
    if(logged){
      y_matrix <- (Xp %*% t(br)) + log(od)
      counts <- as.numeric(Xp %*% beta) + log(od)
    } else {
      y_matrix <- exp(Xp %*% t(br))*od
      counts <- as.numeric(exp(Xp %*% beta))*od
    }

    newdata %>%
      dplyr::mutate(counts = counts,
                    q_lo =  y_matrix %>% matrixStats::rowQuantiles(probs = 0.5 - ci_prob/2),
                    q_hi =  y_matrix %>% matrixStats::rowQuantiles(probs = 0.5 + ci_prob/2)) %>%
      #dplyr::mutate(q_hi = pmin(.data$q_hi, lim_factor*max(fit$data$counts)*od)) %>%
      return()

  }

}
