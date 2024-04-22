predict_cpgam <- function(fit,
                          ci_prob = 2 * pnorm(1) - 1,
                          length.out = 500,
                          lim_factor = 5,
                          nsim = 2e3) {


  if(inherits(fit,"scam")){
    beta <- fit$coefficients.t
    V <- fit$Vp.t
  } else{
    beta <- fit$coefficients
    V <- fit$Vp
  }

  newdata <- dplyr::tibble(time = seq(0, max(fit$data$time),
                                      length.out = length.out),
                           td = pmax(fit$cp, .data$time))

  if(fit$model_type == "case-control"){
    newdata <- tidyr::expand_grid(case = 0:1, newdata) %>%
      dplyr::mutate(td = pmax(.data$time*.data$case,fit$cp))
  }

  Xp = stats::predict(fit, newdata = newdata, type = "lpmatrix")
  br = mgcv::rmvn(nsim, mu = as.vector(beta),V = V)
  y_matrix_exp = exp(Xp %*% t(br))

  newdata %>%
    dplyr::mutate(counts = as.numeric(exp(Xp %*% beta)),
                  q_lo =  y_matrix_exp %>% matrixStats::rowQuantiles(probs = 0.5 - ci_prob/2),
                  q_hi =  y_matrix_exp %>% matrixStats::rowQuantiles(probs = 0.5 + ci_prob/2),
                  q_hi = pmin(.data$q_hi, lim_factor*max(fit$y))
    )

}
