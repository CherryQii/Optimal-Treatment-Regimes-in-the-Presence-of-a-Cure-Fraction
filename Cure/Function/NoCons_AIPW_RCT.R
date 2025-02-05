aipw_U_rct <- function(x, y, a, prob, tao, delta, Gc_tao, p.es_0, p.es_1) {
  g <- rbinom(n = length(a), size = 1, prob = 0.5)
  wts <- a * g * 1 / prob + (1 - a) * (1 - g) * (1 / (1 - prob))
  y <- as.numeric(y > tao)
  val_tao <- wts * y * (1 - delta)
  val_tao <- val_tao / Gc_tao
  
  p.es <- g * p.es_1 + (1 - g) * p.es_0
  val_mi <- (wts - 1) * p.es
  
  val <- mean(val_tao - val_mi)
  return(val)
}

aipw_T_rct <- function(x, y, a, prob, tao, delta, Gc_y, Gc_tao, t.es_0, t.es_1, p.es_0, p.es_1) {
  g <- rbinom(n = length(a), size = 1, prob = 0.5)
  wts <- a * g * 1 / prob + (1 - a) * (1 - g) * (1 / (1 - prob))
  t.es <- g * t.es_1*p.es_0 + (1 - g) * t.es_0*p.es_1
  # .T <- ifelse(wts * y * delta / Gc_y > tao, tao, wts * y * delta / Gc_y)
  .T <- wts * y * delta / Gc_y
  v_t <- mean(.T + (1 - wts) * t.es)
  
  v_u <- aipw_U_rct(x, y, a, prob, tao, delta, Gc_tao, p.es_0, p.es_1)
  val_t <- v_t / (1 - v_u)
  
  return(val_t)
}

AIPW_RCT <- function(data, tao = NULL, ps_formula, Q_T_formula, Q_cure_formula, Q_T_length, Tx) {
  # 1. data load
  T <- data$T
  A <- data$A
  Status <- data$Status
  Tx <- as.matrix(cbind(1, data[, Tx]))
  
  # 2. tao and Gc estimation
  if (is.null(tao)){
    tao <- max(y[Status==1])
  }
  fit_Gc <- survfit(Surv(T, 1-Status)~1, data=data)
  surv_est <- stepfun(fit_Gc$time, c(1,fit_Gc$surv))
  Gc_tao <- surv_est(tao)
  Gc_y <- surv_est(T)
  Gc_y[Gc_y==0]=min(Gc_y[Gc_y!=0])
  
  # 3. propensity score estimation
  fit_ps <- glm(as.formula(ps_formula), family = binomial, data = data)
  prob <- predict(fit_ps, type = "response")
  
  # 4. augment term estimation using gfcure(AFT)
  fit_aft <- gfcure(as.formula(Q_T_formula), as.formula(Q_cure_formula), dist = "lognormal", data = data)
  gamma_p_hat <- -fit_aft$coef[(3+Q_T_length):length(fit_aft$coef)]
  coef_names <- names(gamma_p_hat)
  A_terms <- grepl("\\bA\\b", coef_names) ; no_A_terms <- !A_terms
  logit_p <- Tx%*%gamma_p_hat[no_A_terms]; p.es_0 <- exp(logit_p)/(1+exp(logit_p)) 
  logit_p <- Tx%*%gamma_p_hat[no_A_terms] + Tx%*%gamma_p_hat[A_terms]; p.es_1 <- exp(logit_p)/(1+exp(logit_p))
  gamma_t_hat <- fit_aft$coef[2:(Q_T_length+2)]
  coef_names <- names(gamma_t_hat)
  A_terms <- grepl("\\bA\\b", coef_names) ; no_A_terms <- !A_terms
  t.es_0 <- exp(Tx%*%gamma_t_hat[no_A_terms]) * exp(exp(fit_aft$coef[1])^2/2)
  t.es_1 <- exp(Tx%*%gamma_t_hat[no_A_terms] + 1* Tx %*%gamma_t_hat[A_terms]) * exp(exp(fit_aft$coef[1])^2/2)
  
  result_U = aipw_U_rct(Tx, T, A, prob, tao, Status, Gc_tao, p.es_0,p.es_1)
  result_T = aipw_T_rct(Tx, T, A, prob, tao, Status, Gc_y, Gc_tao,t.es_0, t.es_1, p.es_0, p.es_1)
  
  return(list(EU_RCT = result_U, ET_RCT = result_T))
}
