#################################### AIPW: Max U in deterministic treatment regime ########################
aipw_U <- function(eta, x, y, a, prob, tao, delta, Gc_tao, p.es_0, p.es_1) {
  g <- as.numeric(x %*% eta > 0)
  wts <- a * g * 1 / prob + (1 - a) * (1 - g) * (1 / (1 - prob)) # 逆概率权重
  y <- as.numeric(y > tao) # tao标记
  val_tao <- wts * y * (1 - delta)
  val_tao <- val_tao / Gc_tao
  
  p.es <- g * p.es_1 + (1 - g) * p.es_0
  val_mi <- (wts - 1) * p.es # 标准aipw
  
  val <- mean(val_tao - val_mi)
  return(val)
}
GA_U_aipw <- function(x, y, a, prob, tao, delta, Gc_tao, p.es_0, p.es_1) {
  d <- dim(x)[2]
  start_point <- rep(0, d)
  Domains <- cbind(rep(-1,d),rep(1,d))
  
  est<-genoud(fn=aipw_U, nvars=d, x=x,y=y,a=a,prob=prob,tao=tao,delta=delta,Gc_tao=Gc_tao,p.es_0=p.es_0,p.es_1=p.es_1,
              print.level=0,max=TRUE,pop.size=3000,wait.generations=10,gradient.check=FALSE,
              BFGS=FALSE,P1=50, P2=50, P3=10, P4=50, P5=50, P6=50, P7=50, P8=50, P9=0,Domains=Domains,starting.values=start_point,solution.tolerance=0.0001,optim.method="Nelder-Mead")
  eta <- est$par
  return(eta)
}

AIPW_U_free <- function(data, ps_formula, Q_T_formula, Q_cure_formula, Q_T_length, Tx, tao=NULL){
  # ps_formula: a string that represents the formula for the propensity score model?
  # Tx: the treatment variable name
  T <- data$T
  A <- data$A
  Status <- data$Status
  Tx <- as.matrix(cbind(1, data[, Tx]))
  
  fit_ps <- glm(as.formula(ps_formula), family = binomial(link = "logit"), data = data)
  prob <- predict(fit_ps, type = "response")
  
  if (is.null(tao)){
    tao <- max(y[delta==1])
  }
  fit_Gc <- survfit(Surv(T, 1-Status)~1, data=data)
  surv_est <- stepfun(fit_Gc$time, c(1,fit_Gc$surv))
  Gc_tao <- surv_est(tao)
  
  # augment term estimation using gfcure(AFT)
  fit_aft <- gfcure(as.formula(Q_T_formula), as.formula(Q_cure_formula), dist = "lognormal", data = data)
  gamma_p_hat <- -fit_aft$coef[(3+Q_T_length):length(fit_aft$coef)]
  coef_names <- names(gamma_p_hat)
  A_terms <- grepl("\\bA\\b", coef_names) 
  no_A_terms <- !A_terms
  
  logit_p <- Tx%*%gamma_p_hat[no_A_terms]; p.es_0 <- exp(logit_p)/(1+exp(logit_p)) 
  logit_p <- Tx%*%gamma_p_hat[no_A_terms] + Tx%*%gamma_p_hat[A_terms]; p.es_1 <- exp(logit_p)/(1+exp(logit_p))
  
  # genetic algorithm
  d <- dim(Tx)[2]
  start_point <- rep(0, d)
  Domains <- cbind(rep(-1,d),rep(1,d))
  est<-genoud(fn=aipw_U, nvars=d, x=Tx,y=T,a=A,prob=prob,tao=tao,delta=Status,Gc_tao=Gc_tao,p.es_0=p.es_0,p.es_1=p.es_1,
              print.level=0,max=TRUE,pop.size=3000,wait.generations=10,gradient.check=FALSE,
              BFGS=FALSE,P1=50, P2=50, P3=10, P4=50, P5=50, P6=50, P7=50, P8=50, P9=0,Domains=Domains,starting.values=start_point,solution.tolerance=0.0001,optim.method="Nelder-Mead")
  eta <- est$par
  return(list(eta, est$value))
}