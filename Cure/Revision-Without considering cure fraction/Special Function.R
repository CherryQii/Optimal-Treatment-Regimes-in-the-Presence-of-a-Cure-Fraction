#################### IPW U & IPW T############################
ipw_U <- function(eta, x, y, a, prob, tao, delta, G_c_tao) {
  g <- as.numeric(x %*% eta > 0)
  wts <- a* g * 1 / prob + (1-a) * (1 - g) * (1 / (1 - prob))
  wts <- wts / G_c_tao # 删失逆概率
  y <- as.numeric(y > tao) # tao标记
  val <- mean(wts * y * (1 - delta))
  val
}
ipw_T <- function(eta, x, y, a, prob, tao, delta, G_c, G_c_tao) {
  g <- as.numeric(x %*% eta > 0)
  wts <- a* g * 1 / prob + (1-a) * (1 - g) * (1 / (1 - prob))
  v_t <- mean(wts * y * delta / (G_c))
  v_u <- ipw_U(eta, x, y, a, prob, tao, delta, G_c_tao)
  val <- v_t/(1 - v_u)
  val
}
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
aipw_T <- function(eta, x, y, a, prob, tao, delta, G_c, G_c_tao, t.es_0,t.es_1,p.es_0,p.es_1) {
  g <- as.numeric(x %*% eta > 0)
  wts <- a* g * 1 / prob + (1-a) * (1 - g) * (1 / (1 - prob))
  v_t <- wts * y * delta / (G_c)
  H <- g * t.es_1* p.es_0 + (1-g) * t.es_0 * p.es_1
  val_t <- mean(v_t + (1-wts)*H)
  v_u <- aipw_U(eta, x, y, a, prob, tao, delta, G_c_tao, p.es_0,p.es_1)
  val <- val_t/(1-v_u)
  return(val)
}
GA_T_aipw <- function(x, y, a, prob, tao, delta, G_c, G_c_tao,t.es_0, t.es_1, p.es_0, p.es_1) {
  d <- dim(x)[2]
  start_point <- rep(0, d)
  Domains <- cbind(rep(-1,d),rep(1,d))
  
  est<-genoud(fn=aipw_T, nvars=d, x=x,y=y,a=a,prob=prob,tao=tao,delta=delta,G_c=G_c,G_c_tao=G_c_tao,t.es_0=t.es_0,t.es_1=t.es_1,p.es_0=p.es_0,p.es_1=p.es_1,
              print.level=0,max=TRUE,pop.size=3000,wait.generations=10,gradient.check=FALSE,
              BFGS=FALSE,P1=50, P2=50, P3=10, P4=50, P5=50, P6=50, P7=50, P8=50, P9=0,Domains=Domains,starting.values=start_point,solution.tolerance=0.0001,optim.method="Nelder-Mead")
  eta <- est$par
  return(eta)
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
#################### IPW Max ET without considering cure rate ############################
ipw_ET_without_cure <- function(eta, x, y, a, prob, delta, Gc_y){
  g <- as.numeric(x %*% eta > 0)
  wts <- a* g * 1 / prob + (1-a) * (1 - g) * (1 / (1 - prob))
  wts <- wts / Gc_y
  val <- mean(wts * y * delta)
  val
}
GA_ipw_ET <- function(x, y, a, prob, delta, Gc_y) {
  d <- dim(x)[2]
  start_point <- rep(0, d)
  Domains <- cbind(rep(-1,d),rep(1,d))
  
  est<-genoud(fn=ipw_ET_without_cure, nvars=d, x=x,y=y,a=a,prob=prob,delta=delta,Gc_y=Gc_y,
              print.level=0,max=TRUE,pop.size=3000,wait.generations=10,gradient.check=FALSE,
              BFGS=FALSE,P1=50, P2=50, P3=10, P4=50, P5=50, P6=50, P7=50, P8=50, P9=0,Domains=Domains,starting.values=start_point,solution.tolerance=0.0001,optim.method="Nelder-Mead")
  eta <- est$par
  value <- est$value
  return(list(eta, value))
}

#################### Main ############################
df <- read.csv("seer_esophageal_clean.csv")
attach(df)
# x <- cbind(1,age, stage_r,stage_d)
x <- cbind(1, Age, stage_r,stage_d)
y <- Survival.months
a <- A1
fit_pi <- glm(a~age+Sex+Stage+Histology, family = binomial(link = "logit"),data = df)
prob <- predict(fit_pi, type = "response")
delta <- Status
tao <- 120
fit_Gc <- survfit(Surv(Survival.months, 1-Status)~1, data=df)
surv_est <- stepfun(fit_Gc$time, c(1,fit_Gc$surv))
Gc_tao <- surv_est(tao)
Gc_y <- surv_est(y)
Gc_y[Gc_y==0]=min(Gc_y[Gc_y!=0])
detach(df)
fit_aft2 <- gfcure(Surv(Survival.months, Status)~ Age*A1 + stage_r*A1 + stage_d*A1,~Age*A1 + stage_r*A1 + stage_d*A1,data= df, dist = "lognormal")
gamma_p_hat <- -fit_aft2$coef[10:17]
logit_p <- x%*%gamma_p_hat[c(1,2,4,5)]; p.es_0 <- exp(logit_p)/(1+exp(logit_p))
logit_p <- x%*%gamma_p_hat[c(1,2,4,5)] + x%*%gamma_p_hat[c(3,6,7,8)]; p.es_1 <- exp(logit_p)/(1+exp(logit_p))
gamma_t_hat <- fit_aft2$coef[2:9]
t.es_0 <- exp(x%*%gamma_t_hat[c(1,2,4,5)]) * exp(exp(fit_aft2$coef[1])^2/2)
t.es_1 <- exp(x%*%gamma_t_hat[c(1,2,4,5)] + 1* x %*%gamma_t_hat[c(3,6,7,8)]) * exp(exp(fit_aft2$coef[1])^2/2)


# Max ET without considering cure rate
set.seed(2023)
result_ET_without_cure = GA_ipw_ET(x, y, a, prob, delta, Gc_y)
ipw_ET_without_cure(result_ET_without_cure[[1]], x, y, a, prob, delta, Gc_y)
aipw_U(result_ET_without_cure[[1]],x, y, a, prob, tao, delta, Gc_tao, p.es_0,p.es_1)
aipw_T(result_ET_without_cure[[1]], x, y, a, prob, tao, delta, Gc_y, Gc_tao,t.es_0, t.es_1, p.es_0, p.es_1)
# Max ET & EU considering cure rate
set.seed(2023)
result_EU = GA_U_aipw(x, y, a, prob, tao, delta, Gc_tao, p.es_0, p.es_1)
ipw_ET_without_cure(result_EU, x, y, a, prob, delta, Gc_y)
aipw_U(result_EU,x, y, a, prob, tao, delta, Gc_tao, p.es_0,p.es_1)
aipw_T(result_EU, x, y, a, prob, tao, delta, Gc_y, Gc_tao,t.es_0, t.es_1, p.es_0, p.es_1)
set.seed(2023)
result_ET = GA_T_aipw(x, y, a, prob, tao, delta, Gc_y, Gc_tao, t.es_0, t.es_1, p.es_0, p.es_1)
ipw_ET_without_cure(result_ET, x, y, a, prob, delta, Gc_y)
aipw_U(result_ET,x, y, a, prob, tao, delta, Gc_tao, p.es_0,p.es_1)
aipw_T(result_ET, x, y, a, prob, tao, delta, Gc_y, Gc_tao,t.es_0, t.es_1, p.es_0, p.es_1)

# Consider cure rate and constraint AIPW
