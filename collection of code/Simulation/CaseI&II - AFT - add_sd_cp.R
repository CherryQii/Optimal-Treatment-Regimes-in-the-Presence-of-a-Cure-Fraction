########################## package loading ################################
library(mixcure)
library(rgenoud)
library(parallel);library(doParallel);library(foreach);library(doRNG)
cores_used <- 10;cl <- makeCluster(cores_used);registerDoParallel(cl)
clusterEvalQ(cl, library(rgenoud));clusterEvalQ(cl, library(mixcure));clusterEvalQ(cl, library(gfcure))
# stopCluster(cl)

########################## parameters setting ################################
gamma_p = c(c(-1,0.25,0.25),c(0,1.5,-1.5));gamma_t=c(c(0,0.25,0.25),c(0,1,1))
lambda_c = 0.03
########################## sample ################################
set.seed(2023)
EDA(gamma_p, gamma_t, lambda_c, dist = 'Weibull')
########################## IPW ################################
ipw_U <- function(eta, x, y, a, prob, tao, delta, G_c_tao) {
  g <- as.numeric(x %*% eta > 0)
  wts <- a* g * 1 / prob + (1-a) * (1 - g) * (1 / (1 - prob)) # 逆概率权重
  wts <- wts / G_c_tao # 删失逆概率
  y <- as.numeric(y > tao) # tao标记
  val <- mean(wts * y * (1 - delta))
  val
}
ipw_T <- function(eta, x, y, a, prob, tao, delta, G_c, G_c_tao) {
  g <- as.numeric(x %*% eta > 0)
  wts <- a* g * 1 / prob + (1-a) * (1 - g) * (1 / (1 - prob)) # 逆概率权重
  v_t <- mean(wts * y * delta / (G_c))
  v_u <- ipw_U(eta, x, y, a, prob, tao, delta, G_c_tao)
  val <- v_t/(1 - v_u)
  val
}


# 最后返回一个向量\psi2_i
estimate_psi2_i <- function(i, x, y, delta, prob, Gc_y, tao, Lambda_c) {
  # y: 观察时间
  # delta: 事件指示变量
  # prob: 估计的概率（可能是IPW估计的倾向性得分）
  # S_c: 估计的生存函数
  # tao: 阈值
  # Lambda_c: 累计风险函数估计
  
  # 首先，将Y、J_i_c、Lambda_c排序，以便计算积分部分
  sorted_indices <- order(y)  # 按y排序
  y_sorted <- y[sorted_indices]
  delta_sorted <- delta[sorted_indices]
  prob_sorted <- prob[sorted_indices]
  Lambda_c_sorted <- Lambda_c[sorted_indices]  # 累计风险函数
  
  # 初始化结果变量
  psi2_values <- numeric(length(y))
  
  # 计算积分部分，M_c_i(vector)
  for (j in 1:length(y_sorted)) {
    nu_Y = min(y_sorted[j], tao)
    if (y_sorted[j] <= nu_Y){
      N_c_i_j = ifelse(y_sorted[j] <= nu_Y && delta_sorted[j] == 0, 1, 0)
      integral_part = Lambda_c_stepfun(y_sorted[j])
    } else {
      N_c_i_j = 0
      integral_part = Lambda_c_stepfun(nu_Y)
    }
    M_c_i_j = N_c_i_j - integral_part
    psi2_values[j] = -Gc_y_stepfun(nu_Y) * M_c_i_j / mean(y_sorted[j] > nu_Y)
  }
  
  # 恢复原来的顺序
  psi2_values_final <- psi2_values[order(sorted_indices)]
  
  return(psi2_values_final)
}

estimate_psi2_vector <- function(x, y, delta, prob, S_c, tao) {
  psi_2i_vector <- sapply(1:length(y), function(i) {
    estimate_psi2(x[i,], y[i], delta[i], prob[i], S_c, tao)
  })
  return(psi_2i_vector)
}

sd_ipw_U <- function(eta, x, y, a, prob, tao, delta, G_c_tao, theta, S_c) {
  g <- as.numeric(x %*% eta > 0)
  wts <- a* g * 1 / prob + (1-a) * (1 - g) * (1 / (1 - prob))
  wts2 <- wts^2
  y <- as.numeric(y > tao)
  n <- length(y)
  
  term1 <- wts * y * (1 - delta)/G_c_tao
  term2 <- ipw_U(eta, x, y, a, prob, tao, delta, G_c_tao)
  D1 <- colMeans(-wts2 * y * (1 - delta)/G_c_tao * (2*g-1) * prob * (1 - prob) * x)
  D2 <- -wts * y * (1 - delta)/(G_c_tao)^2 / n
  H1 <- t(x) %*% (prob * (1 - prob) * x) / n
  phi1 <- solve(H1) %*% t(x) * (A - prob) 
  term3 <- sum(D1 %*% phi1)
  term4 <- 0
  for (i in 1:n) {
    phi2_i <- estimate_psi2_i(i, x, y, delta, prob, S_c, tao, theta)
   term4 <- term4 + D2 %*% phi2_i
  }
  sd = 
}
########################## IPW-U estimator ################################
GA_U_ipw <- function(x, y, a, prob, tao, delta, G_c_tao) {
  d <- dim(x)[2]
  start_point <- rep(1, d)
  Domains <- cbind(rep(-1,d),rep(1,d))
  
  est<-genoud(fn=ipw_U, nvars=d, x=x,y=y,a=a,prob=prob,tao=tao,delta=delta,G_c_tao=G_c_tao,
              print.level=0,max=TRUE,pop.size=3000,wait.generations=10,gradient.check=FALSE,
              BFGS=FALSE,P1=50, P2=50, P3=10, P4=50, P5=50, P6=50, P7=50, P8=50, P9=0,Domains=Domains,starting.values=start_point,solution.tolerance=0.0001,optim.method="Nelder-Mead")
  eta <- est$par
  return(eta)
}
sim_U_ipw <- function(n = 500, gamma_p, gamma_t, lambda_c,sim_num = 100, mis_ps = FALSE, dist = 'norm'){
  estv <- 0;mc_estv <- 0;truev <- 0;ratio <- c()
  eta <- c();value <- c();AA <- c()
  for (i in 1:sim_num) {
    dat <- sim_dat(n, gamma_p=gamma_p, gamma_t=gamma_t, lambda_c=lambda_c, dist = dist)
    dat_obs <- dat[[1]]
    dat_true <- dat[[2]]
    
    attach(dat_obs)
    x <- cbind(1, x1, x2)
    
    tao <- max(T_obs[delta == 1])
    fit_Gc <- survfit(Surv(T_obs, 1-delta)~1, data=dat_obs)
    surv <- stepfun(fit_Gc$time, c(1, fit_Gc$surv), right=FALSE, f=0)
    G_c <- surv(T_obs); G_c_tao <- surv(tao)
    G_c[G_c==0] <- min(G_c[G_c!=0])
    # propensity score正确给定
    fit <- glm(A~I(x1^2) + I(x2^2), family = binomial(link = "logit"))
    # if (fit$converged == FALSE){break}
    prob <- predict(fit, type = "response")
    
    if (mis_ps == FALSE){
      eta_hat <- GA_U_ipw(x, T_obs, A, prob, tao, delta, G_c_tao)
    } else {
      # propensity score错误给定
      fit <- glm(A~1, family = binomial(link = "logit"))
      prob <- predict(fit, type = "response")
      # prob <- rep(0.5,n)
      eta_hat <- GA_U_ipw(x, T_obs, A, prob, tao, delta, G_c_tao)
    }
    eta_hat <- eta_hat/sum(abs(eta_hat)) # 一范数归一化
    
    estv = ipw_U(eta_hat, x, T_obs, A, prob,tao,delta,G_c_tao)
    mc_estv <- eval_performance_g(eta_hat, gamma_p. = gamma_p,gamma_t. = gamma_t, dist = dist)[1]
    V_hat_g_opt <- ipw_U(gamma_p[4:6], x, T_obs, A, prob ,tao,delta,G_c_tao)
    truev <- eval_performance_g(gamma_p[4:6], gamma_p. = gamma_p,gamma_t. = gamma_t, dist = dist)[1]
    AA_given <- eval_performance_g(eta_hat, gamma_p. = gamma_p,gamma_t. = gamma_t, dist = dist)[3]
    
    eta <- rbind(eta, eta_hat)
    value <- rbind(value, c(estv, mc_estv, V_hat_g_opt, truev))
    AA <- c(AA, AA_given)
    colnames(value) <- c("V_hat(g_hat_opt)", "V(g_hat_opt)", "V_hat(g_opt)","V(g_opt)")
    
    detach(dat_obs)
    print(i)
  }
  list(eta, value,AA)
}
sim_U_ipw(n=1000, gamma_p, gamma_t, lambda_c,sim_num = 1, mis_ps = FALSE, dist = "norm") # Test
sim_U_ipw(n=1000, gamma_p, gamma_t, lambda_c,sim_num = 1, mis_ps = FALSE, dist = "Weibull") # Test
# 1.ipw-U-true
registerDoRNG(2023);a0 = proc.time()[3]
result_ipwU_t = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  sim_U_ipw(n=1000, gamma_p, gamma_t, lambda_c,sim_num = 50, mis_ps = FALSE, dist="Weibull")};print(proc.time()[3] - a0)

`#  1.ipw-U-false
registerDoRNG(2023);a0 = proc.time()[3]
result_ipwU_f = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  sim_U_ipw(n=1000, gamma_p, gamma_t, lambda_c,sim_num = 50, mis_ps = TRUE, dist="Weibull")};print(proc.time()[3] - a0)

result = break_result(result_ipwU_f)
cat(
  "eta", apply(result[[1]], 2, mean),"\n",
  "eta_sd", apply(result[[1]], 2, sd),"\n",
  "V", apply(result[[2]], 2, mean),"\n",
  "V_sd", apply(result[[2]], 2, sd),"\n",
  "AA:", mean(result[[3]]), "\n",
  "AA_sd", sd(result[[3]]))
saveRDS(result, "result\\Scenario_1\\Weibull\\ipw-U-t.RDS")
saveRDS(result, "result\\Scenario_1\\Weibull\\ipw-U-f.RDS")
########################## IPW-T estimator ################################
GA_T_ipw <- function(x, y, a, prob, tao, delta, G_c, G_c_tao) {
  d <- dim(x)[2]
  start_point <- rep(0, d)
  # start_point <- c(0, -0.5, 0.5)
  Domains <- cbind(rep(-1,d),rep(1,d))
  
  est<-genoud(fn=ipw_T, nvars=d, x=x,y=y,a=a,prob=prob,tao=tao,delta=delta,G_c=G_c, G_c_tao=G_c_tao,
              print.level=0,max=TRUE,pop.size=3000,wait.generations=10,gradient.check=FALSE,
              BFGS=FALSE,P1=50, P2=50, P3=10, P4=50, P5=50, P6=50, P7=50, P8=50, P9=0,Domains=Domains,starting.values=start_point,solution.tolerance=0.0001,optim.method="Nelder-Mead")
  eta <- est$par
  return(eta)
}
sim_T_ipw <- function(n = 500, gamma_p=gamma_p, gamma_t=gamma_t, lambda_c, sim_num = 100, ps_model=TRUE, dist = 'norm'){
  T_hat<-c();T_true<-c();T_hat_g_opt<-c();T_true_g_opt<-c();
  eta <- c();value<-c();
  AA <- c()
  for (i in 1:sim_num) {
    dat <- sim_dat(n, gamma_p=gamma_p, gamma_t=gamma_t, lambda_c=lambda_c, dist = dist)
    dat_obs <- dat[[1]]
    attach(dat_obs)
    x <- cbind(1, x1, x2)
    tao <- max(T_obs[delta == 1])
    fit_Gc <- survfit(Surv(T_obs, 1-delta)~1, data=dat_obs)
    surv <- stepfun(fit_Gc$time, c(1, fit_Gc$surv), right=FALSE, f=0)
    G_c <- surv(T_obs); G_c_tao <- surv(tao);
    G_c[G_c==0] <- min(G_c[G_c!=0])
    # propensity score正确给定
    fit <- glm(A~I(x1^2) + I(x2^2), family = binomial(link = "logit"))
    # if (fit$converged == FALSE){break}
    prob <- predict(fit, type = "response")
    
    if (ps_model == TRUE){
      eta_hat <- GA_T_ipw(x, T_obs, A, prob, tao, delta, G_c, G_c_tao)
    } else {
      # prob = mean(dat_obs$A)
      # prob <- rep(0.25, n)
      fit <- glm(A~1, family = binomial(link = "logit"))
      prob <- predict(fit, type = "response")
      eta_hat <- GA_T_ipw(x, T_obs, A, prob, tao, delta, G_c, G_c_tao)
    }
    eta_hat <- eta_hat/sum(abs(eta_hat))
    
    T_hat = ipw_T(eta_hat, x, T_obs, A, prob,tao,delta,G_c, G_c_tao)
    T_true <- eval_performance_g(eta_hat, gamma_p. = gamma_p, gamma_t. = gamma_t, dist = dist)[2]
    T_hat_g_opt <- ipw_T(gamma_t[4:6], x, T_obs, A, prob ,tao,delta,G_c, G_c_tao)
    T_true_g_opt <- eval_performance_g(gamma_t[4:6], gamma_p. = gamma_p, gamma_t. = gamma_t, dist = dist)[2]
    AA_given <- eval_performance_g(eta_hat, gamma_p. = gamma_p, gamma_t. = gamma_t,dist = dist)[4]
    
    eta <- rbind(eta, eta_hat)
    value <- rbind(value, c(T_hat, T_true, T_hat_g_opt, T_true_g_opt))
    AA <- c(AA, AA_given)
    colnames(value) <- c("V_hat(g_hat_opt)", "V(g_hat_opt)", "V_hat(g_opt)","V(g_opt)")
    
    detach(dat_obs)
    print(i)
  }
  list(eta, value,AA)
}
# r <- sim_T_ipw(n=1000, gamma_p=gamma_p, gamma_t=gamma_t, lambda_c=lambda_c, sim_num = 1, ps_model = TRUE) # test

# 1.ipw-T-true
registerDoRNG(2023);a0 = proc.time()[3]
result_ipwT_t = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  sim_T_ipw(n=1000, gamma_p=gamma_p, gamma_t=gamma_t, lambda_c=lambda_c,sim_num = 50, ps_model = TRUE, dist = "Weibull")
};print(proc.time()[3] - a0)
# 1.ipw-T-false
registerDoRNG(2023);a0 = proc.time()[3]
result_ipwT_f = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  sim_T_ipw(n=1000, gamma_p=gamma_p, gamma_t=gamma_t, lambda_c=lambda_c,sim_num = 50, ps_model = FALSE, dist = "Weibull")
};print(proc.time()[3] - a0)


result = break_result(result_ipwT_f)
cat(
  "eta", apply(result[[1]], 2, mean),"\n",
  "eta_sd", apply(result[[1]], 2, sd),"\n",
  "V", apply(result[[2]], 2, mean),"\n",
  "V_sd", apply(result[[2]], 2, sd),"\n",
  "AA:", mean(result[[3]]), "\n",
  "AA_sd", sd(result[[3]]))
# save as rds
saveRDS(result, "result\\Scenario_1\\Weibull\\ipw-T-t.RDS")
saveRDS(result, "result\\Scenario_1\\Weibull\\ipw-T-f.RDS")
load("result\\Scenario_1\\Weibull\\ipw_t.RDS")
########################## AIPW ################################
aipw_U <- function(eta, x, y, a, prob, tao, delta, G_c_tao ,p.es_0,p.es_1) {
  g <- as.numeric(x %*% eta > 0)
  wts <- a* g * 1 / prob + (1-a) * (1 - g) * (1 / (1 - prob)) # 逆概率权重
  y <- as.numeric(y > tao) # tao标记
  val_tao <- wts * y * (1 - delta)
  val_tao <- val_tao/G_c_tao
  Q <- g * p.es_1 + (1-g) * p.es_0
  val_mi <- (1-wts) * Q # 标准aipw
  
  val <- mean(val_tao + val_mi)
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
########################## AIPW-U estimator ################################
GA_U_aipw <- function(x, y, a, prob, tao, delta, G_c_tao, p.es_0,p.es_1) {
  d <- dim(x)[2]
  start_point <- rep(0, d)
  Domains <- cbind(rep(-1,d),rep(1,d))
  
  est<-genoud(fn=aipw_U, nvars=d, x=x,y=y,a=a,prob=prob,tao=tao,delta=delta,G_c_tao=G_c_tao,p.es_1=p.es_1,p.es_0=p.es_0,
              print.level=0,max=TRUE,pop.size=3000,wait.generations=10,gradient.check=FALSE,
              BFGS=FALSE,P1=50, P2=50, P3=10, P4=50, P5=50, P6=50, P7=50, P8=50, P9=0,Domains=Domains,starting.values=start_point,solution.tolerance=0.0001,optim.method="Nelder-Mead")
  eta <- est$par
  return(eta)
}
sim_U_aipw <- function(n = 500, gamma_p, gamma_t, lambda_c,sim_num = 100, ps_model=TRUE,augment_model=TRUE, dist = 'norm'){
  estv <- 0;mc_estv <- 0;truev <- 0;ratio <- c()
  eta <- c();value <- c();AA <- c()
  for (i in 1:sim_num) {
    dat <- sim_dat(n, gamma_p=gamma_p, gamma_t=gamma_t, lambda_c=lambda_c, dist = dist)
    dat_obs <- dat[[1]]
    dat_true <- dat[[2]]
    
    attach(dat_obs)
    x <- cbind(1, x1, x2)
    
    tao <- max(T_obs[delta == 1])
    fit_Gc <- survfit(Surv(T_obs, 1-delta)~1, data=dat_obs)
    surv <- stepfun(fit_Gc$time, c(1, fit_Gc$surv), right=FALSE, f=0)
    G_c <- surv(T_obs); G_c_tao <- surv(tao);
    G_c[G_c==0] <- min(G_c[G_c!=0])
    
    if (augment_model == TRUE){
      if (dist == 'Weibull') {
        dist0 <- "weibull"
      }
      if (dist == 'norm') {
        dist0 <- "lognormal"
      }
      fit <- gfcure(Surv(T_obs, delta) ~ x1+x2+x1*A+x2*A, ~ x1+x2+x1*A+x2*A, data= dat_obs, dist = dist0)
      gamma_p_hat <- -fit$coef[8:13]
      logit_p <- x%*%gamma_p_hat[1:3]; p.es_0 <- exp(logit_p)/(1+exp(logit_p)) 
      logit_p <- x%*%gamma_p_hat[1:3] + x%*%gamma_p_hat[4:6]; p.es_1 <- exp(logit_p)/(1+exp(logit_p))
    } else {
      # fit <- gfcure(Surv(T_obs, delta)~x1+x2+A, ~x1+x2+A, data= dat_obs, dist = "lognormal")
      # gamma_p_hat <- -c(fit$coef[6:9], 0, 0)
      # fit <- gfcure(Surv(T_obs, delta)~x1+x2+A, ~A, data= dat_obs, dist = "lognormal")
      # gamma_p_hat <- -c(fit$coef[6], 0, 0, fit$coef[7], 0,0)
      # fit <- gfcure(Surv(T_obs, delta)~x1+x2+A, ~x1, data= dat_obs, dist = "lognormal")
      # gamma_p_hat <- -c(fit$coef[6:7], 0, 0, 0,0)
      # logit_p <- x%*%gamma_p_hat[1:3]; p.es_0 <- exp(logit_p)/(1+exp(logit_p)) 
      # logit_p <- x%*%gamma_p_hat[1:3] + x%*%gamma_p_hat[4:6]; p.es_1 <- exp(logit_p)/(1+exp(logit_p))
      fit_f <- mixcure(Surv(T_obs, delta)~x1+x2+A, ~A, data=dat_obs, savedata = TRUE)
      dat_A1 <- dat_obs;dat_A1$A <- 1
      dat_A0 <- dat_obs;dat_A0$A <- 0
      result_1 = predict(fit_f, newdata=dat_A1,times=5)
      result_0 = predict(fit_f, newdata=dat_A0,times=5)
      p.es_1 <- result_1$cure[,1]
      p.es_0 <- result_0$cure[,1]
    }
    
    if (ps_model == TRUE){
      fit <- glm(A~I(x1^2) + I(x2^2), family = binomial(link = "logit"))
      # if (fit$converged == FALSE){break}
      prob <- predict(fit, type = "response")
      eta_hat <- GA_U_aipw(x, T_obs, A, prob, tao, delta, G_c_tao, p.es_0,p.es_1)
    } else {
      # propensity score错误给定
      fit <- glm(A~1, family = binomial(link = "logit"))
      prob <- predict(fit, type = "response")
      # prob= runif(n,0,1)
      eta_hat <- GA_U_aipw(x, T_obs, A, prob, tao, delta, G_c_tao, p.es_0,p.es_1)
    }
    eta_hat <- eta_hat/sum(abs(eta_hat)) # 一范数归一化
    
    estv = aipw_U(eta_hat, x, T_obs, A, prob,tao,delta,G_c_tao,p.es_0,p.es_1)
    mc_estv <- eval_performance_g(eta_hat, gamma_p. = gamma_p,gamma_t. = gamma_t, dist = dist)[1]
    V_hat_g_opt <- aipw_U(gamma_p[4:6], x, T_obs, A, prob ,tao,delta,G_c_tao,p.es_0,p.es_1)
    truev <- eval_performance_g(gamma_p[4:6], gamma_p. = gamma_p,gamma_t. = gamma_t, dist = dist)[1]
    AA_given <- eval_performance_g(eta_hat, gamma_p. = gamma_p,gamma_t. = gamma_t, dist = dist)[3]
    
    eta <- rbind(eta, eta_hat)
    value <- rbind(value, c(estv, mc_estv, V_hat_g_opt, truev))
    AA <- c(AA, AA_given)
    colnames(value) <- c("V_hat(g_hat_opt)", "V(g_hat_opt)", "V_hat(g_opt)","V(g_opt)")
    
    detach(dat_obs)
    print(i)
  }
  list(eta, value,AA)
}
# aipw-U-tt
registerDoRNG(2023);a0 = proc.time()[3]
result_aipwU_tt = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  sim_U_aipw(n=1000, gamma_p=gamma_p, gamma_t=gamma_t, lambda_c=lambda_c,sim_num = 50,ps_model = TRUE,
             augment_model = TRUE, dist="Weibull")
};print(proc.time()[3] - a0)

# aipw-U-tf
registerDoRNG(2023);a0 = proc.time()[3]
result_aipwU_tf = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  sim_U_aipw(n=1000, gamma_p=gamma_p, gamma_t=gamma_t, lambda_c=lambda_c, sim_num = 50,
             ps_model = TRUE,augment_model = FALSE,dist="Weibull")};print(proc.time()[3] - a0)
# aipw-U-ft
registerDoRNG(2023);a0 = proc.time()[3]
result_aipwU_ft = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  sim_U_aipw(n=1000, gamma_p=gamma_p, gamma_t=gamma_t, lambda_c=lambda_c, sim_num = 50,
             ps_model = FALSE,augment_model = TRUE,dist="Weibull")};print(proc.time()[3] - a0)
# aipw-U-ff
registerDoRNG(2023);a0 = proc.time()[3]
result_aipwU_ff = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  sim_U_aipw(n=1000, gamma_p=gamma_p, gamma_t=gamma_t, lambda_c=lambda_c, sim_num = 50,
             ps_model = FALSE,augment_model = FALSE, dist="Weibull")
};print(proc.time()[3] - a0)

result = break_result(result_aipwU_ff)
cat(
  "eta", apply(result[[1]], 2, mean),"\n",
  "eta_sd", apply(result[[1]], 2, sd),"\n",
  "V", apply(result[[2]], 2, mean),"\n",
  "V_sd", apply(result[[2]], 2, sd),"\n",
  "AA:", mean(result[[3]]), "\n",
  "AA_sd", sd(result[[3]]))
saveRDS(result, "result\\Scenario_1\\Weibull\\aipw-U-ff.RDS")

########################## AIPW-T estimator ################################
GA_T_aipw <- function(x, y, a, prob, tao, delta, G_c, G_c_tao,t.es_0,t.es_1, p.es_0, p.es_1) {
  d <- dim(x)[2]
  start_point <- rep(0, d)
  Domains <- cbind(rep(-1,d),rep(1,d))
  
  est<-genoud(fn=aipw_T, nvars=d, x=x,y=y,a=a,prob=prob,tao=tao,delta=delta, G_c=G_c, G_c_tao=G_c_tao,
              t.es_0=t.es_0,t.es_1=t.es_1,p.es_0=p.es_0,p.es_1=p.es_1,
              print.level=0,max=TRUE,pop.size=3000,wait.generations=10,gradient.check=FALSE,
              BFGS=FALSE,P1=50, P2=50, P3=10, P4=50, P5=50, P6=50, P7=50, P8=50, P9=0,Domains=Domains,starting.values=start_point,solution.tolerance=0.0001,optim.method="Nelder-Mead")
  eta <- est$par
  return(eta)
}
sim_T_aipw <- function(n = 1000, gamma_p=gamma_p, gamma_t=gamma_t, lambda_c,
                       sim_num = 100, ps_model=TRUE, augment_model=TRUE,dist = "norm"){
  T_hat<-c();T_true<-c();T_hat_g_opt<-c();T_true_g_opt<-c();
  eta <- c();value<-c();
  AA <- c()
  for (i in 1:sim_num) {
    dat <- sim_dat(n, gamma_p=gamma_p, gamma_t=gamma_t, lambda_c=lambda_c, dist=dist)
    dat_true <- dat[[2]]
    dat_obs <- dat[[1]]
    attach(dat_obs)
    x <- cbind(1, x1, x2)
    tao <- max(T_obs[delta == 1])
    fit_Gc <- survfit(Surv(T_obs, 1-delta)~1, data=dat_obs)
    surv <- stepfun(fit_Gc$time, c(1, fit_Gc$surv), right=FALSE, f=0)
    G_c <- surv(T_obs); G_c_tao <- surv(tao);
    G_c[G_c==0] <- min(G_c[G_c!=0])
    # propensity score正确给定
    fit <- glm(A~I(x1^2) + I(x2^2), family = binomial(link = "logit"))
    # if (fit$converged == FALSE){break}
    prob <- predict(fit, type = "response")

    if (augment_model == TRUE){
      if (dist == 'Weibull') {
        dist0 <- "weibull"
      }
      if (dist == 'norm') {
        dist0 <- "lognormal"
      }
      fit <- gfcure(Surv(T_obs, delta) ~ x1+x2+x1*A+x2*A, ~ x1+x2+x1*A+x2*A, data= dat_obs, dist = dist0)
      gamma_p_hat <- -fit$coef[8:13]
      logit_p <- x%*%gamma_p_hat[1:3]; p.es_0 <- exp(logit_p)/(1+exp(logit_p)) 
      logit_p <- x%*%gamma_p_hat[1:3] + x%*%gamma_p_hat[4:6]; p.es_1 <- exp(logit_p)/(1+exp(logit_p))
      gamma_t_hat <- fit$coef[2:7]
      t.es_0 <- exp(x%*%gamma_t_hat[1:3]) * exp(exp(fit$coef[1])^2/2)
      t.es_1 <- exp(x%*%gamma_t_hat[1:3] + 1* x %*%gamma_t_hat[4:6]) * exp(exp(fit$coef[1])^2/2)
    } else {
      fit_f <- mixcure(Surv(T_obs, delta)~x1+x2+A, ~A, data=dat_obs, savedata = TRUE)
      df_A1 <- dat_obs;df_A1$A <- 1
      df_A0 <- dat_obs;df_A0$A <- 0
      result_1 = predict(fit_f, newdata=df_A1,times=5)
      result_0 = predict(fit_f, newdata=df_A0,times=5)
      p.es_1 <- result_1$cure[,1]
      p.es_0 <- result_0$cure[,1]
      fit_latency_1 = survfit(fit_f$lfit$fit, newdata = df_A1, type = "kaplan-meier")
      fit_latency_0 = survfit(fit_f$lfit$fit, newdata = df_A0, type = "kaplan-meier")
      t.es_1 <-survival:::survmean(fit_latency_1, rmean = 100)$matrix[,'rmean']
      t.es_0 <-survival:::survmean(fit_latency_0, rmean=100)$matrix[,'rmean']
    }
    if (ps_model == TRUE){
      eta_hat <- GA_T_aipw(x, T_obs, A, prob, tao, delta, G_c, G_c_tao, t.es_0, t.es_1, p.es_0, p.es_1)
    } else {
      # prob = mean(dat_obs$A)
      # prob <- rep(0.25, n)
      fit <- glm(A~1, family = binomial(link = "logit"))
      prob <- predict(fit, type = "response")
      eta_hat <- GA_T_aipw(x, T_obs, A, prob, tao, delta, G_c, G_c_tao, t.es_0, t.es_1, p.es_0, p.es_1)
    }
    eta_hat <- eta_hat/sum(abs(eta_hat))
    
    T_hat = aipw_T(eta_hat, x, T_obs, A, prob,tao,delta,G_c, G_c_tao,t.es_0, t.es_1, p.es_0, p.es_1)
    T_true <- eval_performance_g(eta_hat, gamma_p=gamma_p, gamma_t=gamma_t,dist=dist)[2]
    T_hat_g_opt <- aipw_T(gamma_t[4:6], x, T_obs, A, prob ,tao,delta,G_c, G_c_tao,t.es_0, t.es_1, p.es_0, p.es_1)
    T_true_g_opt <- eval_performance_g(gamma_t[4:6], gamma_p=gamma_p, gamma_t=gamma_t,dist=dist)[2]
    AA_given <- eval_performance_g(eta_hat, gamma_p=gamma_p, gamma_t=gamma_t,dist=dist)[3]
    
    eta <- rbind(eta, eta_hat)
    value <- rbind(value, c(T_hat, T_true, T_hat_g_opt, T_true_g_opt))
    AA <- c(AA, AA_given)
    colnames(value) <- c("V_hat(g_hat_opt)", "V(g_hat_opt)", "V_hat(g_opt)","V(g_opt)")
    
    detach(dat_obs)
    print(i)
  }
  list(eta, value,AA)
}
# sim_T_aipw(n=1000, theta=theta, lambda_c=lambda_c,tao_method = "max", sim_num = 1, ps_model = TRUE)
# aipw-T-tt
registerDoRNG(2023);a0 = proc.time()[3]
result_aipwT_tt = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  sim_T_aipw(n=1000, gamma_p=gamma_p, gamma_t=gamma_t, lambda_c=lambda_c, sim_num = 50,
             ps_model = TRUE,augment_model = TRUE,dist="Weibull")};print(proc.time()[3] - a0)
# aipw-T-tf
registerDoRNG(2023)
result_aipwT_tf = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  sim_T_aipw(n=1000, gamma_p=gamma_p, gamma_t=gamma_t  , lambda_c=lambda_c, sim_num = 50,
             ps_model = TRUE,augment_model = FALSE,dist="Weibull")}
# aipw-T-ft
registerDoRNG(2023);a0 = proc.time()[3]
result_aipwT_ft = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  sim_T_aipw(n=1000, gamma_p=gamma_p, gamma_t=gamma_t, lambda_c=lambda_c, sim_num = 50,
             ps_model = FALSE,augment_model = TRUE,dist="Weibull")};print(proc.time()[3] - a0)
# aipw-T-ff
registerDoRNG(2023);a0 = proc.time()[3]
result_aipwT_ff = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  sim_T_aipw(n=1000, gamma_p=gamma_p, gamma_t=gamma_t, lambda_c=lambda_c, sim_num = 50,
             ps_model = FALSE,augment_model = FALSE,dist="Weibull")};print(proc.time()[3] - a0)

result = break_result(result_aipwT_tf)
cat(
  "eta", apply(result[[1]], 2, mean),"\n",
  "eta_sd", apply(result[[1]], 2, sd),"\n",
  "V", apply(result[[2]], 2, mean),"\n",
  "V_sd", apply(result[[2]], 2, sd),"\n",
  "AA:", mean(result[[3]]), "\n",
  "AA_sd", sd(result[[3]]))
saveRDS(result, "result\\Scenario_1\\Weibull\\aipw-T-ff.RDS")
