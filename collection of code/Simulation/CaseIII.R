########################## package loading ################################
library(mixcure);library(rgenoud);library(gfcure)
library(parallel);library(doParallel);library(foreach);library(doRNG)
cores_used <- 10 ;cl <- makeCluster(cores_used);registerDoParallel(cl)
clusterEvalQ(cl, library(rgenoud));clusterEvalQ(cl, library(mixcure));clusterEvalQ(cl, library(gfcure))
# stopCluster(cl)
########################## true value given treatment(eta) #######################
eval_performance_g_stoc <- function(beta, sim_num = 1e5, gamma_p., gamma_t.) {
  x1 <- runif(sim_num,-1.5,1.5)
  x2 <- runif(sim_num,-1,1)
  x <- cbind(1, x1, x2)
  g <- 1 / (1 + exp(x %*% beta))
  A = rbinom(dim(x)[1],size=1,prob=1/(1+exp(x%*%beta)))
  
  logit_p <- x %*% gamma_p[1:3] + A * x %*% gamma_p[4:6]
  p <- exp(logit_p) / (1 + exp(logit_p))
  cure_status_true <- rbinom(sim_num, 1, p)
  EU = mean(cure_status_true)
  
  z <- x
  epsilon <- rnorm(sim_num - sum(cure_status_true), 0 ,1)
  T <- exp(z[cure_status_true==0,]%*%gamma_t[1:3] + A[cure_status_true==0]*z[cure_status_true==0,]%*%gamma_t[4:6]+epsilon)
  ET = mean(T)
  
  result <- c(EU, ET)
  return(result)
}
########################## parameters setting ################################
gamma_p = c(c(-0.75,0.25,0.25),c(0,1.5,-1.5));gamma_t=c(c(0,0.25,0.25),c(0,-1,1))
lambda_c = 0.03
set.seed(2023)
EDA(gamma_p, gamma_t, lambda_c)
################################### stochastic regime & value ###################################
# IPW
ipw_U_stoc <- function(beta, x, y, a, prob, tao, delta, G_c_tao) {
  #   g <- as.numeric(x %*% eta > 0)
  g <- 1 / (1 + exp(x %*% beta))
  wts <- a * g * 1 / prob + (1 - a) * (1 - g) * (1 / (1 - prob)) # 逆概率权重
  
  y <- as.numeric(y > tao) # tao标记
  val_tao <- wts * y * (1 - delta)
  val_tao <- val_tao / G_c_tao
  
  val <- mean(val_tao)
  return(val)
}
ipw_T_stoc <- function(beta, x, y, a, prob, tao, delta, G_c, G_c_tao) {
  #   g <- as.numeric(x %*% eta < 0)
  g <- 1 / (1 + exp(x %*% beta))
  wts <- a * g * 1 / prob + (1 - a) * (1 - g) * (1 / (1 - prob)) # 逆概率权重
  v_t <- mean(wts * y * delta / G_c)
  v_u <- ipw_U_stoc(beta, x, y, a, prob, tao, delta, G_c_tao)
  val_t <- v_t / (1 - v_u)
  
  val <- mean(val_t)
  return(val)
}
cptgrd_U <- function(beta, x, y, a, prob, tao, delta, G_c_tao) {
  d <- length(beta)
  m0 <- ipw_U_stoc(beta, x, y, a, prob, tao, delta, G_c_tao)
  grd <- rep(0, d)
  for (j in 1:d) {
    tmp <- rep(0, d)
    tmp[j] <- 0.01 * (2 * rbinom(1, size = 1, prob = 0.5) - 1)
    grd[j] <- (ipw_U_stoc(beta+tmp, x, y, a, prob, tao, delta, G_c_tao) - m0) / tmp[j] # ???ţ?
  }
  return(grd)
}
cptgrd_T <- function(beta, x, y, a, prob, tao, delta, G_c, G_c_tao) {
  d <- length(beta)
  m0 <- ipw_T_stoc(beta, x, y, a, prob, tao, delta, G_c, G_c_tao)
  grd <- rep(0, d)
  for (j in 1:d) {
    tmp <- rep(0, d)
    tmp[j] <- 0.01 * (2 * rbinom(1, size = 1, prob = 0.5) - 1)
    grd[j] <- (ipw_T_stoc(beta+tmp, x, y, a, prob, tao, delta, G_c, G_c_tao) - m0) / tmp[j] # ???ţ?
  }
  return(grd)
}
# AIPW
aipw_U_stoc <- function(beta, x, y, a, prob, tao, delta, G_c_tao, p.es_0, p.es_1) {
  g <- 1 / (1 + exp(x %*% beta))
  wts <- a * g * 1 / prob + (1 - a) * (1 - g) * (1 / (1 - prob)) # 逆概率权重
  
  y <- as.numeric(y > tao) # tao标记
  val_tao <- wts * y * (1 - delta)
  val_tao <- val_tao / G_c_tao
  
  p.es <- g * p.es_1 + (1 - g) * p.es_0
  val_mi <- (wts - 1) * p.es # 标准aipw
  
  val <- mean(val_tao - val_mi)
  return(val)
}
aipw_T_stoc <- function(beta, x, y, a, prob, tao, delta, G_c, G_c_tao, t.es_0, t.es_1, p.es_0, p.es_1) {
  g <- 1 / (1 + exp(x %*% beta))
  wts <- a * g * 1 / prob + (1 - a) * (1 - g) * (1 / (1 - prob))
  H <- g * t.es_1 * p.es_0 + (1 - g) * t.es_0 * p.es_1
  v_t <- mean(wts * y * delta / G_c + (1-wts) * H)
  
  v_u <- aipw_U_stoc(beta, x, y, a, prob, tao, delta, G_c_tao, p.es_0, p.es_1)
  val_t <- v_t / (1 - v_u)
  return(val_t)
}
cptgrd_U_aipw <- function(beta, x, y, a, prob, tao, delta, G_c_tao, p.es_0, p.es_1) {
  d <- length(beta)
  m0 <- aipw_U_stoc(beta, x, y, a, prob, tao, delta, G_c_tao, p.es_0, p.es_1)
  grd <- rep(0, d)
  for (j in 1:d) {
    tmp <- rep(0, d)
    tmp[j] <- 0.01 * (2 * rbinom(1, size = 1, prob = 0.5) - 1)
    grd[j] <- (aipw_U_stoc(beta+tmp, x, y, a, prob, tao, delta, G_c_tao, p.es_0, p.es_1) - m0) / tmp[j] # ???ţ?
  }
  return(grd)
}
cptgrd_T_aipw <- function(beta, x, y, a, prob, tao, delta, G_c, G_c_tao, t.es_0, t.es_1, p.es_0, p.es_1) {
  d <- length(beta)
  m0 <- aipw_T_stoc(beta, x, y, a, prob, tao, delta, G_c, G_c_tao, t.es_0, t.es_1, p.es_0, p.es_1)
  grd <- rep(0, d)
  for (j in 1:d) {
    tmp <- rep(0, d)
    tmp[j] <- 0.01 * (2 * rbinom(1, size = 1, prob = 0.5) - 1)
    grd[j] <- (aipw_T_stoc(beta+tmp, x, y, a, prob, tao, delta, G_c, G_c_tao, t.es_0, t.es_1, p.es_0, p.es_1) - m0) / tmp[j] # ???ţ?
  }
  return(grd)
}
############################## 1:max E(T|U=0), s.t. E(U)>=q ##############################
########################## IPW: algorithm 1-dual ################################
# dual function
dual_1 <- function(beta,lambda,q, x, y, a, prob, tao, delta,G_c, G_c_tao) {
  val <- -ipw_T_stoc(beta, x, y, a, prob, tao, delta, G_c, G_c_tao) +
    lambda * (q - ipw_U_stoc(beta, x, y, a, prob, tao, delta, G_c_tao))
  return(val)
}
# dual function gradient
cptgrd_d1 <- function(beta,lambda,q, x, y, a, prob, tao, delta, G_c, G_c_tao) {
  grd_u <- cptgrd_U(beta, x, y, a, prob, tao, delta, G_c_tao)
  grd_t <- cptgrd_T(beta, x, y, a, prob, tao, delta, G_c, G_c_tao)
  grd <- -grd_t - lambda * grd_u
  return(grd)
}
solve_d1 <- function(lambda, q, x, y, a, prob, tao, delta, G_c, G_c_tao,
                     stp_sz = 0.01, tol = 1e-3, max_iter = 10000, result="value") {
  beta <- rep(0, length(x[1,])) 
  iter <- 0 
  err <- 999
  while ((err > tol) && (iter < max_iter)) {
    grd <- cptgrd_d1(beta,lambda,q, x, y, a, prob, tao, delta, G_c, G_c_tao)
    beta <- beta - stp_sz * grd
    iter <- iter + 1
    err <- sum(abs(grd))
    # print(iter)
  }
  # return(beta) # argmax
  if (result == "value") {
    return(dual_1(beta,lambda,q, x, y, a, prob, tao, delta, G_c, G_c_tao))
  } else {
    return(beta)
  }
}
ipw_T_d1 <- function(q_cons, n = 500, gamma_p, gamma_t, lambda_c, sim_num = 100){
  EU_hat<-0;ET_hat<-0;EU<-0;ET<-0;
  beta <- c();value <- c();lambda<-c()
  for (i in 1:sim_num) {
    dat <- sim_dat(n, gamma_p=gamma_p, gamma_t=gamma_t, lambda_c=lambda_c)
    dat_obs <- dat[[1]]
    dat_true <- dat[[2]]
    
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

    lambda_opt <- optimize(solve_d1, interval = c(0, 100), q = q_cons,
                           x = x, y = T_obs, a = A, prob = prob, tao = tao, delta = delta, G_c=G_c , G_c_tao=G_c_tao,
                           maximum = TRUE, tol = 0.1)$maximum
    # print(lambda_opt)
    beta_hat <- solve_d1(lambda_opt, q = q_cons, x = x, y = T_obs, a = A, prob = prob, tao = tao, delta = delta, G_c=G_c ,G_c_tao=G_c_tao,
                         result = "beta")
    # 不再标准化系数
    
    # 以下均为对偶最优解，beta_hat的value_hat & value
    EU_hat <- ipw_U_stoc(beta_hat, x, T_obs, A, prob,tao,delta,G_c_tao)
    ET_hat <- ipw_T_stoc(beta_hat, x, T_obs, A, prob,tao,delta,G_c, G_c_tao)
    temp_eval <- eval_performance_g_stoc(beta_hat, gamma_p. = gamma_p,gamma_t. = gamma_t)
    EU <- temp_eval[1]
    ET <- temp_eval[2]
    
    beta <- rbind(beta, beta_hat)
    value <- rbind(value, c(EU_hat,EU,ET_hat,ET))
    lambda <- rbind(lambda, lambda_opt)
    colnames(value) <- c("EU_hat", "EU", "ET_hat","ET")
    
    detach(dat_obs)
    print(i)
  }
  list(beta, value, lambda)
}
ipw_T_d1(q_cons=0.40, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 1)
registerDoRNG(2023);a0 = proc.time()[3]
result_ipwT_d1_38 = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  ipw_T_d1(q_cons=0.38, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 25)
};print(proc.time()[3] - a0)
registerDoRNG(2023);a0 = proc.time()[3]
result_ipwT_d1_40 = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  ipw_T_d1(q_cons=0.40, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 25)
};print(proc.time()[3] - a0)
registerDoRNG(2023);a0 = proc.time()[3]
result_ipwT_d1_42 = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  ipw_T_d1(q_cons=0.42, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 25)
};print(proc.time()[3] - a0)
registerDoRNG(2023);a0 = proc.time()[3]
result_ipwT_d1_44 = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  ipw_T_d1(q_cons=0.44, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 25)
};print(proc.time()[3] - a0)
registerDoRNG(2023);a0 = proc.time()[3]
result_ipwT_d1_46 = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  ipw_T_d1(q_cons=0.46, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 25)
};print(proc.time()[3] - a0)

result = break_result(result_ipwT_d1_46)
cat(
  "eta", apply(result[[1]], 2, mean),"\n",
  "eta_sd", apply(result[[1]], 2, sd),"\n",
  "V", apply(result[[2]], 2, mean),"\n",
  "V_sd", apply(result[[2]], 2, sd),"\n",
  "AA:", mean(result[[3]]), "\n",
  "AA_sd", sd(result[[3]]))
saveRDS(result, "D:\\论文\\Cure Model + ITR\\模拟\\复现整理版 - 副本\\T_ipw_d1_46.RDS")
########################## IPW: algorithm 2-penalty ################################
cptgrd_p1 <- function(beta,lambda,q, x, y, a, prob, tao, delta,G_c, G_c_tao){
  grd_t <- cptgrd_T(beta, x, y, a, prob, tao, delta, G_c, G_c_tao)
  # 罚函数分段求导
  if (ipw_U_stoc(beta, x, T_obs, A, prob,tao,delta,G_c_tao)>=q){
    grd_u = 0
  } else{
    grd_u <- cptgrd_U(beta, x, y, a, prob, tao, delta, G_c_tao)
  }
  grd = -grd_t - grd_u*lambda # 正负号（原始代码是负的。感觉是错的）
  return(grd)
}
solve_p1 <- function(lambda, q, x, y, a, prob, tao, delta, G_c, G_c_tao,
                     stp_sz = 0.01, tol = 1e-3, max_iter = 10000, result="beta") {
  beta <- rep(0, length(x[1,])) # 初始点设为原点
  iter <- 0;err <- 999
  while ((err > tol) && (iter < max_iter)) {
    grd <- cptgrd_p1(beta,lambda,q, x, y, a, prob, tao, delta, G_c, G_c_tao)
    beta <- beta - stp_sz * grd
    iter <- iter + 1
    err <- sum(abs(grd))
    # print(iter)
  }
  # return(beta) # argmax
  if (result == "value") {
    return(dual_1(beta,lambda,q, x, y, a, prob, tao, delta, G_c, G_c_tao))
  } else {
    return(beta)
  }
}
MaxT_ipw_p1 <- function(q_cons, n = 500, gamma_p, gamma_t, lambda_c,sim_num = 100){
  EU_hat<-0;ET_hat<-0;EU<-0;ET<-0;
  beta <- c();value <- c();lambda<-c()
  for (i in 1:sim_num) {
    dat <- sim_dat(n, gamma_p=gamma_p, gamma_t=gamma_t, lambda_c=lambda_c)
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

    # 罚函数
    lambda0 = 2
    c=1.5
    # lambda_seq = 0.25*lambda0*c^(0:13)
    # lambda_seq = 0.25*lambda0*c^(0:13)
    lambda_seq = c(0.1, 0.5, 1, 2, 8, 32, 64, 128, 256, 512)
    i=0; EU_hat=0
    while((EU_hat<q_cons) & (i<=9)){
      i=i+1
      lambda_temp = lambda_seq[i]
      beta_hat  = solve_p1(lambda_temp, q = q_cons, x = x, y = T_obs, a = A, prob = prob, tao = tao, delta = delta, 
                           G_c=G_c, G_c_tao = G_c_tao,result = "beta")
      EU_hat <- ipw_U_stoc(beta_hat, x, T_obs, A, prob,tao,delta,G_c_tao)
    }
    # 不再标准化系数
    
    # 以下均为对偶最优解，beta_hat的value_hat & value
    EU_hat <- ipw_U_stoc(beta_hat, x, T_obs, A, prob,tao,delta,G_c_tao)
    ET_hat <- ipw_T_stoc(beta_hat, x, T_obs, A, prob,tao,delta,G_c, G_c_tao)
    temp_eval <- eval_performance_g_stoc(beta_hat, gamma_p. = gamma_p,gamma_t. = gamma_t)
    EU <- temp_eval[1]
    ET <- temp_eval[2]
    
    beta <- rbind(beta, beta_hat)
    value <- rbind(value, c(EU_hat,EU,ET_hat,ET))
    lambda <- rbind(lambda, lambda_temp)
    colnames(value) <- c("EU_hat", "EU", "ET_hat","ET")
    
    detach(dat_obs)
    print(i)
  }
  list(beta, value, lambda)
}
MaxT_ipw_p1(q_cons=0.01, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 1)
registerDoRNG(2023);a0 = proc.time()[3]
result_T_ipw_p1_0 = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  MaxT_ipw_p1(q_cons=1e-5, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 10)};print(proc.time()[3] - a0)

registerDoRNG(2023);a0 = proc.time()[3]
result_T_ipw_p1_38 = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  MaxT_ipw_p1(q_cons=0.38, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 25)};print(proc.time()[3] - a0)

registerDoRNG(2023);a0 = proc.time()[3]
result_T_ipw_p1_40 = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  MaxT_ipw_p1(q_cons=0.40, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 25)};print(proc.time()[3] - a0)

registerDoRNG(2023);a0 = proc.time()[3]
result_T_ipw_p1_42 = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  MaxT_ipw_p1(q_cons=0.42, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 25)};print(proc.time()[3] - a0)

registerDoRNG(2023);a0 = proc.time()[3]
result_T_ipw_p1_44 = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  MaxT_ipw_p1(q_cons=0.44, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 25)};print(proc.time()[3] - a0)

registerDoRNG(2023);a0 = proc.time()[3]
result_T_ipw_p1_46 = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  MaxT_ipw_p1(q_cons=0.46, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 25)};print(proc.time()[3] - a0)

result = break_result(result_T_ipw_p1_46)
cat(
  "beta", apply(result[[1]], 2, mean),"\n",
  "beta_sd", apply(result[[1]], 2, sd),"\n",
  "V", apply(result[[2]], 2, mean),"\n",
  "V_sd", apply(result[[2]], 2, sd),"\n",
  "lambda:", mean(result[[3]]))
saveRDS(result, "D:\\论文\\Cure Model + ITR\\模拟\\复现整理版 - 副本\\T_ipw_p1_46.RDS")
########################## AIPW: algorithm 1-dual ################################
# dual function
dual_aipw_1 <- function(beta,lambda,q, x, y, a, prob, tao, delta, G_c, G_c_tao, t.es_0, t.es_1, p.es_0, p.es_1) {
  val <- -aipw_T_stoc(beta, x, y, a, prob, tao, delta, G_c, G_c_tao, t.es_0, t.es_1, p.es_0, p.es_1) +
    lambda * (q - aipw_U_stoc(beta, x, y, a, prob, tao, delta, G_c_tao , p.es_0, p.es_1))
  return(val)
}
# dual function gradient
cptgrd_aipw_d1 <- function(beta,lambda,q, x, y, a, prob, tao, delta, G_c, G_c_tao, t.es_0, t.es_1, p.es_0, p.es_1) {
  grd_u <- cptgrd_U_aipw(beta, x, y, a, prob, tao, delta, G_c_tao, p.es_0, p.es_1)
  grd_t <- cptgrd_T_aipw(beta, x, y, a, prob, tao, delta, G_c, G_c_tao, t.es_0, t.es_1, p.es_0, p.es_1)
  grd <- -grd_t - lambda * grd_u
  return(grd)
}
solve_aipw_d1 <- function(lambda, q, x, y, a, prob, tao, delta, G_c, G_c_tao, t.es_0, t.es_1, p.es_0, p.es_1,
                     stp_sz = 0.01, tol = 1e-3, max_iter = 10000, result="value") {
  beta <- rep(0, length(x[1,]))
  iter <- 0 
  err <- 999 
  while ((err > tol) && (iter < max_iter)) {
    grd <- cptgrd_aipw_d1(beta,lambda,q, x, y, a, prob, tao, delta, G_c, G_c_tao, t.es_0, t.es_1, p.es_0, p.es_1)
    beta <- beta - stp_sz * grd
    iter <- iter + 1
    err <- sum(abs(grd))
    # print(iter)
  }
  # return(beta) # argmax
  if (result == "value") {
    return(dual_aipw_1(beta,lambda,q, x, y, a, prob, tao, delta, G_c, G_c_tao, t.es_0, t.es_1, p.es_0, p.es_1))
  } else {
    return(beta)
  }
}
aipw_T_d1 <- function(q_cons, n = 500, gamma_p, gamma_t, lambda_c, sim_num = 100){
  EU_hat<-0;ET_hat<-0;EU<-0;ET<-0;
  beta <- c();value <- c();lambda<-c()
  for (i in 1:sim_num) {
    dat <- sim_dat(n, gamma_p=gamma_p, gamma_t=gamma_t, lambda_c=lambda_c)
    dat_obs <- dat[[1]]
    dat_true <- dat[[2]]
    
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
    
    lambda_opt <- optimize(solve_aipw_d1, interval = c(0, 100), q = q_cons,
                           x = x, y = T_obs, a = A, prob = prob, tao = tao, delta = delta, G_c=G_c , G_c_tao=G_c_tao,
                           t.es_0=t.es_0, t.es_1=t.es_1, p.es_0=p.es_0, p.es_1=p.es_1,maximum = TRUE, tol = 0.1)$maximum
    # print(lambda_opt)
    beta_hat <- solve_aipw_d1(lambda_opt, q = q_cons, x = x, y = T_obs, a = A, prob = prob, tao = tao, delta = delta, G_c=G_c, G_c_tao=G_c_tao,
                              t.es_0=t.es_0, t.es_1=t.es_1, p.es_0=p.es_0, p.es_1=p.es_1,result = "beta")
    # 不再标准化系数
    
    # 以下均为对偶最优解，beta_hat的value_hat & value
    EU_hat <- aipw_U_stoc(beta_hat, x, T_obs, A, prob,tao,delta,G_c_tao, p.es_0, p.es_1)
    ET_hat <- aipw_T_stoc(beta_hat, x, T_obs, A, prob,tao,delta,G_c, G_c_tao, t.es_0, t.es_1, p.es_0, p.es_1)
    temp_eval <- eval_performance_g_stoc(beta_hat, gamma_p. = gamma_p,gamma_t. = gamma_t)
    EU <- temp_eval[1]
    ET <- temp_eval[2]
    
    beta <- rbind(beta, beta_hat)
    value <- rbind(value, c(EU_hat,EU,ET_hat,ET))
    lambda <- rbind(lambda, lambda_opt)
    colnames(value) <- c("EU_hat", "EU", "ET_hat","ET")
    
    detach(dat_obs)
    print(i)
  }
  list(beta, value, lambda)
}


registerDoRNG(2023);a0 = proc.time()[3]
result_aipwT_d1_38 = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  aipw_T_d1(q_cons=0.38, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 25)};print(proc.time()[3] - a0)

registerDoRNG(2023);a0 = proc.time()[3]
result_aipwT_d1_40 = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  aipw_T_d1(q_cons=0.40, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 25)};print(proc.time()[3] - a0)

registerDoRNG(2023);a0 = proc.time()[3]
result_aipwT_d1_42 = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  aipw_T_d1(q_cons=0.42, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 25)};print(proc.time()[3] - a0)

registerDoRNG(2023);a0 = proc.time()[3]
result_aipwT_d1_44 = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  aipw_T_d1(q_cons=0.44, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 25)};print(proc.time()[3] - a0)

registerDoRNG(2023);a0 = proc.time()[3]
result_aipwT_d1_46 = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  aipw_T_d1(q_cons=0.46, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 25)};print(proc.time()[3] - a0)

result = break_result(result_aipwT_d1_46)
cat(
  "eta", apply(result[[1]], 2, mean),"\n",
  "eta_sd", apply(result[[1]], 2, sd),"\n",
  "V", apply(result[[2]], 2, mean),"\n",
  "V_sd", apply(result[[2]], 2, sd),"\n",
  "AA:", mean(result[[3]]), "\n",
  "AA_sd", sd(result[[3]]))
saveRDS(result, "D:\\论文\\Cure Model + ITR\\模拟\\复现整理版 - 副本\\T_aipw_d1_46.RDS")
########################## AIPW: algorithm 2-penalty ################################
cptgrd_p1 <- function(beta,lambda,q, x, y, a, prob, tao, delta, G_c, G_c_tao, t.es_0, t.es_1, p.es_0, p.es_1){
  grd_t <- cptgrd_T_aipw(beta, x, y, a, prob, tao, delta, G_c, G_c_tao, t.es_0, t.es_1, p.es_0, p.es_1)
  # 罚函数分段求导
  if (aipw_U_stoc(beta, x, T_obs, A, prob,tao,delta,G_c_tao,p.es_0,p.es_1)>=q){
    grd_u = 0
  } else{
    grd_u <- cptgrd_U_aipw(beta, x, y, a, prob, tao, delta, G_c_tao, p.es_0, p.es_1)
  }
  grd = -grd_t - grd_u*lambda # 正负号（原始代码是负的。感觉是错的）
  return(grd)
}
solve_p1 <- function(lambda, q, x, y, a, prob, tao, delta, G_c, G_c_tao, t.es_0, t.es_1, p.es_0, p.es_1,
                     stp_sz = 0.01, tol = 1e-3, max_iter = 10000, result="beta") {
  beta <- rep(0, length(x[1,])) # 初始点设为原点
  iter <- 0
  err <- 999
  while ((err > tol) && (iter < max_iter)) {
    grd <- cptgrd_p1(beta,lambda,q, x, y, a, prob, tao, delta, G_c, G_c_tao, t.es_0, t.es_1, p.es_0, p.es_1)
    beta <- beta - stp_sz * grd
    iter <- iter + 1
    err <- sum(abs(grd))
    # print(iter)
  }
  # return(beta) # argmax
  if (result == "value") {
    return(dual_1(beta,lambda,q, x, y, a, prob, tao, delta, G_c, G_c_tao, t.es_0, t.es_1, p.es_0, p.es_1))
  } else {
    return(beta)
  }
}
MaxT_aipw_p1 <- function(q_cons, n = 500, gamma_p, gamma_t, lambda_c,sim_num = 100){
  EU_hat<-0;ET_hat<-0;EU<-0;ET<-0;
  beta <- c();value <- c();lambda<-c()
  for (i in 1:sim_num) {
    dat <- sim_dat(n, gamma_p=gamma_p, gamma_t=gamma_t, lambda_c=lambda_c)
    dat_obs <- dat[[1]]
    dat_true <- dat[[2]]
    
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
    
    fit <- gfcure(Surv(T_obs, delta) ~ x1+x2+x1*A+x2*A, ~ x1+x2+x1*A+x2*A, data= dat_obs, dist = "lognormal")
    gamma_p_hat <- -fit$coef[8:13]
    logit_p <- x%*%gamma_p_hat[1:3]; p.es_0 <- exp(logit_p)/(1+exp(logit_p)) 
    logit_p <- x%*%gamma_p_hat[1:3] + x%*%gamma_p_hat[4:6]; p.es_1 <- exp(logit_p)/(1+exp(logit_p))
    gamma_t_hat <- fit$coef[2:7]
    t.es_0 <- exp(x%*%gamma_t_hat[1:3]) * exp(exp(fit$coef[1])^2/2)
    t.es_1 <- exp(x%*%gamma_t_hat[1:3] + 1* x %*%gamma_t_hat[4:6]) * exp(exp(fit$coef[1])^2/2)

    # 罚函数
    lambda0 = 2
    c=1.5
    lambda_seq = 0.25*lambda0*c^(0:14)
    i=0; EU_hat=0
    while((EU_hat<q_cons) & (i<=14)){
      i=i+1
      lambda_temp = lambda_seq[i]
      beta_hat  = solve_p1(lambda_temp, q = q_cons, x = x, y = T_obs, a = A, prob = prob, tao = tao, delta = delta, 
                           G_c=G_c, G_c_tao=G_c_tao , t.es_0 = t.es_0, t.es_1 = t.es_1, p.es_0 = p.es_0, p.es_1 = p.es_1, 
                           result = "beta")
      EU_hat <- aipw_U_stoc(beta_hat, x, T_obs, A, prob,tao,delta, G_c_tao, p.es_0, p.es_1)
    }
    # 不再标准化系数
    
    # 以下均为对偶最优解，beta_hat的value_hat & value
    EU_hat <- aipw_U_stoc(beta_hat, x, T_obs, A, prob,tao, delta, G_c_tao ,p.es_0,p.es_1)
    ET_hat <- aipw_T_stoc(beta_hat, x, T_obs, A, prob,tao,delta , G_c, G_c_tao, t.es_0, t.es_1, p.es_0, p.es_1)
    temp_eval <- eval_performance_g_stoc(beta_hat, gamma_p. = gamma_p,gamma_t. = gamma_t)
    EU <- temp_eval[1]
    ET <- temp_eval[2]
    
    beta <- rbind(beta, beta_hat)
    value <- rbind(value, c(EU_hat,EU,ET_hat,ET))
    lambda <- rbind(lambda, lambda_temp)
    colnames(value) <- c("EU_hat", "EU", "ET_hat","ET")
    
    detach(dat_obs)
    print(i)
  }
  list(beta, value, lambda)
}
# MaxT_aipw_p1(q_cons=0.40, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 1)
registerDoRNG(2023);a0 = proc.time()[3]
result_T_aipw_p1_38 = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  MaxT_aipw_p1(q_cons=0.38, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 25)};print(proc.time()[3] - a0)

registerDoRNG(2023);a0 = proc.time()[3]
result_T_aipw_p1_40 = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  MaxT_aipw_p1(q_cons=0.40, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 25)};print(proc.time()[3] - a0)

registerDoRNG(2023);a0 = proc.time()[3]
result_T_aipw_p1_42 = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  MaxT_aipw_p1(q_cons=0.42, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 25)};print(proc.time()[3] - a0)

registerDoRNG(2023);a0 = proc.time()[3]
result_T_aipw_p1_44 = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  MaxT_aipw_p1(q_cons=0.44, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 25)};print(proc.time()[3] - a0)

registerDoRNG(2023);a0 = proc.time()[3]
result_T_aipw_p1_46 = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  MaxT_aipw_p1(q_cons=0.46, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 25)};print(proc.time()[3] - a0)

result = break_result(result_T_aipw_p1_46)
cat(
  "beta", apply(result[[1]], 2, mean),"\n",
  "beta_sd", apply(result[[1]], 2, sd),"\n",
  "V", apply(result[[2]], 2, mean),"\n",
  "V_sd", apply(result[[2]], 2, sd),"\n",
  "lambda:", mean(result[[3]]))
saveRDS(result, "D:\\论文\\Cure Model + ITR\\模拟\\复现整理版 - 副本\\T_aipw_p1_46.RDS")



############################## 2: max E(U), s.t. E(T|U=0)>=q##############################
########################## IPW: algorithm 1-dual ################################
# dual function
dual_2 <- function(beta,lambda,q, x, y, a, prob, tao, delta, G_c, G_c_tao) {
  val <- -ipw_U_stoc(beta, x, y, a, prob, tao, delta, G_c_tao) +
    lambda * (q - ipw_T_stoc(beta, x, y, a, prob, tao, delta, G_c, G_c_tao))
  return(val)
}
# dual function gradient
cptgrd_d2 <- function(beta,lambda,q, x, y, a, prob, tao, delta, G_c, G_c_tao) {
  grd_u <- cptgrd_U(beta, x, y, a, prob, tao, delta, G_c_tao)
  grd_t <- cptgrd_T(beta, x, y, a, prob, tao, delta, G_c, G_c_tao)
  grd <- -grd_u - lambda * grd_t
  return(grd)
}
solve_d2 <- function(lambda, q, x, y, a, prob, tao, delta, G_c, G_c_tao,
                     stp_sz = 0.01, tol = 1e-3, max_iter = 10000, result="value") {
  beta <- rep(0, length(x[1,])) # 初始点设为原点
  iter <- 0 # 迭代次数
  err <- 999 # 迭代误差
  # stp_sz: 步长
  # tol:误差容忍度
  while ((err > tol) && (iter < max_iter)) {
    grd <- cptgrd_d2(beta,lambda,q, x, y, a, prob, tao, delta, G_c, G_c_tao)
    beta <- beta - stp_sz * grd
    iter <- iter + 1
    err <- sum(abs(grd))
    # print(iter)
  }
  # return(beta) # argmax
  if (result == "value") {
    return(dual_2(beta,lambda,q, x, y, a, prob, tao, delta, G_c, G_c_tao))
  } else {
    return(beta)
  }
}
ipw_U_d2 <- function(q_cons, n = 500, gamma_p, gamma_t, lambda_c,sim_num = 100){
  EU_hat<-0;ET_hat<-0;EU<-0;ET<-0;
  beta <- c();value <- c();lambda<-c()
  for (i in 1:sim_num) {
    dat <- sim_dat(n, gamma_p=gamma_p, gamma_t=gamma_t, lambda_c=lambda_c)
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
    
    lambda_opt <- optimize(solve_d2, interval = c(0, 100), q = q_cons,
                           x = x, y = T_obs, a = A, prob = prob, tao = tao, delta = delta, G_c=G_c , G_c_tao=G_c_tao,
                           maximum = TRUE, tol = 0.1)$maximum
    # print(lambda_opt)
    beta_hat <- solve_d2(lambda_opt, q = q_cons, x = x, y = T_obs, a = A, prob = prob, tao = tao, delta = delta, G_c=G_c ,G_c_tao=G_c_tao,
                         result = "beta")
    # 不再标准化系数
    
    # 以下均为对偶最优解，beta_hat的value_hat & value
    EU_hat <- ipw_U_stoc(beta_hat, x, T_obs, A, prob,tao,delta,G_c_tao)
    ET_hat <- ipw_T_stoc(beta_hat, x, T_obs, A, prob,tao,delta,G_c, G_c_tao)
    temp_eval <- eval_performance_g_stoc(beta_hat, gamma_p. = gamma_p,gamma_t. = gamma_t)
    EU <- temp_eval[1]
    ET <- temp_eval[2]
    
    beta <- rbind(beta, beta_hat)
    value <- rbind(value, c(EU_hat,EU,ET_hat,ET))
    lambda <- rbind(lambda, lambda_opt)
    colnames(value) <- c("EU_hat", "EU", "ET_hat","ET")
    
    detach(dat_obs)
    print(i)
  }
  list(beta, value, lambda)
}

ipw_U_d2(q_cons=2.6, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 1)

registerDoRNG(2023);a0 = proc.time()[3]
result_ipwU_d2_20 = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  ipw_U_d2(q_cons=2.0, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 25)};print(proc.time()[3] - a0)

registerDoRNG(2023);a0 = proc.time()[3]
result_ipwU_d2_23 = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  ipw_U_d2(q_cons=2.3, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 25)};print(proc.time()[3] - a0)

registerDoRNG(2023);a0 = proc.time()[3]
result_ipwU_d2_26 = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  ipw_U_d2(q_cons=2.6, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 25)};print(proc.time()[3] - a0)

registerDoRNG(2023);a0 = proc.time()[3]
result_ipwU_d2_29 = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  ipw_U_d2(q_cons=2.9, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 25)};print(proc.time()[3] - a0)

registerDoRNG(2023);a0 = proc.time()[3]
result_ipwU_d2_32 = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  ipw_U_d2(q_cons=3.2, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 25)};print(proc.time()[3] - a0)

result = break_result(result_ipwU_d2_32)
cat(
  "eta", apply(result[[1]], 2, mean),"\n",
  "eta_sd", apply(result[[1]], 2, sd),"\n",
  "V", apply(result[[2]], 2, mean),"\n",
  "V_sd", apply(result[[2]], 2, sd),"\n",
  "AA:", mean(result[[3]]), "\n",
  "AA_sd", sd(result[[3]]))
saveRDS(result, "D:\\论文\\Cure Model + ITR\\模拟\\复现整理版 - 副本\\U_ipw_d2_29.RDS")
########################## IPW: algorithm 2-penalty ################################
cptgrd_p2 <- function(beta,lambda,q, x, y, a, prob, tao, delta, G_c, G_c_tao){
  grd_u <- cptgrd_U(beta, x, y, a, prob, tao, delta, G_c_tao)
  # 罚函数分段求导
  if (ipw_T_stoc(beta, x, T_obs, A, prob,tao,delta,G_c, G_c_tao)>=q){
    grd_t = 0
  } else{
    grd_t <- cptgrd_T(beta, x, y, a, prob, tao, delta, G_c, G_c_tao)
  }
  grd = -grd_u - grd_t*lambda # 正负号（原始代码是负的。感觉是错的）
  return(grd)
}
solve_p2 <- function(lambda, q, x, y, a, prob, tao, delta, G_c, G_c_tao,
                     stp_sz = 0.01, tol = 1e-3, max_iter = 10000, result="beta") {
  beta <- rep(0, length(x[1,])) # 初始点设为原点
  iter <- 0 # 迭代次数
  err <- 999 # 迭代误差
  # stp_sz: 步长
  # tol:误差容忍度
  while ((err > tol) && (iter < max_iter)) {
    grd <- cptgrd_p2(beta,lambda,q, x, y, a, prob, tao, delta, G_c, G_c_tao)
    beta <- beta - stp_sz * grd
    iter <- iter + 1
    err <- sum(abs(grd))
    # print(iter)
  }
  # return(beta) # argmax
  if (result == "value") {
    return(dual_1(beta,lambda,q, x, y, a, prob, tao, delta, G_c, G_c_tao))
  } else {
    return(beta)
  }
}
MaxU_ipw_p2 <- function(q_cons, n = 500, gamma_p, gamma_t, lambda_c,sim_num = 100){
  EU_hat<-0;ET_hat<-0;EU<-0;ET<-0;
  beta <- c();value <- c();lambda<-c()
  for (i in 1:sim_num) {
    dat <- sim_dat(n, gamma_p=gamma_p, gamma_t=gamma_t, lambda_c=lambda_c)
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
    
    # 罚函数
    lambda0 = 2
    c=1.5
    # lambda_seq = 0.25*lambda0*c^(0:13)
    lambda_seq = 0.25*lambda0*c^(0:13)
    i=0; EU_hat=0
    while((ET_hat<q_cons) & (i<=13)){
      i=i+1
      lambda_temp = lambda_seq[i]
      beta_hat  = solve_p2(lambda_temp, q = q_cons, x = x, y = T_obs, a = A, prob = prob, tao = tao, delta = delta, 
                           G_c=G_c, G_c_tao = G_c_tao,result = "beta")
      ET_hat <- ipw_T_stoc(beta_hat, x, T_obs, A, prob,tao,delta,G_c, G_c_tao)
    }
    # 不再标准化系数
    
    # 以下均为对偶最优解，beta_hat的value_hat & value
    EU_hat <- ipw_U_stoc(beta_hat, x, T_obs, A, prob,tao,delta,G_c_tao)
    ET_hat <- ipw_T_stoc(beta_hat, x, T_obs, A, prob,tao,delta,G_c, G_c_tao)
    temp_eval <- eval_performance_g_stoc(beta_hat, gamma_p. = gamma_p,gamma_t. = gamma_t)
    EU <- temp_eval[1]
    ET <- temp_eval[2]
    
    beta <- rbind(beta, beta_hat)
    value <- rbind(value, c(EU_hat,EU,ET_hat,ET))
    lambda <- rbind(lambda, lambda_temp)
    colnames(value) <- c("EU_hat", "EU", "ET_hat","ET")
    
    detach(dat_obs)
    print(i)
  }
  list(beta, value, lambda)
}
MaxU_ipw_p2(q_cons=0.001, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 1)

registerDoRNG(2023);a0 = proc.time()[3]
result_U_ipw_p1_20 = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  MaxU_ipw_p2(q_cons=2.0, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 25)};print(proc.time()[3] - a0)
registerDoRNG(2023);a0 = proc.time()[3]
result_U_ipw_p1_23 = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  MaxU_ipw_p2(q_cons=2.3, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 25)};print(proc.time()[3] - a0)
registerDoRNG(2023);a0 = proc.time()[3]
result_U_ipw_p1_26 = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  MaxU_ipw_p2(q_cons=2.6, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 25)};print(proc.time()[3] - a0)
registerDoRNG(2023);a0 = proc.time()[3]
result_U_ipw_p1_29 = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  MaxU_ipw_p2(q_cons=2.9, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 25)};print(proc.time()[3] - a0)
registerDoRNG(2023);a0 = proc.time()[3]
result_U_ipw_p1_32 = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  MaxU_ipw_p2(q_cons=3.2, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 25)};print(proc.time()[3] - a0)

result = break_result(result_U_ipw_p1_32)
cat(
  "beta", apply(result[[1]], 2, mean),"\n",
  "beta_sd", apply(result[[1]], 2, sd),"\n",
  "V", apply(result[[2]], 2, mean),"\n",
  "V_sd", apply(result[[2]], 2, sd),"\n",
  "lambda:", mean(result[[3]]))
saveRDS(result, "D:\\论文\\Cure Model + ITR\\模拟\\复现整理版 - 副本\\U_ipw_p1_32.RDS")
########################## AIPW: algorithm 1-dual ################################
# dual function
dual_aipw_d2 <- function(beta,lambda,q, x, y, a, prob, tao, delta, G_c, G_c_tao , t.es_0, t.es_1, p.es_0, p.es_1) {
  val <- -aipw_U_stoc(beta, x, y, a, prob, tao, delta, G_c_tao, p.es_0, p.es_1) +
    lambda * (q - aipw_T_stoc(beta, x, y, a, prob, tao, delta, G_c, G_c_tao , t.es_0, t.es_1, p.es_0, p.es_1))
  return(val)
}
# dual function gradient
cptgrd_aipw_d2 <- function(beta,lambda,q, x, y, a, prob, tao, delta, G_c, G_c_tao , t.es_0, t.es_1, p.es_0, p.es_1) {
  grd_u <- cptgrd_U_aipw(beta, x, y, a, prob, tao, delta, G_c_tao, p.es_0, p.es_1)
  grd_t <- cptgrd_T_aipw(beta, x, y, a, prob, tao, delta, G_c, G_c_tao , t.es_0, t.es_1, p.es_0, p.es_1)
  grd <- -grd_u - lambda * grd_t
  return(grd)
}

# Given \lambda, 梯度下降法求解最优 \beta
# 可以优化？
solve_aipw_d2 <- function(lambda, q, x, y, a, prob, tao, delta, G_c, G_c_tao , t.es_0, t.es_1, p.es_0, p.es_1,
                     stp_sz = 0.01, tol = 1e-3, max_iter = 10000, result="value") {
  beta <- rep(0, length(x[1,])) # 初始点设为原点
  iter <- 0 # 迭代次数
  err <- 999 # 迭代误差
  # stp_sz: 步长
  # tol:误差容忍度
  while ((err > tol) && (iter < max_iter)) {
    grd <- cptgrd_aipw_d2(beta,lambda,q, x, y, a, prob, tao, delta, G_c, G_c_tao , t.es_0, t.es_1, p.es_0, p.es_1)
    beta <- beta - stp_sz * grd
    iter <- iter + 1
    err <- sum(abs(grd))
    # print(iter)
  }
  # return(beta) # argmax
  if (result == "value") {
    return(dual_aipw_d2(beta,lambda,q, x, y, a, prob, tao, delta, G_c, G_c_tao , t.es_0, t.es_1, p.es_0, p.es_1))
  } else {
    return(beta)
  }
}
aipw_U_d2 <- function(q_cons, n = 500, gamma_p, gamma_t, lambda_c, sim_num = 100){
  EU_hat<-0;ET_hat<-0;EU<-0;ET<-0;
  beta <- c();value <- c();lambda<-c()
  for (i in 1:sim_num) {
    dat <- sim_dat(n, gamma_p=gamma_p, gamma_t=gamma_t, lambda_c=lambda_c)
    dat_obs <- dat[[1]]
    dat_true <- dat[[2]]
    
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
    
    fit <- gfcure(Surv(T_obs, delta) ~ x1+x2+x1*A+x2*A, ~ x1+x2+x1*A+x2*A, data= dat_obs, dist = "lognormal")
    gamma_p_hat <- -fit$coef[8:13]
    logit_p <- x%*%gamma_p_hat[1:3]; p.es_0 <- exp(logit_p)/(1+exp(logit_p)) 
    logit_p <- x%*%gamma_p_hat[1:3] + x%*%gamma_p_hat[4:6]; p.es_1 <- exp(logit_p)/(1+exp(logit_p))
    gamma_t_hat <- fit$coef[2:7]
    t.es_0 <- exp(x%*%gamma_t_hat[1:3]) * exp(exp(fit$coef[1])^2/2)
    t.es_1 <- exp(x%*%gamma_t_hat[1:3] + 1* x %*%gamma_t_hat[4:6]) * exp(exp(fit$coef[1])^2/2)
    
    lambda_opt <- optimize(solve_aipw_d2, interval = c(0, 100), q = q_cons,
                           x = x, y = T_obs, a = A, prob = prob, tao = tao, delta = delta, G_c=G_c , G_c_tao=G_c_tao,
                           t.es_0=t.es_0, t.es_1=t.es_1, p.es_0=p.es_0, p.es_1=p.es_1,maximum = TRUE, tol = 0.1)$maximum
    # print(lambda_opt)
    beta_hat <- solve_aipw_d2(lambda_opt, q = q_cons, x = x, y = T_obs, a = A, prob = prob, tao = tao, delta = delta, G_c=G_c, G_c_tao=G_c_tao,
                              t.es_0=t.es_0, t.es_1=t.es_1, p.es_0=p.es_0, p.es_1=p.es_1,result = "beta")
    # 不再标准化系数
    
    # 以下均为对偶最优解，beta_hat的value_hat & value
    EU_hat <- aipw_U_stoc(beta_hat, x, T_obs, A, prob,tao,delta,G_c_tao, p.es_0, p.es_1)
    ET_hat <- aipw_T_stoc(beta_hat, x, T_obs, A, prob,tao,delta,G_c, G_c_tao, t.es_0, t.es_1, p.es_0, p.es_1)
    temp_eval <- eval_performance_g_stoc(beta_hat, gamma_p. = gamma_p,gamma_t. = gamma_t)
    EU <- temp_eval[1]
    ET <- temp_eval[2]
    
    beta <- rbind(beta, beta_hat)
    value <- rbind(value, c(EU_hat,EU,ET_hat,ET))
    lambda <- rbind(lambda, lambda_opt)
    colnames(value) <- c("EU_hat", "EU", "ET_hat","ET")
    
    detach(dat_obs)
    print(i)
  }
  list(beta, value, lambda)
}
aipw_U_d2(q_cons=2.6, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 1)

registerDoRNG(2023);a0 = proc.time()[3]
result_aipwU_d2_20 = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  aipw_U_d2(q_cons=2.0, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 25)};print(proc.time()[3] - a0)

registerDoRNG(2023);a0 = proc.time()[3]
result_aipwU_d2_23 = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  aipw_U_d2(q_cons=2.3, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 25)};print(proc.time()[3] - a0)

registerDoRNG(2023);a0 = proc.time()[3]
result_aipwU_d2_26 = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  aipw_U_d2(q_cons=2.6, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 25)};print(proc.time()[3] - a0)

registerDoRNG(2023);a0 = proc.time()[3]
result_aipwU_d2_29 = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  aipw_U_d2(q_cons=2.9, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 25)};print(proc.time()[3] - a0)

registerDoRNG(2023);a0 = proc.time()[3]
result_aipwU_d2_32 = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  aipw_U_d2(q_cons=3.2, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 25)};print(proc.time()[3] - a0)

result = break_result(result_aipwU_d2_32)
cat(
  "beta", apply(result[[1]], 2, mean),"\n",
  "beta_sd", apply(result[[1]], 2, sd),"\n",
  "V", apply(result[[2]], 2, mean),"\n",
  "V_sd", apply(result[[2]], 2, sd),"\n",
  "lambda:", mean(result[[3]]))
saveRDS(result, "D:\\论文\\Cure Model + ITR\\模拟\\复现整理版 - 副本\\U_aipw_d2_29.RDS")

########################## AIPW: algorithm 2-penalty ################################
cptgrd_p2 <- function(beta,lambda,q, x, y, a, prob, tao, delta, G_c, G_c_tao, t.es_0, t.es_1, p.es_0, p.es_1){
  grd_u <- cptgrd_U_aipw(beta, x, y, a, prob, tao, delta, G_c_tao, p.es_1, p.es_0)
  # 罚函数分段求导
  if (aipw_T_stoc(beta, x, y, a, prob, tao, delta,G_c, G_c_tao,t.es_1,t.es_0,p.es_1,p.es_0)>=q){
    grd_t = 0
  } else{
    grd_t <- cptgrd_T_aipw(beta, x, y, a, prob, tao, delta, G_c, G_c_tao, t.es_0, t.es_1, p.es_0, p.es_1)
  }
  grd = -grd_u - grd_t*lambda # 正负号（原始代码是负的。感觉是错的）
  return(grd)
}
solve_p2 <- function(lambda, q, x, y, a, prob, tao, delta, G_c, G_c_tao, t.es_0, t.es_1, p.es_0, p.es_1,
                     stp_sz = 0.01, tol = 1e-3, max_iter = 10000, result="beta") {
  beta <- rep(0, length(x[1,])) # 初始点设为原点
  iter <- 0
  err <- 999
  while ((err > tol) && (iter < max_iter)) {
    grd <- cptgrd_p2(beta,lambda,q, x, y, a, prob, tao, delta, G_c, G_c_tao, t.es_0, t.es_1, p.es_0, p.es_1)
    beta <- beta - stp_sz * grd
    iter <- iter + 1
    err <- sum(abs(grd))
    # print(iter)
  }
  # return(beta) # argmax
  if (result == "value") {
    return(dual_1(beta,lambda,q, x, y, a, prob, tao, delta, G_c, G_c_tao, t.es_0, t.es_1, p.es_0, p.es_1))
  } else {
    return(beta)
  }
}
MaxU_aipw_p2 <- function(q_cons, n = 500, gamma_p, gamma_t, lambda_c,sim_num = 100){
  EU_hat<-0;ET_hat<-0;EU<-0;ET<-0;
  beta <- c();value <- c();lambda<-c()
  for (i in 1:sim_num) {
    dat <- sim_dat(n, gamma_p=gamma_p, gamma_t=gamma_t, lambda_c=lambda_c)
    dat_obs <- dat[[1]]
    dat_true <- dat[[2]]
    
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
    
    fit <- gfcure(Surv(T_obs, delta) ~ x1+x2+x1*A+x2*A, ~ x1+x2+x1*A+x2*A, data= dat_obs, dist = "lognormal")
    gamma_p_hat <- -fit$coef[8:13]
    logit_p <- x%*%gamma_p_hat[1:3]; p.es_0 <- exp(logit_p)/(1+exp(logit_p)) 
    logit_p <- x%*%gamma_p_hat[1:3] + x%*%gamma_p_hat[4:6]; p.es_1 <- exp(logit_p)/(1+exp(logit_p))
    gamma_t_hat <- fit$coef[2:7]
    t.es_0 <- exp(x%*%gamma_t_hat[1:3]) * exp(exp(fit$coef[1])^2/2)
    t.es_1 <- exp(x%*%gamma_t_hat[1:3] + 1* x %*%gamma_t_hat[4:6]) * exp(exp(fit$coef[1])^2/2)
    
    # 罚函数
    lambda0 = 2
    c=1.5
    lambda_seq = 0.25*lambda0*c^(0:14)
    i=0; EU_hat=0
    while((ET_hat<q_cons) & (i<=14)){
      i=i+1
      lambda_temp = lambda_seq[i]
      beta_hat  = solve_p2(lambda_temp, q = q_cons, x = x, y = T_obs, a = A, prob = prob, tao = tao, delta = delta, 
                           G_c=G_c, G_c_tao=G_c_tao , t.es_0 = t.es_0, t.es_1 = t.es_1, p.es_0 = p.es_0, p.es_1 = p.es_1, 
                           result = "beta")
      ET_hat <- aipw_T_stoc(beta_hat, x, T_obs, A, prob,tao,delta, G_c, G_c_tao, t.es_0, t.es_1 ,p.es_0, p.es_1)
    }
    # 不再标准化系数
    
    # 以下均为对偶最优解，beta_hat的value_hat & value
    EU_hat <- aipw_U_stoc(beta_hat, x, T_obs, A, prob,tao, delta, G_c_tao ,p.es_0,p.es_1)
    ET_hat <- aipw_T_stoc(beta_hat, x, T_obs, A, prob,tao,delta , G_c, G_c_tao, t.es_0, t.es_1, p.es_0, p.es_1)
    temp_eval <- eval_performance_g_stoc(beta_hat, gamma_p. = gamma_p,gamma_t. = gamma_t)
    EU <- temp_eval[1]
    ET <- temp_eval[2]
    
    beta <- rbind(beta, beta_hat)
    value <- rbind(value, c(EU_hat,EU,ET_hat,ET))
    lambda <- rbind(lambda, lambda_temp)
    colnames(value) <- c("EU_hat", "EU", "ET_hat","ET")
    
    detach(dat_obs)
    print(i)
  }
  list(beta, value, lambda)
}
MaxU_aipw_p2(q_cons=3.0, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 1)

registerDoRNG(2023);a0 = proc.time()[3]
result_U_aipw_p2_20 = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  MaxU_aipw_p2(q_cons=2.0, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 25)};print(proc.time()[3] - a0)
registerDoRNG(2023);a0 = proc.time()[3]
result_U_aipw_p2_23 = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  MaxU_aipw_p2(q_cons=2.3, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 25)};print(proc.time()[3] - a0)
registerDoRNG(2023);a0 = proc.time()[3]
result_U_aipw_p2_26 = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  MaxU_aipw_p2(q_cons=2.6, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 25)};print(proc.time()[3] - a0)
registerDoRNG(2023);a0 = proc.time()[3]
result_U_aipw_p2_29 = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  MaxU_aipw_p2(q_cons=2.9, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 25)};print(proc.time()[3] - a0)
registerDoRNG(2023);a0 = proc.time()[3]
result_U_aipw_p2_32 = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  MaxU_aipw_p2(q_cons=3.2, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 25)};print(proc.time()[3] - a0)

result = break_result(result_U_aipw_p2_32)
cat(
  "beta", apply(result[[1]], 2, mean),"\n",
  "beta_sd", apply(result[[1]], 2, sd),"\n",
  "V", apply(result[[2]], 2, mean),"\n",
  "V_sd", apply(result[[2]], 2, sd),"\n",
  "lambda:", mean(result[[3]]))
saveRDS(result, "D:\\论文\\Cure Model + ITR\\模拟\\复现整理版 - 副本\\U_aipw_p2_32.RDS")
