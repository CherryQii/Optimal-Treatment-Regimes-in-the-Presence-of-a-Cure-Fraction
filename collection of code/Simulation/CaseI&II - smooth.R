# ########################## package loading ################################
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
  # g <- as.numeric(x %*% eta > 0)
  c0 <- 4^(1/3)
  etaX <- x %*% eta
  h <- c0*length(a)^(- 1/3)*sd(etaX)
  g <- pnorm(etaX/h, mean = 0, sd = 1)
  
  wts <- a* g * 1 / prob + (1-a) * (1 - g) * (1 / (1 - prob)) # 逆概率权重
  wts <- wts / G_c_tao # 删失逆概率
  y <- as.numeric(y > tao) # tao标记
  val <- mean(wts * y * (1 - delta))
  val
}
ipw_T <- function(eta, x, y, a, prob, tao, delta, G_c, G_c_tao) {
  # g <- as.numeric(x %*% eta > 0)
  c0 <- 4^(1/3)
  etaX <- x %*% eta
  h <- c0*length(a)^(- 1/3)*sd(etaX)
  g <- pnorm(etaX/h, mean = 0, sd = 1)
  wts <- a* g * 1 / prob + (1-a) * (1 - g) * (1 / (1 - prob)) # 逆概率权重
  v_t <- mean(wts * y * delta / (G_c))
  v_u <- ipw_U(eta, x, y, a, prob, tao, delta, G_c_tao)
  val <- v_t/(1 - v_u)
  val
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

# 1.ipw-U-true
registerDoRNG(2023);a0 = proc.time()[3]
result_ipwU_t_s = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  sim_U_ipw(n=1000, gamma_p, gamma_t, lambda_c,sim_num = 10, mis_ps = FALSE, dist="norm")};print(proc.time()[3] - a0)

`#  1.ipw-U-false
registerDoRNG(2023);a0 = proc.time()[3]
result_ipwU_f_s = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
  sim_U_ipw(n=1000, gamma_p, gamma_t, lambda_c,sim_num = 10, mis_ps = TRUE, dist="norm")};print(proc.time()[3] - a0)

result = break_result(result_ipwU_t_s)
cat(
  "eta", apply(result[[1]], 2, mean),"\n",
  "eta_sd", apply(result[[1]], 2, sd),"\n",
  "V", apply(result[[2]], 2, mean),"\n",
  "V_sd", apply(result[[2]], 2, sd),"\n",
  "AA:", mean(result[[3]]), "\n",
  "AA_sd", sd(result[[3]]))
