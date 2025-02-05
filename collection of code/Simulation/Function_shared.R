break_result <- function(result){
  Val_result <- c()
  eta_result <- c()
  for (i in 1:10){
    eta_result <- as.data.frame(rbind(eta_result, result[i,1][[1]]))
    Val_result <- as.data.frame(rbind(Val_result, result[i,2][[1]]))
  }
  AA <- unlist(result[,3])
  
  return(list(eta_result, Val_result, AA))
}
sim_dat <- function(n=1000,gamma_pi=c(-1,0.8,0.8),gamma_p, gamma_t,lambda_c, dist = 'norm'){
  # Params:
  # n:sample size
  # gamma_pi: params of treatment propensity_score
  # gamma_p: params of cure rate(main effect + treatment effect)
  # gamma_t: params of uncured survival time(main effect + treatment effect)
  # lambda_c: rexp(n,lambda_c), censoring rate
  x1 <- runif(n,-1.5,1.5)
  x2 <- runif(n,-1,1)
  
  # Propensity score: P(A=1|X)
  x <- cbind(1,x1^2,x2^2)
  logit_pi <- x%*%gamma_pi
  .pi <- exp(logit_pi)/(1+exp(logit_pi)) 
  A <- rbinom(n,1,.pi)
  
  # cure rate: p(x,a)
  x <- cbind(1,x1,x2)
  logit_p <- x%*%gamma_p[1:3] + A*x%*%gamma_p[4:6]
  p <- exp(logit_p)/(1+exp(logit_p)) 
  cure_status_true <- rbinom(n,1,p) # 1为cured，0为uncured
  
  dat <- as.data.frame(cbind(x1,x2,A))
  
  # generate survival time
  # uncured portion: follow varying coefficient exponential distribution. T_obs = min(T,C)
  # cured portion: T=inf. So, T_obs = C
  z=x
  # lambda <- 1 * exp(z[cure_status_true==0,]%*%gamma_t[1:3] + A[cure_status_true==0]*z[cure_status_true==0,]%*%gamma_t[4:6])
  if (dist == 'Weibull'){ # extreme value distribution
    epsilon <- -0.1*log(-log(runif(n-sum(cure_status_true))))}
  if (dist == 'norm'){  
  epsilon <- rnorm(n-sum(cure_status_true), 0 ,1)}
  
  T <- exp(z[cure_status_true==0,]%*%gamma_t[1:3] + A[cure_status_true==0]*z[cure_status_true==0,]%*%gamma_t[4:6]+epsilon)
  # tau=50;T <- ifelse(T>tau, tau, T)
  C <- rexp(n,lambda_c) # 与协变量、A无关，随机删失
  # plot(density(T));lines(density(C))
  
  dat[cure_status_true==1,'T_obs'] <- C[cure_status_true==1]
  dat[cure_status_true==0,'T_obs'] <- ifelse(T<C[cure_status_true==0],T,C[cure_status_true==0])
  dat[cure_status_true==1,'delta'] <- 0
  dat[cure_status_true==0,'delta'] <- ifelse(T<C[cure_status_true==0],1,0)
  
  dat_true <- as.data.frame(cbind(.pi,p,cure_status_true))
  colnames(dat_true) <- c('.pi','p','cure_status_true')
  
  # dat:(x1,x2,A,T_obs,delta)
  # dat_true:(.pi,p,cure_status_true, p.es_0,p.es_1,t.es_0,t.es_1)
  # uncured survival time
  list(dat,dat_true,T)
}
########################## true value given treatment(eta) #######################
eval_performance_g <- function(eta, size = 1e5, gamma_p., gamma_t., RCT=FALSE, dist = 'norm') {
  x1 <- runif(size, -1.5, 1.5)
  x2 <- runif(size, -1, 1)
  x <- cbind(1, x1, x2)
  A <- as.numeric(x %*% eta > 0)
  if (RCT==TRUE){
    A <- rbinom(length(x1),1, 0.5)
  }
  
  A_U_true <- as.numeric(x %*% gamma_p.[4:6] > 0)
  AA_U <- mean(A_U_true == A)
  A_T_true <- as.numeric(x %*% gamma_t.[4:6] > 0) # 只有特定情况下才成立
  AA_T <- mean(A_T_true == A)
  
  logit_p <- x %*% gamma_p[1:3] + A * x %*% gamma_p[4:6]
  p <- exp(logit_p) / (1 + exp(logit_p))
  cure_status_true <- rbinom(size, 1, p)
  EU = mean(cure_status_true)
  
  z <- x
  # lambda <- 1 * exp((z[cure_status_true == 0, ] %*% gamma_t[1:3] + A[cure_status_true == 0] * z[cure_status_true == 0, ] %*% gamma_t[4:6]))
  # T <- rexp(size - sum(cure_status_true), lambda)
  if (dist == 'Weibull') {
    epsilon <- -0.1*log(-log(runif(size-sum(cure_status_true))))}
  if (dist == 'norm') {
    epsilon <- rnorm(size - sum(cure_status_true), 0, 1)}
  
  T <- exp(z[cure_status_true==0,]%*%gamma_t[1:3] + A[cure_status_true==0]*z[cure_status_true==0,]%*%gamma_t[4:6]+epsilon)
  # T[T > 50] <- 50
  ET = mean(T)
  
  result <- c(EU, ET, AA_U, AA_T)
  return(result)
}
EDA <- function(gamma_p, gamma_t, lambda_c, RCT = FALSE, dist = 'norm'){
  # input the dat_generating data
  dat <- sim_dat(n=1e5,gamma_p=gamma_p,lambda_c=lambda_c, gamma_t = gamma_t, dist = dist)
  dat_obs <- dat[[1]]
  dat_true <- dat[[2]]
  
  censoring_rate <- mean(dat_obs$delta==0)
  EU_sample <- mean(dat_true$cure_status_true)
  ET_sample <- mean(dat[[3]])
  
  result <- eval_performance_g(RCT=TRUE,gamma_p[4:6], gamma_p. = gamma_p,gamma_t. = gamma_t, size=1e6, dist = dist)
  cure_rct <- result[1]
  T_rct <- result[2]
  
  result <- eval_performance_g(gamma_p[4:6], gamma_p. = gamma_p,gamma_t. = gamma_t, size=1e6, dist = dist)
  cure_opt1 <- result[1]
  T_opt1 <- result[2]
  
  result <- eval_performance_g(gamma_t[4:6], gamma_p. = gamma_p,gamma_t. = gamma_t, size = 1e6, dist = dist)
  cure_opt2 <- result[1]
  T_opt2 <- result[2]
  # list with name
  return(list(censoring_rate = censoring_rate, EU_sample = EU_sample, ET_sample = ET_sample, cure_rct = cure_rct, T_rct = T_rct, cure_opt1 = cure_opt1, T_opt1 = T_opt1, cure_opt2 = cure_opt2, T_opt2 = T_opt2))
}
