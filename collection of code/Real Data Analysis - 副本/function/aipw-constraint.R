aipw_U_stoc <- function(beta, x, y, a, prob, tao, delta, Gc_tao, p.es_0, p.es_1) {
    g <- 1 / (1 + exp(x %*% beta))
    wts <- a * g * 1 / prob + (1 - a) * (1 - g) * (1 / (1 - prob)) # 逆概率权重
    y <- as.numeric(y > tao) # tao标记
    val_tao <- wts * y * (1 - delta)
    val_tao <- val_tao / Gc_tao

    p.es <- g * p.es_1 + (1 - g) * p.es_0
    val_mi <- (wts - 1) * p.es # 标准aipw

    val <- mean(val_tao - val_mi)
    return(val)
}
aipw_T_stoc <- function(beta, x, y, a, prob, tao, delta, Gc_tao, Gc_y, t.es_0, t.es_1, p.es_0, p.es_1) {
  g <- 1 / (1 + exp(x %*% beta))
  wts <- a * g * 1 / prob + (1 - a) * (1 - g) * (1 / (1 - prob))
  H <- g * t.es_1 * p.es_0 + (1 - g) * t.es_0 * p.es_1
  v_t <- mean(wts * y * delta / Gc_y + (1-wts) * H)
  
  v_u <- aipw_U_stoc(beta, x, y, a, prob, tao, delta, Gc_tao, p.es_0, p.es_1)
  val_t <- v_t / (1 - v_u)
  return(val_t)
}
cptgrd_U <- function(beta, x, y, a, prob, tao, delta, Gc_tao, p.es_0, p.es_1) {
  d <- length(beta)
  m0 <- aipw_U_stoc(beta, x, y, a, prob, tao, delta, Gc_tao,p.es_0, p.es_1)
  grd <- rep(0, d)
  for (j in 1:d) {
    tmp <- rep(0, d)
    tmp[j] <- 0.05 * (2 * rbinom(1, size = 1, prob = 0.5) - 1)
    grd[j] <- (aipw_U_stoc(beta+tmp, x, y, a, prob, tao, delta, Gc_tao,p.es_0, p.es_1) - m0) / tmp[j] # ???ţ?
  }
  return(grd)
}
cptgrd_T <- function(beta, x, y, a, prob, tao, delta, Gc_tao, Gc_y, t.es_0, t.es_1, p.es_0, p.es_1) {
  d <- length(beta)
  m0 <- aipw_T_stoc(beta, x, y, a, prob, tao, delta, Gc_tao, Gc_y, t.es_0, t.es_1, p.es_0, p.es_1)
  grd <- rep(0, d)
  for (j in 1:d) {
    tmp <- rep(0, d)
    tmp[j] <- 0.01 * (2 * rbinom(1, size = 1, prob = 0.5) - 1)
    grd[j] <- (aipw_T_stoc(beta+tmp, x, y, a, prob, tao, delta, Gc_tao, Gc_y, t.es_0, t.es_1, p.es_0, p.es_1) - m0) / tmp[j] # ???ţ?
  }
  return(grd)
}

####################################### 无约束 ######################################
# Adam
solve_U_aipw <- function(x, y, a, prob, tao, delta, Gc_tao, p.es_0, p.es_1,
                    stp_sz = 0.01, tol = 1e-3, max_iter = 10000, beta1 = 0.9, beta2 = 0.999,
                    num_starts=5,range=c(-10,10)
                    ) {
  best_beta <- NULL;best_obj_value <- -Inf
  for (i in 1:num_starts) {
      if(i==1){
        beta <- rep(0, length(x[1,]))
      } else{
        beta <- runif(length(x[1,]),range[1],range[2])
      }
      iter <- 0;err <- 999
      lr <- stp_sz
      v <- rep(0, length(beta));v_hat <- rep(0, length(beta))
      m <- rep(0, length(beta));m_hat <- rep(0, length(beta))
      while ((err > tol) && (iter < max_iter)) {
        grd <- cptgrd_U(beta, x, y, a, prob, tao, delta, Gc_tao, p.es_0, p.es_1)
        if (sum(grd)==0){
          norm_grd <- rep(0, length(beta))
        }else{norm_grd <- grd / sqrt(sum(grd^2))}
      
        for (j in 1:length(grd)) {
            v[j] <- beta2 * v[j] + (1 - beta2) * norm_grd[j]^2
            m[j] <- beta1 * m[j] + (1 - beta1) * norm_grd[j]
            v_hat[j] <- v[j] / (1 - beta2^(iter + 1))
            m_hat[j] <- m[j] / (1 - beta1^(iter + 1))
            beta[j] <- beta[j] + lr * m_hat[j] / (sqrt(v_hat[j]) + 1e-8)
        }
        beta_new <- beta

        beta <- beta_new
        iter <- iter + 1
        err <- sum(abs(grd))
      }
      obj_value <- aipw_U_stoc(beta, x, y, a, prob, tao, delta, Gc_tao, p.es_0, p.es_1)
      # browser()
      if (obj_value > best_obj_value) {
          best_obj_value <- obj_value
          best_beta <- beta
      }
      print(iter);print(err)
  }
  
  return(best_beta)
}
solve_T_aipw <- function(x, y, a, prob, tao, delta, Gc_tao, Gc_y, t.es_0, t.es_1, p.es_0, p.es_1,
                    stp_sz = 0.001, tol = 1e-3, max_iter = 20000,beta1 = 0.9, beta2 = 0.999,
                    num_starts=5,range=c(-10,10)
                    ) {
  best_beta <- NULL;best_obj_value <- -Inf
  for (i in 1:num_starts) {
      if(i==1){
        beta <- rep(0, length(x[1,]))
      } else{
        beta <- runif(length(x[1,]),range[1],range[2])
      }
      iter <- 0;err <- 999
      lr <- stp_sz
      v <- rep(0, length(beta));v_hat <- rep(0, length(beta))
      m <- rep(0, length(beta));m_hat <- rep(0, length(beta))
      while ((err > tol) && (iter < max_iter)) {
        grd <- cptgrd_T(beta, x, y, a, prob, tao, delta, Gc_tao, Gc_y, t.es_0, t.es_1, p.es_0, p.es_1)
        if (sum(grd)==0){
          norm_grd <- rep(0, length(beta))
        }else{norm_grd <- grd / sqrt(sum(grd^2))}
        for (j in 1:length(grd)) {
            v[j] <- beta2 * v[j] + (1 - beta2) * norm_grd[j]^2
            m[j] <- beta1 * m[j] + (1 - beta1) * norm_grd[j]
            v_hat[j] <- v[j] / (1 - beta2^(iter + 1))
            m_hat[j] <- m[j] / (1 - beta1^(iter + 1))
            beta[j] <- beta[j] + lr * m_hat[j] / (sqrt(v_hat[j]) + 1e-8)
        }
        beta_new <- beta

        beta <- beta_new
        iter <- iter + 1
        err <- sum(abs(grd))
      }
      obj_value <- aipw_T_stoc(beta, x, y, a, prob, tao, delta, Gc_tao, Gc_y, t.es_0, t.es_1, p.es_0, p.es_1)
      if (obj_value > best_obj_value) {
          best_obj_value <- obj_value
          best_beta <- beta
      }
  }
  print(iter);print(err)
  return(best_beta)
}
################################ max U, s.t. T ##############################
aipw_T_stoc_q <- function(q, beta, x, y, a, prob, tao, delta, Gc_tao, Gc_y, t.es_0, t.es_1, p.es_0, p.es_1) {
  g <- 1 / (1 + exp(x %*% beta))
  wts <- a * g * 1 / prob + (1 - a) * (1 - g) * (1 / (1 - prob)) # 逆概率权重
  t.es <- g * t.es_1*p.es_0 + (1 - g) * t.es_0 * p.es_1
  # .T <- ifelse(wts * y * delta / Gc_y > tao, tao, wts * y * delta / Gc_y)
  .T <- wts * y * delta / Gc_y
  v_t <- mean(.T + (1 - wts) * t.es)

  v_u <- aipw_U_stoc(beta, x, y, a, prob, tao, delta, Gc_tao, p.es_0, p.es_1)
  val_t <- v_t / (1 - v_u)

  val <- -(val_t - q)^2
  # val <- mean(val_t)
  return(val)
}
cptgrd_T_q <- function(q, beta, x, y, a, prob, tao, delta, Gc_tao, Gc_y, t.es_0, t.es_1, p.es_0, p.es_1) {
  d <- length(beta)
  m0 <- aipw_T_stoc_q(q, beta, x, y, a, prob, tao, delta, Gc_tao, Gc_y, t.es_0, t.es_1, p.es_0, p.es_1)
  grd <- rep(0, d)
  for (j in 1:d) {
    tmp <- rep(0, d)
    tmp[j] <- 0.01 * (2 * rbinom(1, size = 1, prob = 0.5) - 1)
    grd[j] <- (aipw_T_stoc_q(q, beta+tmp, x, y, a, prob, tao, delta, Gc_tao, Gc_y, t.es_0, t.es_1, p.es_0, p.es_1) - m0) / tmp[j] # ???ţ?
  }
  return(grd)
}
cptgrd_p1 <- function(beta,lambda,q, x, y, a, prob, tao, delta, Gc_tao, Gc_y, t.es_0, t.es_1, p.es_0, p.es_1){
  grd_u <- cptgrd_U(beta, x, y, a, prob, tao, delta, Gc_tao, p.es_0, p.es_1)
  # 罚函数分段求导
  if (aipw_T_stoc(beta, x, y, a, prob, tao, delta , Gc_tao, Gc_y, t.es_0, t.es_1, p.es_0, p.es_1)>=q){
    grd_t = 0
  } else{
    # grd_t <- cptgrd_T(beta, x, y, a, prob, tao, delta, Gc_tao, Gc_y)
    grd_t <- cptgrd_T_q(q, beta, x, y, a, prob, tao, delta, Gc_tao, Gc_y, t.es_0, t.es_1, p.es_0, p.es_1)
  }
  grd = grd_u + grd_t*lambda # 正负号（原始代码是负的。感觉是错的）
  return(grd)
}
# Given \lambda, 梯度下降法求解最优 \beta
solve_p1 <- function(
  lambda, q, x, y, a, prob, tao, delta, Gc_tao, Gc_y, t.es_0, t.es_1, p.es_0, p.es_1,
                     stp_sz = 0.01, tol = 1e-3, max_iter = 20000,
                     beta1 = 0.9, beta2 = 0.999,num_starts=5,range=c(-10,10)) {
  best_beta <- NULL;best_obj_value <- -Inf
  for (i in 1:num_starts) {
      if(i==1){
        beta <- rep(0, length(x[1,]))
      } else{
        beta <- runif(length(x[1,]),range[1],range[2])
      }
      iter <- 0;err <- 999
      lr <- stp_sz
      v <- rep(0, length(beta));v_hat <- rep(0, length(beta))
      m <- rep(0, length(beta));m_hat <- rep(0, length(beta))
      while ((err > tol) && (iter < max_iter)) {
        grd <- cptgrd_p1(beta,lambda,q, x, y, a, prob, tao, delta, Gc_tao, Gc_y, t.es_0, t.es_1, p.es_0, p.es_1)
        if (sum(grd)==0){
          norm_grd <- rep(0, length(beta))
        }else{norm_grd <- grd / sqrt(sum(grd^2))}
        for (j in 1:length(grd)) {
            v[j] <- beta2 * v[j] + (1 - beta2) * norm_grd[j]^2
            m[j] <- beta1 * m[j] + (1 - beta1) * norm_grd[j]
            v_hat[j] <- v[j] / (1 - beta2^(iter + 1))
            m_hat[j] <- m[j] / (1 - beta1^(iter + 1))
            beta[j] <- beta[j] + lr * m_hat[j] / (sqrt(v_hat[j]) + 1e-8)
        }
        beta_new <- beta

        beta <- beta_new
        iter <- iter + 1
        err <- sum(abs(grd))
      }
      obj_value <- aipw_U_stoc(beta, x, y, a, prob, tao, delta, Gc_tao, p.es_0, p.es_1)
      if (obj_value > best_obj_value) {
          best_obj_value <- obj_value
          best_beta <- beta
      }
  }
  # print(iter);print(err)
  return(best_beta)
}

solve_penalty_1 <- function(lambda, q_cons,x, y, a, prob, tao, delta, Gc_tao, Gc_y, t.es_0, t.es_1, p.es_0, p.es_1,
                            stp_sz = 0.001, tol = 1e-3, max_iter = 10000,num_starts=5,range=c(-10,10)){
  # lambda0 = 2
  lambda0 = lambda
  c=2
  lambda_seq = lambda0*c^(0:10)
  i=0; EU_hat=0;ET_hat <- 0
  while((ET_hat<q_cons) & (i<=0)){
    i=i+1
    lambda_temp = lambda_seq[i]
    
    beta_hat  = solve_p1(lambda_temp, q = q_cons, x = x, y = y, a = a, prob = prob, tao = tao, delta = delta, t.es_0, t.es_1, p.es_0, p.es_1,
                         Gc_tao=Gc_tao, Gc_y=Gc_y, stp_sz = 0.001, tol = 1e-3, max_iter = 10000,num_starts=5,range=c(-10,10))
    ET_hat <- aipw_T_stoc(beta_hat, x, y, a, prob, tao, delta, Gc_tao, Gc_y, t.es_0, t.es_1, p.es_0, p.es_1)
  }
  return(beta_hat)
}
################################ max T, s.t. U ##############################
aipw_U_stoc_q <- function(q_u, beta, x, y, a, prob, tao, delta, Gc_tao, p.es_0, p.es_1) {
  g <- 1 / (1 + exp(x %*% beta))
  wts <- a * g * 1 / prob + (1 - a) * (1 - g) * (1 / (1 - prob)) # 逆概率权重
  y <- as.numeric(y > tao) # tao标记
  val_tao <- wts * y * (1 - delta)
  val_tao <- val_tao / Gc_tao
  
  p.es <- g * p.es_1 + (1 - g) * p.es_0
  val_mi <- (wts - 1) * p.es # 标准aipw
  
  val <- mean(val_tao - val_mi)
  val <- -(val - q_u)^2
  return(val)
}
cptgrd_U_q <- function(q_u, beta, x, y, a, prob, tao, delta, Gc_tao,p.es_0, p.es_1) {
  d <- length(beta)
  m0 <- aipw_U_stoc_q(q_u, beta, x, y, a, prob, tao, delta, Gc_tao,p.es_0, p.es_1)
  grd <- rep(0, d)
  for (j in 1:d) {
    tmp <- rep(0, d)
    tmp[j] <- 0.01 * (2 * rbinom(1, size = 1, prob = 0.5) - 1)
    grd[j] <- (aipw_U_stoc_q(q_u, beta+tmp, x, y, a, prob, tao, delta, Gc_tao,p.es_0, p.es_1) - m0) / tmp[j] # ???ţ?
  }
  return(grd)
}
cptgrd_p2 <- function(beta, lambda, q, x, y, a, prob, tao, delta, Gc_tao, Gc_y, t.es_0, t.es_1, p.es_0, p.es_1){
  grd_t <- cptgrd_T(beta, x, y, a, prob, tao, delta, Gc_tao, Gc_y, t.es_0, t.es_1, p.es_0, p.es_1)
  # 罚函数分段求导
  if (aipw_U_stoc(beta, x, y, a, prob, tao, delta , Gc_tao, p.es_0, p.es_1)>=q){
    grd_u = 0
  } else{
    grd_u <- cptgrd_U_q(q, beta, x, y, a, prob, tao, delta, Gc_tao,p.es_0, p.es_1)
  }
  grd = grd_t + grd_u*lambda # 正负号（原始代码是负的。感觉是错的）
  return(grd)
}
# Given \lambda, 梯度下降法求解最优 \beta
solve_p2 <- function(lambda, q, x, y, a, prob, tao, delta, Gc_tao, Gc_y, t.es_0, t.es_1, p.es_0, p.es_1,
                     stp_sz = 0.001, tol = 1e-3, max_iter = 10000,
                     beta1 = 0.9, beta2 = 0.999,num_starts=5,range=c(-10,10)) {
  best_beta <- NULL;best_obj_value <- -Inf
  for (i in 1:num_starts) {
      if(i==1){
        beta <- rep(0, length(x[1,]))
      } else{
        beta <- runif(length(x[1,]),range[1],range[2])
      }
      iter <- 0;err <- 999
      lr <- stp_sz
      v <- rep(0, length(beta));v_hat <- rep(0, length(beta))
      m <- rep(0, length(beta));m_hat <- rep(0, length(beta))
      while ((err > tol) && (iter < max_iter)) {
        grd <- cptgrd_p2(beta,lambda,q, x, y, a, prob, tao, delta, Gc_tao, Gc_y, t.es_0, t.es_1, p.es_0, p.es_1)
        norm_grd <- grd / sqrt(sum(grd^2))
        for (j in 1:length(grd)) {
            v[j] <- beta2 * v[j] + (1 - beta2) * norm_grd[j]^2
            m[j] <- beta1 * m[j] + (1 - beta1) * norm_grd[j]
            v_hat[j] <- v[j] / (1 - beta2^(iter + 1))
            m_hat[j] <- m[j] / (1 - beta1^(iter + 1))
            beta[j] <- beta[j] + lr * m_hat[j] / (sqrt(v_hat[j]) + 1e-8)
        }
        beta_new <- beta

        beta <- beta_new
        iter <- iter + 1
        err <- sum(abs(grd))
      }
      obj_value <- aipw_T_stoc(beta, x, y, a, prob, tao, delta, Gc_tao, Gc_y, t.es_0, t.es_1, p.es_0, p.es_1)
      if (obj_value > best_obj_value) {
          best_obj_value <- obj_value
          best_beta <- beta
      }
  }
  # print(iter);print(err)
  return(best_beta)
}

solve_penalty_2 <- function(lambda, q_cons,x, y, a, prob, tao, delta, Gc_tao, Gc_y, t.es_0, t.es_1, p.es_0, p.es_1,
                            stp_sz = 0.001, tol = 1e-3, max_iter = 50000,
                            beta1 = 0.9, beta2 = 0.999,num_starts=5,range=c(-10,10)){
  # lambda0 = 2
  lambda0 = lambda
  c=2
  lambda_seq = lambda0*c^(0:10)
  i=0; EU_hat=0;ET_hat <- 0
  while((EU_hat<q_cons) & (i<=0)){
    i=i+1
    lambda_temp = lambda_seq[i]
    beta_hat  = solve_p2(lambda_temp, q = q_cons, x = x, y = y, a = a, prob = prob, tao = tao, delta = delta, t.es_0, t.es_1, p.es_0, p.es_1,
                         Gc_tao=Gc_tao, Gc_y=Gc_y, stp_sz = stp_sz, tol = tol, max_iter = max_iter,
                         num_starts=num_starts,range=range)
    EU_hat <- aipw_U_stoc(beta_hat, x, y, a, prob, tao, delta, Gc_tao, p.es_0, p.es_1)
  }
  return(beta_hat)
}

# 拉格朗日
# dual function
dual_1 <- function(beta,lambda,q, x, y, a, prob, tao, delta, Gc_tao, Gc_y) {
  val <- -ipw_T_stoc(beta, x, y, a, prob, tao, delta, Gc_tao, Gc_y) +
    lambda * (q - ipw_U_stoc(beta, x, y, a, prob, tao, delta, Gc_tao))
  return(val)
}
# cptgrd_d1(c(0,0,0,0,0),1,0.26,x, y, a, prob, tao, delta, Gc_tao, Gc_y )
# dual function gradient
cptgrd_d1 <- function(beta,lambda,q, x, y, a, prob, tao, delta, Gc_tao, Gc_y) {
  grd_u <- cptgrd_U(beta, x, y, a, prob, tao, delta, Gc_tao)
  grd_t <- cptgrd_T(beta, x, y, a, prob, tao, delta, Gc_tao, Gc_y)
  grd <- -grd_t - lambda * grd_u
  return(grd)
}

# Given \lambda, 梯度下降法求解最优 \beta
solve_d1 <- function(lambda, q, x, y, a, prob, tao, delta, Gc_tao=Gc_tao,Gc_y=Gc_y,
                     stp_sz = 0.01, tol = 1e-3, max_iter = 10000, result="value") {
  beta <- rep(0, length(x[1,])) # 初始点设为原点
  iter <- 0 # 迭代次数
  err <- 999 # 迭代误差
  # stp_sz: 步长
  # tol:误差容忍度
  while ((err > tol) && (iter < max_iter)) {
    grd <- cptgrd_d1(beta,lambda,q, x, y, a, prob, tao, delta, Gc_tao=Gc_tao,Gc_y=Gc_y)
    beta <- beta - stp_sz * grd
    iter <- iter + 1
    err <- sum(abs(grd))
    # print(iter)
  }
  # return(beta) # argmax
  if (result == "value") {
    return(dual_1(beta,lambda,q, x, y, a, prob, tao, delta, Gc_tao=Gc_tao,Gc_y=Gc_y))
  } else {
    return(beta)
  }
}
solve_lagragian_2 <- function(q_cons,x, y, a, prob, tao, delta, Gc_tao, Gc_y,
                              stp_sz = 0.001, tol = 1e-3, max_iter = 50000){
  lambda_opt <- optimize(solve_d1, interval = c(0, 100), q = q_cons,
                           x = x, y = y, a = a, prob = prob, tao = tao, delta = delta, Gc_tao=Gc_tao,Gc_y=Gc_y,
                           maximum = TRUE, tol = 0.1)$maximum
  beta_hat <- solve_d1(lambda_opt, q = q_cons, x = x, y = y, a = a, prob = prob, tao = tao, delta = delta, Gc_tao=Gc_tao,Gc_y=Gc_y,
                         result = "beta")
  return(beta_hat)
}
