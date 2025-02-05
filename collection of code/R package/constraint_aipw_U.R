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

solve_penalty_1 <- function(sigma, q_cons,data, tao=NUll, ps_formula, Q_T_formula, Q_cure_formula, Q_T_length, Tx,
                            stp_sz = 0.001, tol = 1e-3, max_iter = 10000,num_starts=5,range=c(-10,10)){
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
  Gc_y <- surv_est(T)
  Gc_y[Gc_y==0]=min(Gc_y[Gc_y!=0])
  
  # augment term estimation using gfcure(AFT)
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
  
  
  if (is.null(sigma)){
    sigma = 0.5
  }
  lambda0 = sigma
  c=2
  lambda_seq = lambda0*c^(0:10)
  i=0; EU_hat=0;ET_hat <- 0
  while((ET_hat<q_cons) && (i<=0)){
    i=i+1
    lambda_temp = lambda_seq[i]
    
    beta_hat  = solve_p1(lambda_temp, q = q_cons, x = Tx, y = T, a = A, prob = prob, tao = tao, delta = Status, t.es_0, t.es_1, p.es_0, p.es_1,
                         Gc_tao=Gc_tao, Gc_y=Gc_y, stp_sz = stp_sz, tol = tol, max_iter = max_iter,num_starts=num_starts,range=range)
    ET_hat <- aipw_T_stoc(beta_hat, Tx, T, A, prob, tao, Status, Gc_tao, Gc_y, t.es_0, t.es_1, p.es_0, p.es_1)
    EU_hat <- aipw_U_stoc(beta_hat, Tx, T, A, prob, tao, Status, Gc_tao, p.es_0, p.es_1)
  }
  return(list(beta_hat, EU_hat, ET_hat))
}

AIPW_U_stT <- function(n_bootstrap = 100, q_t, data, tao = NULL, ps_formula, Q_T_formula, 
                       Q_cure_formula, Q_T_length, Tx, stp_sz = 0.001, tol = 1e-3, 
                       max_iter = 10000, range = c(-10, 10), num_starts = 5) {
  
  ################## overall estimation ########
  start_time <- proc.time()
  result_overall <- solve_penalty_1(sigma = 0.5, q_cons = q_t, 
                            data = data, tao = tao, ps_formula = ps_formula, 
                            Q_T_formula = Q_T_formula, Q_cure_formula = Q_cure_formula, 
                            Q_T_length = Q_T_length, Tx = Tx, stp_sz = stp_sz, 
                            tol = tol, max_iter = max_iter, range = range, num_starts = num_starts)
  elapsed_time <- proc.time() - start_time
  print(paste("Overall estimation done. Time spent:", elapsed_time[3], "seconds"))
  
  ################## bootsrap ##################
  # Initialize storage for results
  beta_hat_list <- matrix(NA, nrow = n_bootstrap, ncol = length(Tx)+1)  # Assuming beta_hat has the same dimension as A
  eu_list <- numeric(n_bootstrap)
  et_list <- numeric(n_bootstrap)
  
  # Loop for Bootstrap resampling
  for (i in 1:n_bootstrap) {
    start_time <- proc.time()
    # Resampling
    df_bootstrap <- data[sample(1:nrow(data), replace = TRUE), ]
    
    # Call solve_penalty_1 function to get estimates
    result <- solve_penalty_1(sigma = 0.5, q_cons = q_t, 
                              data = df_bootstrap, tao = tao, ps_formula = ps_formula, 
                              Q_T_formula = Q_T_formula, Q_cure_formula = Q_cure_formula, 
                              Q_T_length = Q_T_length, Tx = Tx, stp_sz = stp_sz, 
                              tol = tol, max_iter = max_iter, range = range, num_starts = num_starts)
    
    # Store beta_hat, EU_hat, ET_hat
    beta_hat_list[i, ] <- result[[1]]  # Ensure you're using the correct element
    eu_list[i] <- result[[2]]
    et_list[i] <- result[[3]]
    print(i)
    elapsed_time <- proc.time() - start_time
    print(paste("Bootsrap Time spent:", elapsed_time[3], "seconds"))
  }
  
  # Calculate standard errors
  beta_hat_se <- apply(beta_hat_list, 2, sd)
  eu_se <- sd(eu_list)
  et_se <- sd(et_list)
  
  # Return estimates and standard errors
  return(list(beta_hat = result_overall[[1]], EU_hat = result_overall[[2]], 
              ET_hat = result_overall[[3]], beta_hat_se = beta_hat_se, 
              EU_se = eu_se, ET_se = et_se))
}
