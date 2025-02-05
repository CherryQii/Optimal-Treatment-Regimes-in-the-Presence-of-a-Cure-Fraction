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

