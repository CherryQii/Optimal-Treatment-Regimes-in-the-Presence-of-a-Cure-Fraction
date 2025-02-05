aipw_U_rct <- function(x, y, a, prob, tao, delta, Gc_tao, p.es_0, p.es_1) {
  g <- rbinom(n = length(a), size = 1, prob = 0.5)
  wts <- a * g * 1 / prob + (1 - a) * (1 - g) * (1 / (1 - prob)) # 逆概率权重
  y <- as.numeric(y > tao) # tao标记
  val_tao <- wts * y * (1 - delta)
  val_tao <- val_tao / Gc_tao

  p.es <- g * p.es_1 + (1 - g) * p.es_0
  val_mi <- (wts - 1) * p.es # 标准aipw

  val <- mean(val_tao - val_mi)
  return(val)
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
