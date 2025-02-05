aipw_T_rct <- function(x, y, a, prob, tao, delta, Gc_y, Gc_tao, t.es_0, t.es_1, p.es_0, p.es_1) {
    g <- rbinom(n = length(a), size = 1, prob = 0.5)
    wts <- a * g * 1 / prob + (1 - a) * (1 - g) * (1 / (1 - prob)) # 逆概率权重
    t.es <- g * t.es_1*p.es_0 + (1 - g) * t.es_0*p.es_1
    # .T <- ifelse(wts * y * delta / Gc_y > tao, tao, wts * y * delta / Gc_y)
    .T <- wts * y * delta / Gc_y
    v_t <- mean(.T + (1 - wts) * t.es)

    v_u <- aipw_U_rct(x, y, a, prob, tao, delta, Gc_tao, p.es_0, p.es_1)
    val_t <- v_t / (1 - v_u)

    return(val_t)
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
