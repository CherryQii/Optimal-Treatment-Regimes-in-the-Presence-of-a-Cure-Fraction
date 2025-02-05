rm(list=ls())
library(survival);library(survminer);library(mixcure);library(ROCR);library(rgenoud);library(gfcure)
set.seed(2023)
df <- read.csv("seer_esophageal.csv")
df <- read.csv("8.19.csv")
# read .Rda
load("样本数据.Rda")
df = df_store
#################################### data clean ########################
# Stage
table(df$Stage) # 1=Localized,2=Regional,7=Distant
df <- df[df$Stage!=9 & df$Stage!=14,]
df <- df[df$Status != 8 , ]
# histology
table(df$Histology)# 2:squamous cell; 5:adenomas & adenocarcinomas; others
df$histology <- 1
df$histology[df$Histology==5] <- 2
df$histology[(df$Histology!=2) & (df$Histology !=5)] <- 3
table(df$histology)
# surgery performed
df <- df[df$surgery==0,]
df$A1 <- 0
df$A1[df$Surg.Rad.Seq!=0] <- 1

# convert to dummy variable
df$stage_r <- as.numeric(df$Stage==2) # stage_r=stage2
df$stage_d <- as.numeric(df$Stage==7) # stage_d = stage3
df$histology_sq <- as.numeric(df$Histology==2)
df$histology_ad <- as.numeric(df$Histology==5)
# age
df$age_1 <- (df$Age<=45)
df$age_2 <- (df$Age>45 & df$Age<=65)
df$age_3 <- (df$Age>65)
df$age_class <- 0
df$age_class[df$age_2==1] <- 1
df$age_class[df$age_3==1] <- 2
table(df$age_class)
df$age <- (df$Age-min(df$Age))/(max(df$Age) - min(df$Age))
# remove the patients with 0 survival time
df = df[df$Survival.months!=0,]

table(df$Stage);table(df$Sex);table(df$A1)# 1:Distant, 2:Localized, 3:Regional
colnames(df)
dim(df)
################################### overall K-M plot/estimated overall cure rate #########################
exp(-summary(km)$cumhaz)
km <- survfit(Surv(Survival.months, Status) ~ 1, data = df)
# add a horizontal line at 0.454 in ggsurvplot
p = ggsurvplot(km, data = df,xlab = "Months", ylab = "Survival Probability", conf.int = FALSE, color="red")
p$plot = p$plot + geom_hline(yintercept = km$surv[length(km$surv)], linetype = "dashed", color = "blue")
print(p)
summary(km)
km$surv[length(km$surv)] # estimated overall cure rate
sum(1-df$Status)
sum(1-df$Status)/7080 # sample censoring rate
max(df$Survival.months[df$Status==1])
sum(df$Survival.months>190)
################################### fig-1: K-M plot #########################
km1 <- survfit(Surv(Survival.months, Status) ~ A1, data = df)
ggsurvplot(km1, data = df,legend.title="Treatment",legend.labs=c("A=0","A=1"),xlab = "Months", ylab = "Survival Probability", conf.int = TRUE)
###############################################################
attach(df)
# x <- cbind(1,age, stage_r,stage_d)
x <- cbind(1, Age, stage_r,stage_d)
y <- Survival.months
a <- A1
############################### PS model ###############################
fit_pi <- glm(a~Age+Sex+Stage+Histology, family = binomial(link = "logit"),data = df)
# fit_pi <- glm(a~Age+Sex+stage_r+stage_d+histology_sq+histology_ad, family = binomial(link = "logit"),data = df)
prob <- predict(fit_pi, type = "response")
pred <- prediction(prob, a);auc <- performance(pred, "auc");auc <- unlist(slot(auc, "y.values"));auc
delta <- Status

# tao <- max(y[delta==1]) # 182
tao <- 120
sum(y>tao)
detach(df)
##################### estimated survival function of censoring time ########################
fit_Gc <- survfit(Surv(Survival.months, 1-Status)~1, data=df)
surv_est <- stepfun(fit_Gc$time, c(1,fit_Gc$surv))
Gc_tao <- surv_est(tao)
Gc_y <- surv_est(y)
Gc_y[Gc_y==0]=min(Gc_y[Gc_y!=0])
#################### without considering cure rate ############################
ipw_ET_without_cure <- function(eta, x, y, a, prob, delta, Gc_y){
  g <- as.numeric(x %*% eta > 0)
  wts <- a* g * 1 / prob + (1-a) * (1 - g) * (1 / (1 - prob)) # 逆概率权重
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
set.seed(2024)
result_ET_without_cure = GA_ipw_ET(x, y, a, prob, delta, Gc_y)
result_ET_without_cure
ipw_ET_without_cure(result_ET_without_cure[[1]], x, y, a, prob, delta, Gc_y)
ipw_U(result_ET_without_cure[[1]], x, y, a, prob, tao, delta, Gc_tao)
ipw_T(result_ET_without_cure[[1]], x, y, a, prob, tao, delta, Gc_y, Gc_tao)

set.seed(2024)
result_EU = GA_U_ipw(x, y, a, prob, tao, delta, Gc_tao)
ipw_ET_without_cure(result_EU, x, y, a, prob, delta, Gc_y)
ipw_U(result_EU, x, y, a, prob, tao, delta, Gc_tao)
ipw_T(result_EU, x, y, a, prob, tao, delta, Gc_y, Gc_tao)

set.seed(2024)
result_ET = GA_T_ipw(x, y, a, prob, tao, delta, Gc_y, Gc_tao)
ipw_ET_without_cure(result_ET, x, y, a, prob, delta, Gc_y)
ipw_U(result_ET, x, y, a, prob, tao, delta, Gc_tao)
ipw_T(result_ET, x, y, a, prob, tao, delta, Gc_y, Gc_tao)

ipw_U_rct(x, y, a, prob, tao, delta, Gc_tao)
ipw_T_rct(x, y, a, prob, tao, delta, Gc_y, Gc_tao)
ipw_U_samp(x, y, a, prob, tao, delta, Gc_tao)
########################### AFT model ############################
fit_aft1 <- gfcure(Surv(Survival.months, Status)~Sex*A1 + Age*A1 + stage_r*A1 + stage_d*A1+histology_sq*A1 + histology_ad*A1,
              ~Sex*A1 + Age*A1 + stage_r*A1 + stage_d*A1+histology_sq*A1 + histology_ad*A1,
              data= df, dist = "lognormal")
fit_aft2 <- gfcure(Surv(Survival.months, Status)~ Age*A1 + stage_r*A1 + stage_d*A1,
                  ~Age*A1 + stage_r*A1 + stage_d*A1,
                  data= df, dist = "lognormal")
gamma_p_hat <- -fit_aft2$coef[10:17]
logit_p <- x%*%gamma_p_hat[c(1,2,4,5)]; p.es_0 <- exp(logit_p)/(1+exp(logit_p)) 
logit_p <- x%*%gamma_p_hat[c(1,2,4,5)] + x%*%gamma_p_hat[c(3,6,7,8)]; p.es_1 <- exp(logit_p)/(1+exp(logit_p))
gamma_t_hat <- fit_aft2$coef[2:9]
t.es_0 <- exp(x%*%gamma_t_hat[c(1,2,4,5)]) * exp(exp(fit_aft2$coef[1])^2/2)
t.es_1 <- exp(x%*%gamma_t_hat[c(1,2,4,5)] + 1* x %*%gamma_t_hat[c(3,6,7,8)]) * exp(exp(fit_aft2$coef[1])^2/2)

eta_scale <- function(eta){
  eta <- eta/sum(abs(eta))
  return(eta)
}


source("function\\aipw-U.R")
source("function\\aipw-T.R")
# RCT
set.seed(2023)
aipw_U_rct(x, y, a, prob, tao, delta, Gc_tao, p.es_0,p.es_1)
aipw_T_rct(x, y, a, prob, tao, delta, Gc_y, Gc_tao,t.es_0, t.es_1, p.es_0, p.es_1)
# D-aipw-U: maxU
set.seed(2023)
eta1 <- GA_U_aipw(x, y, a, prob, tao, delta, Gc_tao, p.es_0,p.es_1)
aipw_U(eta1,x, y, a, prob, tao, delta, Gc_tao, p.es_0,p.es_1)
aipw_T(eta1, x, y, a, prob, tao, delta, Gc_y, Gc_tao,t.es_0, t.es_1, p.es_0, p.es_1)
# D-aipw-T: maxT
set.seed(2023)
eta2 <- GA_T_aipw(x, y, a, prob, tao, delta, Gc_y, Gc_tao,t.es_0,t.es_1, p.es_0, p.es_1)
aipw_U(eta2,x, y, a, prob, tao, delta, Gc_tao, p.es_0,p.es_1)
aipw_T(eta2, x, y, a, prob, tao, delta, Gc_y, Gc_tao,t.es_0, t.es_1, p.es_0, p.es_1)

sum(as.numeric(x %*% eta1 > 0))

apply(rbind(eta1,eta2),1,eta_scale)


source("function\\aipw-constraint.R")
# max U
set.seed(2023)
beta1 <- solve_U_aipw(x, y, a, prob, tao, delta, Gc_tao,p.es_0, p.es_1,
                      stp_sz = 0.001,tol = 1e-3, max_iter = 20000, range=c(-15,15),num_starts = 25)
beta1
aipw_U_stoc(beta1,x, y, a, prob, tao, delta, Gc_tao, p.es_0, p.es_1)
aipw_T_stoc(beta1,x, y, a, prob, tao, delta, Gc_tao, Gc_y, t.es_0, t.es_1, p.es_0, p.es_1)
# max T
set.seed(2023)
beta2 <- solve_T_aipw(x, y, a, prob, tao, delta, Gc_tao, Gc_y, t.es_0, t.es_1, p.es_0, p.es_1,
                      stp_sz = 0.01,tol = 1e-3, max_iter = 20000,range=c(-10,10),num_starts = 25)
aipw_U_stoc(beta2,x, y, a, prob, tao, delta, Gc_tao, p.es_1, p.es_0)
aipw_T_stoc(beta2,x, y, a, prob, tao, delta, Gc_tao, Gc_y, t.es_0, t.es_1, p.es_0, p.es_1)

# write.csv(rbind(eta1,eta2,eta11,eta22,beta1,beta2),
#           file="no_constraint_result.csv")
########################################## Constraint ############################################
## penalty
set.seed(2023)
beta_3_54 <- solve_penalty_1(2, q_cons=54,x, y, a, prob, tao, delta, Gc_tao, Gc_y,t.es_0, t.es_1, p.es_0, p.es_1,
                             stp_sz = 0.01, tol = 1e-3, max_iter = 20000,range=c(-10,10),num_starts = 5)
beta_3_55 <- solve_penalty_1(0.1, q_cons=55,x, y, a, prob, tao, delta, Gc_tao, Gc_y,t.es_0, t.es_1, p.es_0, p.es_1,
                             stp_sz = 0.01, tol = 1e-3, max_iter = 20000,range=c(-10,10),num_starts = 5)
beta_3_56 <- solve_penalty_1(0.1, q_cons=56,x, y, a, prob, tao, delta, Gc_tao, Gc_y,t.es_0, t.es_1, p.es_0, p.es_1,
                             stp_sz = 0.01, tol = 1e-3, max_iter = 20000,range=c(-10,10),num_starts = 5)
aipw_U_stoc(beta_3_54, x, y, a, prob, tao, delta, Gc_tao, p.es_0,p.es_1)
aipw_T_stoc(beta_3_54, x, y, a, prob, tao, delta, Gc_tao, Gc_y, t.es_0, t.es_1, p.es_0, p.es_1)

beta_3_57 <- solve_penalty_1(0.1, q_cons=57,x, y, a, prob, tao, delta, Gc_tao, Gc_y,t.es_0, t.es_1, p.es_0, p.es_1,
                             stp_sz = 0.01, tol = 1e-3, max_iter = 20000,range=c(-10,10),num_starts = 5)
aipw_U_stoc(beta_3_54, x, y, a, prob, tao, delta, Gc_tao, p.es_0,p.es_1)
aipw_T_stoc(beta_3_54, x, y, a, prob, tao, delta, Gc_tao, Gc_y, t.es_0, t.es_1, p.es_0, p.es_1)

beta_3_58 <- solve_penalty_1(0.1, q_cons=58,x, y, a, prob, tao, delta, Gc_tao, Gc_y,t.es_0, t.es_1, p.es_0, p.es_1,
                             stp_sz = 0.01, tol = 1e-3, max_iter = 20000,range=c(-10,10),num_starts = 5)
aipw_U_stoc(beta_3_58, x, y, a, prob, tao, delta, Gc_tao, p.es_0,p.es_1)
aipw_T_stoc(beta_3_58, x, y, a, prob, tao, delta, Gc_tao, Gc_y, t.es_0, t.es_1, p.es_0, p.es_1)

beta_3_59 <- solve_penalty_1(0.1, q_cons=59,x, y, a, prob, tao, delta, Gc_tao, Gc_y,t.es_0, t.es_1, p.es_0, p.es_1,
                             stp_sz = 0.01, tol = 1e-3, max_iter = 20000,range=c(-10,10),num_starts = 5)
aipw_U_stoc(beta_3_59, x, y, a, prob, tao, delta, Gc_tao, p.es_0,p.es_1)
aipw_T_stoc(beta_3_59, x, y, a, prob, tao, delta, Gc_tao, Gc_y, t.es_0, t.es_1, p.es_0, p.es_1)

beta_3_60 <- solve_penalty_1(0.1, q_cons=60,x, y, a, prob, tao, delta, Gc_tao, Gc_y,t.es_0, t.es_1, p.es_0, p.es_1,
                             stp_sz = 0.01, tol = 1e-3, max_iter = 20000,range=c(-10,10),num_starts = 5)
aipw_U_stoc(beta_3_60, x, y, a, prob, tao, delta, Gc_tao, p.es_0,p.es_1)
aipw_T_stoc(beta_3_60, x, y, a, prob, tao, delta, Gc_tao, Gc_y, t.es_0, t.es_1, p.es_0, p.es_1)

beta_3_61 <- solve_penalty_1(0.1, q_cons=61,x, y, a, prob, tao, delta, Gc_tao, Gc_y,t.es_0, t.es_1, p.es_0, p.es_1,
                             stp_sz = 0.01, tol = 1e-3, max_iter = 20000,range=c(-10,10),num_starts = 5)
aipw_U_stoc(beta_3_61, x, y, a, prob, tao, delta, Gc_tao, p.es_0,p.es_1)
aipw_T_stoc(beta_3_61, x, y, a, prob, tao, delta, Gc_tao, Gc_y, t.es_0, t.es_1, p.es_0, p.es_1)


set.seed(2023)
beta_3_list <- NULL
for (q_t in seq(50,60,2)){
  beta_3_temp <- solve_penalty_1(0.5, q_cons=q_t,x, y, a, prob, tao, delta, Gc_tao, Gc_y,t.es_0, t.es_1, p.es_0, p.es_1,
                                 stp_sz = 0.01, tol = 1e-3, max_iter = 20000,range=c(-10,10),num_starts = 1)
  beta_3_list <- rbind(beta_3_list, beta_3_temp)
}

beta_3_50 <- solve_penalty_1(0.5, q_cons=50,x, y, a, prob, tao, delta, Gc_tao, Gc_y,t.es_0, t.es_1, p.es_0, p.es_1,
                             stp_sz = 0.01, tol = 1e-3, max_iter = 20000,range=c(-10,10),num_starts = 1)
beta_3_52 <- solve_penalty_1(0.5, q_cons=52,x, y, a, prob, tao, delta, Gc_tao, Gc_y,t.es_0, t.es_1, p.es_0, p.es_1,
                             stp_sz = 0.01, tol = 1e-3, max_iter = 20000,range=c(-10,10),num_starts = 1)
aipw_U_stoc(beta_3_50, x, y, a, prob, tao, delta, Gc_tao, p.es_0,p.es_1)
aipw_T_stoc(beta_3_52, x, y, a, prob, tao, delta, Gc_tao, Gc_y, t.es_0, t.es_1, p.es_0, p.es_1)

beta_3_list <- rbind(beta_3_42,beta_3_44,beta_3_46,beta_3_48)
apply(beta_3_list, 1, aipw_U_stoc,x, y, a, prob, tao, delta, Gc_tao, p.es_0,p.es_1)
apply(beta_3_list, 1, aipw_T_stoc,x, y, a, prob, tao, delta, Gc_tao, Gc_y, t.es_0, t.es_1, p.es_0, p.es_1)
write.csv(beta_3_list, file="D:\\论文\\Cure Model + ITR\\数据\\Esophageal Cancer\\beta_3(st T).csv")

#################################### max T, s.t. U ##########################
## penalty
# 0.03 * 1000=30
set.seed(2023)
beta_4_0.47 <- solve_penalty_2(1500, q_cons=0.47,x, y, a, prob, tao, delta, Gc_tao, Gc_y,t.es_0, t.es_1, p.es_0, p.es_1,
                               stp_sz = 0.001, tol = 1e-3, max_iter = 20000,range=c(-10,10),num_starts = 5)
aipw_U_stoc(beta_4_0.47, x, y, a, prob, tao, delta, Gc_tao, p.es_0,p.es_1)
aipw_T_stoc(beta_4_0.47, x, y, a, prob, tao, delta, Gc_tao, Gc_y, t.es_0, t.es_1, p.es_0, p.es_1)

beta_4_0.48 <- solve_penalty_2(20000, q_cons=0.48,x, y, a, prob, tao, delta, Gc_tao, Gc_y,t.es_0, t.es_1, p.es_0, p.es_1,
                               stp_sz = 0.001, tol = 1e-3, max_iter = 20000,range=c(-10,10),num_starts = 1)
aipw_U_stoc(beta_4_0.48, x, y, a, prob, tao, delta, Gc_tao, p.es_0,p.es_1)
aipw_T_stoc(beta_4_0.48, x, y, a, prob, tao, delta, Gc_tao, Gc_y, t.es_0, t.es_1, p.es_0, p.es_1)
beta_4_0.49 <- solve_penalty_2(20000, q_cons=0.49,x, y, a, prob, tao, delta, Gc_tao, Gc_y,t.es_0, t.es_1, p.es_0, p.es_1,
                               stp_sz = 0.001, tol = 1e-3, max_iter = 20000,range=c(-10,10),num_starts = 1)
aipw_U_stoc(beta_4_0.49, x, y, a, prob, tao, delta, Gc_tao, p.es_0,p.es_1)
aipw_T_stoc(beta_4_0.49, x, y, a, prob, tao, delta, Gc_tao, Gc_y, t.es_0, t.es_1, p.es_0, p.es_1)
beta_4_0.50 <- solve_penalty_2(20000, q_cons=0.50,x, y, a, prob, tao, delta, Gc_tao, Gc_y,t.es_0, t.es_1, p.es_0, p.es_1,
                               stp_sz = 0.001, tol = 1e-3, max_iter = 20000,range=c(-10,10),num_starts = 1)
aipw_U_stoc(beta_4_0.50, x, y, a, prob, tao, delta, Gc_tao, p.es_0,p.es_1)
aipw_T_stoc(beta_4_0.50, x, y, a, prob, tao, delta, Gc_tao, Gc_y, t.es_0, t.es_1, p.es_0, p.es_1)
beta_4_0.51 <- solve_penalty_2(20000, q_cons=0.51,x, y, a, prob, tao, delta, Gc_tao, Gc_y,t.es_0, t.es_1, p.es_0, p.es_1,
                               stp_sz = 0.001, tol = 1e-3, max_iter = 20000,range=c(-10,10),num_starts = 1)
aipw_U_stoc(beta_4_0.51, x, y, a, prob, tao, delta, Gc_tao, p.es_0,p.es_1)
aipw_T_stoc(beta_4_0.51, x, y, a, prob, tao, delta, Gc_tao, Gc_y, t.es_0, t.es_1, p.es_0, p.es_1)
beta_4_0.52 <- solve_penalty_2(20000, q_cons=0.52,x, y, a, prob, tao, delta, Gc_tao, Gc_y,t.es_0, t.es_1, p.es_0, p.es_1,
                               stp_sz = 0.001, tol = 1e-3, max_iter = 20000,range=c(-10,10),num_starts = 1)
aipw_U_stoc(beta_4_0.55, x, y, a, prob, tao, delta, Gc_tao, p.es_0,p.es_1)
aipw_T_stoc(beta_4_0.52, x, y, a, prob, tao, delta, Gc_tao, Gc_y, t.es_0, t.es_1, p.es_0, p.es_1)
beta_4_0.53 <- solve_penalty_2(20000, q_cons=0.53,x, y, a, prob, tao, delta, Gc_tao, Gc_y,t.es_0, t.es_1, p.es_0, p.es_1,
                               stp_sz = 0.001, tol = 1e-3, max_iter = 20000,range=c(-10,10),num_starts = 1)
beta_4_0.54 <- solve_penalty_2(20000, q_cons=0.54,x, y, a, prob, tao, delta, Gc_tao, Gc_y,t.es_0, t.es_1, p.es_0, p.es_1,
                               stp_sz = 0.001, tol = 1e-3, max_iter = 20000,range=c(-10,10),num_starts = 1)
beta_4_0.55 <- solve_penalty_2(20000, q_cons=0.55,x, y, a, prob, tao, delta, Gc_tao, Gc_y,t.es_0, t.es_1, p.es_0, p.es_1,
                               stp_sz = 0.001, tol = 1e-3, max_iter = 20000,range=c(-10,10),num_starts = 1)

############################### figure: d-2 ###############################
# beta_3 <- rbind(beta_3_54,beta_3_55, beta_3_56,beta_3_57,beta_3_58,beta_3_59,beta_3_60)
# write.csv(beta_3, file="D:\\论文\\Cure Model + ITR\\数据\\beta_3.csv")
beta_3 <- read.csv("D:\\论文\\Cure Model + ITR\\数据\\beta_3.csv")
x <- cbind(1,age, stage_r,stage_d);y <- Survival.months;a <- A1
beta_3 <- beta_3[,2:5]
df1 <- NULL
df1$U <- apply(beta_3, 1, aipw_U_stoc, x = x, y = y, a = a, prob = prob, tao = tao, delta = delta, Gc_tao = Gc_tao, p.es_1=p.es_1, p.es_0=p.es_0)
df1$T <- apply(beta_3, 1, aipw_T_stoc, x = x, y = y, a = a, prob = prob, tao = tao, delta = delta,Gc_tao = Gc_tao, Gc_y = Gc_y, t.es_0=t.es_0, t.es_1=t.es_1,p.es_0=p.es_0,p.es_1=p.es_1)
df1$x <- c(54,55,56,57,58,59,60)
df1 <- as.data.frame(df1);row.names(df1)<-NULL
df1
U_range <- range(df1$U)
T_range <- range(df1$T)
transformation <- function(x) {
  scaled_x <- (x - U_range[1]) / (U_range[2] - U_range[1])  # Scale U to [0, 1]
  mapped_x <- scaled_x * (T_range[2] - T_range[1]) + T_range[1]  # Map to T range
  return(mapped_x)
}

ylim.T <- c(45,64)
ylim.U <- c(0.4,0.55)
b <- diff(ylim.T)/diff(ylim.U)
a <- ylim.T[1] - b*ylim.U[1]
# plot double coordinate
library(ggplot2)
academic_palette <- scale_color_brewer(palette = "Set2")
ggplot(df1, aes(x = x))+
  theme_bw()+ # bg
  geom_line(aes(y = T, color='T'), size = 1.5)+
  geom_line(aes(y = a + U * b, color='U'), size = 1.5, linetype='dashed')+
  scale_color_discrete(labels=c(expression(tilde(M)[aipw](hat(g)[alpha])),expression(tilde(V)[aipw](hat(g)[alpha]))))+
  scale_x_continuous(expression(q[t]), breaks = seq(54,61,by=1)) +
  scale_y_continuous(
    name = expression(tilde(M)[aipw](hat(g)[alpha])),  # Label for the left y-axis
    limits = c(35,64),
    sec.axis = sec_axis(~ (.-a)/b, name = expression(tilde(V)[aipw](hat(g)[alpha]))),  # Define the right y-axis for U
  )+
  theme(axis.title.y = element_text(angle = 0,vjust=0.5)) +
  theme(axis.title.y.right = element_text(angle = 0,vjust=0.5))+
  geom_point(aes(y = T, color = "T"), size = 4.5, shape = 24) +
  geom_point(aes(y = a + U * b, color = "U"), size = 4.5, shape = 21)+
  geom_hline(yintercept = 39.691, linetype = "dotted",size=1.5)+
  geom_hline(yintercept = a + 0.404*b, linetype = "dotdash",size=1.5)+
  annotate("text", x = 60, y = 38, label = expression(tilde(M)[aipw](g^RCT)), hjust = 1) + # Label for the dotted line
  annotate("text", x = 60, y = a + 0.393 * b, label = expression(tilde(V)[aipw](g^RCT)), hjust = 1) + # Label for the dotdash line
  theme(legend.position = "top",legend.title = element_blank())+
  guides(colour = guide_legend(override.aes = list(shape = c(24, 21))))
    
  
# ggplot(df1, aes(x = x))+
#   geom_line(aes(y = T,color = 'T'), size = 1)+
#   geom_line(aes(y = a + U * b, color = "U"), size = 1, linetype='dashed') +
#   theme_bw() +
#   # labs( color = c("Estimator",'asd')) +
#   scale_x_continuous(expression(q[t]), breaks = seq(54,61,by=1)) +
#   scale_y_continuous(
#     name = expression(hat(M)[aipw](tilde(g)[alpha])),  # Label for the left y-axis
#     limits = c(35,64),
#     sec.axis = sec_axis(~ (.-a)/b, name = expression(hat(V)[aipw](tilde(g)[alpha]))),  # Define the right y-axis for U
#   )+
#   geom_hline(yintercept = 39.691, linetype = "solid", color = "#085708")+
#   geom_hline(yintercept = a + 0.404*b, linetype = "dashed", color = "red")+
#   geom_point(aes(y = T, color = "T"), size = 3, shape = 24, fill = alpha("#085708", 0.3), color='#085708') +
#   geom_point(aes(y = a + U * b, color = "U"), size = 3, shape = 21, fill = alpha("red",0.3))+
#   academic_palette+
#   theme(legend.position = "top",legend.title = element_blank()) +
#   theme(axis.title.y = element_text(angle = 0,vjust=0.5)) +
#   theme(axis.title.y.right = element_text(angle = 0,vjust=0.5))+
#   guides(colour = guide_legend(override.aes = list(shape = c(24, 21), fill = c(alpha("#085708", 0.3), alpha("red", 0.3)))))

# d-3
x <- cbind(1,age, stage_r,stage_d);y <- Survival.months;a <- A1
# beta_4 <- rbind(beta_4_0.47,beta_4_0.48,beta_4_0.49,beta_4_0.50,beta_4_0.51,beta_4_0.52,beta_4_0.53,beta_4_0.54,beta_4_0.55)
# write.csv(beta_4, file="D:\\论文\\Cure Model + ITR\\数据\\beta_4.csv")
beta_4 <- read.csv("D:\\论文\\Cure Model + ITR\\数据\\beta_4.csv")
beta_4 <- beta_4[, 2:5]
df2 <- NULL
df2$U <- apply(beta_4, 1, aipw_U_stoc, x = x, y = y, a = a, prob = prob, tao = tao, delta = delta, Gc_tao = Gc_tao, p.es_1=p.es_1, p.es_0=p.es_0)
df2$T <- apply(beta_4, 1, aipw_T_stoc, x = x, y = y, a = a, prob = prob, tao = tao, delta = delta,Gc_tao = Gc_tao, Gc_y = Gc_y, t.es_0=t.es_0, t.es_1=t.es_1,p.es_0=p.es_0,p.es_1=p.es_1)
df2$x <- c(0.47,0.48,0.49,0.50,0.51,0.52,0.53,0.54,0.55)
df2 <- as.data.frame(df2);row.names(df2)<-NULL
df2
ylim.T <- c(45,64)
ylim.U <- c(0.4,0.55)
b <- diff(ylim.U)/diff(ylim.T)
a <- ylim.U[1] - b*ylim.T[1]
# plot double coordinate
academic_palette <- scale_color_brewer(palette = "Set2")
ggplot(df2, aes(x = x)) +
  geom_line(aes(y = U, color = "U"), size = 1.5, linetype = 'dashed') +  # Switched the order
  geom_line(aes(y = a + T*b, color = "T"), size = 1.5) +  # Switched the order
  scale_color_discrete(labels=c(expression(tilde(M)[aipw](hat(g)[beta])),expression(tilde(V)[aipw](hat(g)[beta]))))+
  theme_bw() +
  theme(legend.position = "top",legend.title = element_blank()) +
  scale_x_continuous(expression(q[u]), breaks = seq(0.4,0.69,by=0.02)) +
  scale_y_continuous(
    name = expression(tilde(V)[aipw](hat(g)[beta])),  # Label for the left y-axis
    limits = c(0.35,0.55),
    sec.axis = sec_axis(~ (.-a)/b, name = expression(tilde(M)[aipw](hat(g)[beta]))),  # Define the right y-axis for U
  )+
  theme(axis.title.y = element_text(angle = 0,vjust=0.5)) +
  theme(axis.title.y.right = element_text(angle = 0,vjust=0.5))+
  geom_hline(yintercept = 0.404, linetype = "dotdash",size=1.5)+
  geom_hline(yintercept = a + 39.691*b, linetype = "dotted",size=1.5)+
  annotate("text", x = 0.55, y = a + 39 * b, label = expression(tilde(M)[aipw](g^RCT)), hjust = 1) + # Label for the dotted line
  annotate("text", x = 0.55, y = 0.398, label = expression(tilde(V)[aipw](g^RCT)), hjust = 1) + # Label for the dotdash line
  geom_point(aes(y = U, color = "U"), size = 4.5, shape = 21,label = expression(tilde(V)[aipw](tilde(g)[beta]))) +
  geom_point(aes(y = a + T * b, color = "T"), size = 4.5, shape = 24,label = expression(tilde(M)[aipw](tilde(g)[beta])))+
  guides(colour = guide_legend(override.aes = list(shape = c(24, 21))))

############################### 分配概率图 #################################
Stage <- as.factor(Stage);levels(Stage)=c('Localized','Regional','Distant')

beta_3_maxU <- x %*% unlist(beta_3[1,])
g_opt_prob <- 1 / (1 + exp(beta_3_maxU))
df_3 <- as.data.frame(cbind(g_opt_prob, Age, Stage));df_3$Stage <- Stage
ggplot(df_3, aes(x = Age, y = g_opt_prob, shape = Stage, color = Stage)) +
  geom_point(size=2)+  # 添加散点图
  labs(color = "Tumor Stage", shape='Tumor Stage',y='Pr(A=1|X)')+
  scale_shape_manual(values = c(1, 9, 3)) 

beta_3_maxT <- x %*% unlist(beta_3[7,])
g_opt_prob <- 1 / (1 + exp(beta_3_maxT))
df_3 <- as.data.frame(cbind(g_opt_prob, Age, Stage));df_3$Stage <- Stage
ggplot(df_3, aes(x = Age, y = g_opt_prob, shape = Stage, color = Stage)) +
  geom_point(size=2)+  # 添加散点图
  labs(color = "Tumor Stage", shape='Tumor Stage', y='Pr(A=1|X)')+
  scale_shape_manual(values = c(1, 9, 3)) 
  # theme(axis.title.y = element_text(angle = 0,vjust=0.5))

beta_3_qt_57 <- x %*% unlist(beta_3[5,])
g_opt_prob <- 1 / (1 + exp(beta_3_qt_57))
df_3 <- as.data.frame(cbind(g_opt_prob, Age, Stage));df_3$Stage <- Stage
ggplot(df_3, aes(x = Age, y = g_opt_prob, shape = Stage, color = Stage)) +
  geom_point(size=2)+  # 添加散点图
  labs(color = "Tumor Stage", shape='Tumor Stage',y='Pr(A=1|X)')+
  scale_shape_manual(values = c(1, 9, 3)) 

beta_4_qu_52 <- x %*% unlist(beta_4[5,])
g_opt_prob <- 1 / (1 + exp(beta_4_qu_52))
df_3 <- as.data.frame(cbind(g_opt_prob, Age, Stage));df_3$Stage <- Stage
ggplot(df_3, aes(x = Age, y = g_opt_prob, shape = Stage, color = Stage)) +
  geom_point(size=2)+  # 添加散点图
  labs(color = "Tumor Stage", shape='Tumor Stage',y='Pr(A=1|X)')+
  scale_shape_manual(values = c(1, 9, 3)) 
