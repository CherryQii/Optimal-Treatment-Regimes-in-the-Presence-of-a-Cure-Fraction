library(survival);library(rgenoud);library(gfcure);library(dplyr)
#################################### data load ########################
df <- read.csv("seer_esophageal_clean.csv")
df <- df %>%
  rename(
    T = Survival.months,
    Status = Status,
    A = A1
  )
#################################### Propensity Score ########################
fit_pi <- glm(A~age+Sex+Stage+Histology, family = binomial(link = "logit"),data = df)
prob <- predict(fit_pi, type = "response")
pred <- prediction(prob, df$A);auc <- performance(pred, "auc");auc <- unlist(slot(auc, "y.values"));auc
############################  initialization  #########################
ps_formula <- "A ~ Age + Sex + Stage + Histology"
Q_T_formula <- "Surv(T, Status) ~ Age * A + stage_r * A + stage_d * A"
Q_cure_formula <- "~ Age * A + stage_r * A + stage_d * A"
Tx = c("Age", "stage_r", "stage_d")
Q_T_length = 7

source("Function\\shared_function.R")

source("Function\\constraint_aipw_U.R")
#################################### AIPW_RCT ########################
source("Function\\NoCons_AIPW_RCT.R")
set.seed(2023)
AIPW_RCT(df, tao = 120, ps_formula, Q_T_formula, Q_cure_formula, Q_T_length, Tx)


#################################### AIPW_U_free ########################
# set.seed(2023)
# source("Function\\AIPW_U_free.R")
# AIPW_U_free(df, ps_formula, Q_T_formula, Q_cure_formula, Q_T_length,Tx, tao=120)
#################################### AIPW_maxU_stT ########################
set.seed(2023)
library(foreach);library(doParallel);library(doRNG);library(parallel)
# detectCores()
cl <- makeCluster(10)
registerDoParallel(cl)
clusterEvalQ(cl, library(rgenoud));clusterEvalQ(cl, library(survival));clusterEvalQ(cl, library(gfcure))
registerDoRNG(2023)
q_t_list <- seq(54, 60, 1)
results <- foreach(q_t = q_t_list, .combine = 'list', .multicombine = TRUE) %dopar% {
  print(paste("Processing q_t =", q_t))
  
  # 调用 AIPW_U_stT 进行估计
  AIPW_U_stT(
    n_bootstrap = 100,
    q_t = q_t,
    data = df,
    tao = 120,
    ps_formula = ps_formula,
    Q_T_formula = Q_T_formula,
    Q_cure_formula = Q_cure_formula,
    Q_T_length = Q_T_length,
    Tx = Tx,
    stp_sz = 0.001,
    tol = 1e-3,
    max_iter = 10000,
    range = c(-10, 10),
    num_starts = 1
  )
}

# 停止并行环境
stopCluster(cl)

# 查看结果
names(results) <- as.character(q_t_list) # 给结果列表添加 q_t 的键值
print(results)
saveRDS(results, "bootsrap.RData")


q_t_list = seq(54,60,1); results = list()
for (q_t in q_t_list) {
  print(paste("Processing q_t =", q_t))
  
  # 调用 solve_penalty_1 进行估计
  result <- AIPW_U_stT(n_bootstrap = 100, q_t = q_t, data = df, tao = 120, ps_formula, Q_T_formula, 
                       Q_cure_formula, Q_T_length, Tx, stp_sz = 0.001, tol = 1e-3, 
                       max_iter = 10000, range = c(-10, 10), num_starts = 1)
  
  # 将结果存入列表，以 q_t 为键
  results[[as.character(q_t)]] <- result
}

results_matrix <- do.call(rbind, lapply(names(results), function(q_t) {
  cbind(
    q_t = as.numeric(q_t),  # 添加 q_t 列
    t(results[[q_t]]$beta_hat),  # beta_hat 值（假设为向量）
    EU_hat = results[[q_t]]$EU_hat, 
    ET_hat = results[[q_t]]$ET_hat, 
    t(results[[q_t]]$beta_hat_se),  # beta_hat 的标准误差
    EU_se = results[[q_t]]$EU_se, 
    ET_se = results[[q_t]]$ET_se
  )
}))

# 为矩阵设置列名
colnames(results_matrix) <- c(
  "q_t",
  paste0("beta_hat_", seq_along(results[[1]]$beta_hat)), 
  "EU_hat",
  "ET_hat",
  paste0("beta_hat_se_", seq_along(results[[1]]$beta_hat)),
  "EU_se",
  "ET_se"
)

alpha <- 0.1
z <- qnorm(1 - alpha / 2)  # 1.96 for 95% CI

# 计算 EU 和 ET 的置信区间
results_matrix <- cbind(
  results_matrix,
  EU_CI_lower = results_matrix[, "EU_hat"] - z * results_matrix[, "EU_se"],
  EU_CI_upper = results_matrix[, "EU_hat"] + z * results_matrix[, "EU_se"],
  ET_CI_lower = results_matrix[, "ET_hat"] - z * results_matrix[, "ET_se"],
  ET_CI_upper = results_matrix[, "ET_hat"] + z * results_matrix[, "ET_se"]
)

# 查看结果
print(results_matrix)

write.csv(results_matrix, file = "result/results_matrix_maxU_stT.csv", row.names = FALSE)
results_matrix <- read.csv("result/results_matrix_maxU_stT.csv")
############################### figure: d-2 ###############################
results_matrix <- as.data.frame(results_matrix)
library(ggplot2)

# 确保 df1 已正确创建
df1 <- NULL
df1$U <- results_matrix$EU_hat
df1$low_U <- results_matrix$EU_CI_lower
df1$up_U <- results_matrix$EU_CI_upper
df1$T <- results_matrix$ET_hat
df1$low_T <- results_matrix$ET_CI_lower
df1$up_T <- results_matrix$ET_CI_upper
df1$x <- c(54, 55, 56, 57, 58, 59, 60)
df1 <- as.data.frame(df1)
row.names(df1) <- NULL

# 计算 U 和 T 的范围
U_range <- range(df1$U)
T_range <- range(df1$T)

# 设置 y 轴的范围
ylim.T <- c(45, 70)  # 左轴范围
ylim.U <- c(0.2, 0.7)  # 右轴范围

# 确定变换参数
b <- diff(ylim.T) / diff(ylim.U)
a <- ylim.T[1] - b * ylim.U[1]

# 开始绘图
ggplot(df1, aes(x = x)) +
  theme_bw() +  # 背景
  # 绘制 T 和 U 的线条
  geom_line(aes(y = T, color = 'T'), size = 1.5) +
  geom_line(aes(y = a + U * b, color = 'U'), size = 1.5, linetype = 'dashed') +
  # 添加 U 的置信区间
  geom_ribbon(aes(ymin = a + low_U * b, ymax = a + up_U * b), fill = "blue", alpha = 0.2) +
  # 颜色图例和轴标签
  scale_color_discrete(labels = c(
    expression(tilde(M)[aipw](hat(g)[alpha])), 
    expression(tilde(V)[aipw](hat(g)[alpha]))
  )) +
  scale_x_continuous(expression(q[t]), breaks = seq(54, 61, by = 1)) +
  scale_y_continuous(
    name = expression(tilde(M)[aipw](hat(g)[alpha])),  # 左轴标签
    limits = ylim.T,
    sec.axis = sec_axis(~ (.-a)/b, name = expression(tilde(V)[aipw](hat(g)[alpha])))  # 右轴标签
  ) +
  # 设置轴标签样式
  theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5),
    axis.title.y.right = element_text(angle = 0, vjust = 0.5)
  ) +
  # 添加点和标记
  geom_point(aes(y = T, color = "T"), size = 4.5, shape = 24) +
  geom_point(aes(y = a + U * b, color = "U"), size = 4.5, shape = 21) +
  geom_hline(yintercept = 39.691, linetype = "dotted", size = 1.5) +
  geom_hline(yintercept = a + 0.404 * b, linetype = "dotdash", size = 1.5) +
  annotate("text", x = 60, y = 38, label = expression(tilde(M)[aipw](g^RCT)), hjust = 1) +
  annotate("text", x = 60, y = a + 0.393 * b, label = expression(tilde(V)[aipw](g^RCT)), hjust = 1) +
  # 图例设置
  theme(legend.position = "top", legend.title = element_blank()) +
  guides(colour = guide_legend(override.aes = list(shape = c(24, 21))))

library(ggplot2)

ggplot(df1, aes(x = x)) +
  theme_bw() +  # 背景
  # 绘制 T 和 U 的线条
  geom_line(aes(y = T, color = 'T'), size = 1.5) +
  geom_line(aes(y = a + U * b, color = 'U'), size = 1.5, linetype = 'dashed') +
  # 使用 geom_errorbar 添加 U 的置信区间
  geom_errorbar(
    aes(ymin = a + low_U * b, ymax = a + up_U * b, color = 'U'),
    width = 0.2, size = 1
  ) +
  # 颜色图例和轴标签
  scale_color_discrete(labels = c(
    expression(tilde(M)[aipw](hat(g)[alpha])), 
    expression(tilde(V)[aipw](hat(g)[alpha]))
  )) +
  scale_x_continuous(expression(q[t]), breaks = seq(54, 61, by = 1)) +
  scale_y_continuous(
    name = expression(tilde(M)[aipw](hat(g)[alpha])),  # 左轴标签
    limits = ylim.T,
    sec.axis = sec_axis(~ (.-a)/b, name = expression(tilde(V)[aipw](hat(g)[alpha])))  # 右轴标签
  ) +
  # 设置轴标签样式
  theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5),
    axis.title.y.right = element_text(angle = 0, vjust = 0.5)
  ) +
  # 添加点和标记
  geom_point(aes(y = T, color = "T"), size = 4.5, shape = 24) +
  geom_point(aes(y = a + U * b, color = "U"), size = 4.5, shape = 21) +
  geom_hline(yintercept = 39.691, linetype = "dotted", size = 1.5) +
  geom_hline(yintercept = a + 0.404 * b, linetype = "dotdash", size = 1.5) +
  annotate("text", x = 60, y = 38, label = expression(tilde(M)[aipw](g^RCT)), hjust = 1) +
  annotate("text", x = 60, y = a + 0.393 * b, label = expression(tilde(V)[aipw](g^RCT)), hjust = 1) +
  # 图例设置
  theme(legend.position = "top", legend.title = element_blank()) +
  guides(colour = guide_legend(override.aes = list(shape = c(24, 21))))
