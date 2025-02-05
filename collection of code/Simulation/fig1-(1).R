library(mixcure);library(rgenoud);library(gfcure)
set.seed(2023)

# maxT, st. U
q_u = c(1e-5,0.36,0.37, 0.38,0.39,0.40,0.41,0.42,0.43,0.44,0.45,0.46);i = 1 # 12个约束点
result_maxT_stU <- matrix(NA, nrow=length(q_u), ncol = 4)
for (q in q_u){
    print(q_u[i])
    
    registerDoRNG(2023);a0 = proc.time()[3]
    result_T_ipw_p1 = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
        MaxT_ipw_p1(q_cons=q, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 1)};print(proc.time()[3] - a0)
    result_maxT_stU[i, ] = colMeans(break_result(result_T_ipw_p1)[[2]])
    
    i = i+1
}
result_maxT_stU = result_maxT_stU[sort(result_maxT_stU[,1], index.return=TRUE)$ix,]

x = result_maxT_stU[,2]; y = result_maxT_stU[,4]
plot(x, y,col = "darkred",xlab="",ylab="",cex=1.5,ylim = c(1.8, 3.5))
text(x, y, labels = 1:length(x), pos = 1, cex = 1, col = "black")
# plot(y1, x1, col = "darkred", lwd = 2,pch=c(1,2,3,4,5) ,xlim =c(0.35, 0.46) , ylim = c(1.8, 3.5), xlab="",ylab="",cex=2)
# Uncured mean survival time versus Cure rate
title(main = "", xlab = expression(V(g[alpha])), ylab = expression(M(g[alpha])))
lines(x, y, col = "dodgerblue3", lwd = 2, lty=4, pch=3) # 0.375
abline(v = y[1], col = "black", lwd = 2, lty=3)
      # 水平虚线
abline(h = x[[5]], col = "black", lwd = 2, lty=3)
      # legend with math expression
legend("topright", c("t-max", expression(q[u]==0.375), expression(q[u]==0.4), expression(q[u]==0.425), "u-max"), col = "darkred", pch=c(1,2,3,4,5))
      # legend("topright", c("t-max", "q=3", "q=4.5", "q=6","u-max "), col = "darkred", pch=c(1,2,3,4,5), lwd = 2)
      
saveRDS(result_maxT_stU, "D:\\论文\\Cure Model + ITR\\模拟\\复现整理版 - 副本\\result_maxT_stU.RDS")
saveRDS(result_maxU_stT, "D:\\论文\\Cure Model + ITR\\模拟\\复现整理版 - 副本\\result_maxU_stT.RDS")
# maxU, st. T
q_t = c(1e-5, 1.9, 2.0, 2.1, 2.2, 2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2);i = 1 # 15个约束点
result_maxU_stT <- matrix(NA, nrow=length(q_t), ncol = 4)
for (q in q_t){
    print(q_t[i])
    
    registerDoRNG(2023);a0 = proc.time()[3]
    result_U_ipw_p2 = foreach(i=1:10, .combine='rbind', .multicombine = TRUE) %dopar% {
        MaxU_ipw_p2(q_cons=q, n=1000, gamma_p, gamma_t, lambda_c, sim_num = 1)};print(proc.time()[3] - a0)
    result_maxU_stT[i, ] = colMeans(break_result(result_U_ipw_p2)[[2]])
    
    i = i+1
}
result_maxU_stT
y = result_maxU_stT[,2]; x = result_maxU_stT[,4]
plot(x, y,col = "darkred",xlab="",ylab="",cex=1.5, ylim=c(0.34,0.44))
text(x, y, labels = 1:length(x), pos = 1, cex = 1, col = "black")
title(main = "", xlab = expression(M(g[alpha])), ylab = expression(V(g[alpha])))
lines(x,y, col = "dodgerblue3", lwd = 2, lty=4, pch=3) # 0.375
# 画三条垂直虚线
# abline(v = 3, col = "black", lwd = 2, lty=3)
# abline(v = 4.5, col = "black", lwd = 2, lty=3)
# abline(v = 6, col = "black", lwd = 2, lty=3)
abline(v = y[[1]], col = "black", lwd = 2, lty=3)
# 水平虚线
abline(h = x[[5]], col = "black", lwd = 2, lty=3)
# legend with math expression
legend("topright", c("t-max", expression(q[u]==0.375), expression(q[u]==0.4), expression(q[u]==0.425), "u-max"), col = "darkred", pch=c(1,2,3,4,5))
# legend("topright", c("t-max", "q=3", "q=4.5", "q=6","u-max "), col = "darkred", pch=c(1,2,3,4,5), lwd = 2)

###################################################################
library(ggplot2)

# 创建数据框
data <- data.frame(
    Cure_Rate = unlist(y),
    Mean_Survival_Time = unlist(x),
    Group = c("t-max", "q[u] == 0.375", "q[u] == 0.4", "q[u] == 0.425", "u-max")
)
ggplot(data, aes(x = Cure_Rate, y = Mean_Survival_Time, color = Group)) +
    geom_point(shape = c(1, 2, 3, 4, 5), size = 3) +  # 设置不同点的形状和大小
    geom_line(linetype = "dashed", size = 1) # 添加折线，并设置线型和大小
# 创建ggplot图
ggplot(data, aes(x = Cure_Rate, y = Mean_Survival_Time, color = Group)) +
    geom_point(shape = c(1, 2, 3, 4, 5), size = 3) +  # 设置不同点的形状和大小
    geom_line(linetype = "dashed", size = 1) +  # 添加折线，并设置线型和大小
    theme_minimal() +  # 使用简洁主题
    labs(
        title = "Uncured Mean Survival Time versus Cure Rate",
        x = expression(V(hat(eta))),
        y = expression(M(hat(eta)))
    ) +  # 设置标题和轴标签
    scale_color_manual(
        values = c("darkred", "dodgerblue3", "dodgerblue3", "dodgerblue3", "darkred"),
        breaks = c("t-max", "q[u] == 0.375", "q[u] == 0.4", "q[u] == 0.425", "u-max"),
        labels = c("t-max", expression(q[u] == 0.375), expression(q[u] == 0.4), expression(q[u] == 0.425), "u-max")
    ) +  # 手动设置颜色和图例标签
    geom_vline(xintercept = y[1], linetype = "dashed", color = "black") +  # 添加垂直虚线
    geom_hline(yintercept = x[5], linetype = "dashed", color = "black") +  # 添加水平虚线
    scale_x_continuous(limits = c(0.35, 0.46)) +  # 设置x轴的限制
    scale_y_continuous(limits = c(2, 7.5)) +  # 设置y轴的限制
    theme(legend.position = "top", legend.title = element_blank())  # 设置图例位置和样式


###################################################################
# plot(density(T_list[[4]]), col = "darkred", lwd = 2, xlim = c(0, 8), ylim = c(0, 0.5), main = "Density of uncured survival time", xlab = "T", ylab = "")
# # 体现趋势
# lines(density(T_list[[1]]), col = "black", lwd = 2, lty=4, pch=3) # 0.375
# lines(density(T_list[[2]]), col = "ivory4", lwd = 2, lty=5)
# lines(density(T_list[[3]]), col = "dodgerblue3", lwd = 2, lty=6)
# # legend
# legend("topright", c("q=0", "p = 0.375", "p = 0.4", "p = 0.425"), col = c("darkred", "black", "ivory4", "dodgerblue3"), lty = c(1,4,5,6), lwd = 2)
# ################################################################### 

# sum(T_list[[2]]==0)

# mean(T_list[[1]]); sd(T_list[[1]])
# mean(T_list[[2]]); sd(T_list[[2]])
# mean(T_list[[3]]); sd(T_list[[3]])
# mean(T_list[[4]]); sd(T_list[[4]])
# mean(T_list[[5]]); sd(T_list[[5]])

# t1 <- as.matrix(B0[,-1])
# plot(t1[6, ]/t1[5,], type = 'b', lty = 1, lwd = 2, pch = 15, cex.axis = 0.7, ylim = c(0.8,1.5),
#      xaxt = 'n', ylab = 'Relative Squared Risk', xlab = 'p', cex.lab = 1.1, col = 'darkred',
#      main = substitute(paste('K = ', B), list(B=B)))
# xy <- par('usr')
# axis(1, at = c(1:5), labels = c(10,50,100,200,300), cex.axis = 0.7)
# #lines(t1[1, ]/t1[5,], type = 'b', lty = 5, pch = 6, lwd = 2, cex=1, col = 'darkred')
# lines(t1[2, ]/t1[5,], type = 'b', lty = 5, pch = 1, lwd = 2, cex=1, col = 'burlywood3')
# lines(t1[3, ]/t1[5,], type = 'b', lty = 4, pch = 17, lwd = 2, col = 'ivory4')
# #lines(t3[4, ], type = 'b', lty = 1, pch = 16, lwd = 1, col = '#91B493')
# lines(t1[4, ]/t1[5,], type = 'b', lty = 1, pch = 16, lwd = 1.5, col = '#91B493')
# lines(t1[5, ]/t1[5,], type = 'l', lty = 1, pch = 16, lwd = 1.5, col= 'black')
# lines(t1[1, ]/t1[5,], type = 'b', lty = 5, pch = 3, lwd = 2, cex=1, col = 'dodgerblue3')