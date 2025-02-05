rm(list=ls())
library(survival);library(rgenoud);library(gfcure);library(dplyr)
set.seed(2023)
df <- read.csv("seer_esophageal.csv")
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

df$age <- (df$Age-min(df$Age))/(max(df$Age) - min(df$Age))
# remove the patients with 0 survival time
df = df[df$Survival.months!=0,]

write.csv(df,"seer_esophageal_clean.csv",row.names = F)
