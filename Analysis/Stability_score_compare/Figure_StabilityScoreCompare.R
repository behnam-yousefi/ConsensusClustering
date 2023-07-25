## Figure 3

rm(list = ls())
setwd("~/Desktop/R_Root/ConsensusClustering")
library(ggplot2)


dir = "Data/results/All_g_pam/"
set = list.files(dir)

names = c("PCC_PAM_RS", "PCC_PAM_LS", "PCC_PAM_PAC", "PCC_PAM_DA", "PCC_PAM_CMavg", 
          "nullPCC_PAM_RS", "nullPCC_PAM_LS", "nullPCC_PAM_PAC", "nullPCC_PAM_DA", "nullPCC_PAM_CMavg",
          "M3C_PAM_Ent_RCSI", "M3C_PAM_Ent_Pval", "Sil", "K_real", "sd_max", "rep")

results = data.frame()
is_correct = data.frame()
for (i in set){
  x = readRDS(paste0(dir,i))
  x = x[,names]
  results = rbind(results, x)
  
  y = cbind(x[,14:16],
            x[,1:13] - x[,14] == 0)
  is_correct = rbind(is_correct, y)
}
dim(results)
dim(is_correct)

# Figure 2a-e
k = 5
methods = c("Sil", "PCC_PAM_CMavg", "PCC_PAM_DA", "M3C_PAM_Ent_Pval", 
            "M3C_PAM_Ent_RCSI", "PCC_PAM_PAC", "PCC_PAM_LS")
ACC = data.frame()
for (s in 1:10){
  x = is_correct[is_correct$K_real == k & is_correct$sd_max == s, ]
  counts = sort(apply(x[4:16], 2, sum), decreasing = TRUE)
  acc = counts[methods] / nrow(x)
  acc = data.frame(acc, method = methods, sd = s/10)
  ACC = rbind(ACC, acc)
}

ACC_smoothed = ACC
for (m in unique(ACC$method)){
  for (sd in (2:9)/10){
    ACC_smoothed[ACC_smoothed$method == m & ACC_smoothed$sd == sd, "acc"] = mean(c(
      ACC[ACC$method == m & ACC$sd == (sd), "acc"],
      ACC[ACC$method == m & ACC$sd == as.character(sd + .1), "acc"],
      ACC[ACC$method == m & ACC$sd == as.character(sd - .1), "acc"]))
  }
}

ACC_smoothed$method = factor(ACC_smoothed$method, 
                             levels = ACC_smoothed$method[1:7])

col = c("orange", "yellow", "purple", "lightblue", "blue", "green", "red")
plt = ggplot(ACC_smoothed) +
  geom_line(aes(x = sd, y = acc, color = method), size = .8) + 
  geom_point(aes(x = sd, y = acc, color = method), size = 1.5) + 
  ylim(c(0,1)) + 
  scale_color_manual(values=col) + 
  theme_classic(base_size = 10) + theme(legend.position = "right")
plt

pdf("Figures/Figure_2d.pdf", height = 4, width = 5)
plt
dev.off()

# Figure 2f: Mean accuracy
methods = c("Sil", "PCC_PAM_CMavg", "PCC_PAM_DA", "M3C_PAM_Ent_Pval", 
            "M3C_PAM_Ent_RCSI", "PCC_PAM_PAC", "PCC_PAM_LS")
ACC = data.frame()
for (s in 1:10){
  x = is_correct[is_correct$sd_max == s, ]
  counts = sort(apply(x[4:16], 2, sum), decreasing = TRUE)
  acc = counts[methods] / nrow(x)
  acc = data.frame(acc, method = methods, sd = s/10)
  ACC = rbind(ACC, acc)
}

ACC_smoothed = ACC
for (m in unique(ACC$method)){
  for (sd in (2:9)/10){
    ACC_smoothed[ACC_smoothed$method == m & ACC_smoothed$sd == sd, "acc"] = mean(c(
      ACC[ACC$method == m & ACC$sd == (sd), "acc"],
      ACC[ACC$method == m & ACC$sd == as.character(sd + .1), "acc"],
      ACC[ACC$method == m & ACC$sd == as.character(sd - .1), "acc"]))
  }
}

ACC_smoothed$method = factor(ACC_smoothed$method, 
                             levels = ACC_smoothed$method[1:7])

plt = ggplot(ACC_smoothed) +
  geom_line(aes(x = sd, y = acc, color = method), size = .8) + 
  geom_point(aes(x = sd, y = acc, color = method), size = 1.5) + 
  ylim(c(0,1)) +
  scale_color_manual(values=col) + 
  theme_classic(base_size = 10) + theme(legend.position = "right")
plt

pdf("Figures/Figure_2f.pdf", height = 4, width = 5)
plt
dev.off()
