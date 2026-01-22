
####Compare 2 ROC-AUC####
library(pROC)
sessionInfo()

AUC1 <- read.csv("/Users/quincy/Documents/Research/Hopkins/Paraspine/ANALYSIS/Radiomics/Exam5_Para/Paraspinal_Paper/R/Predictions_ALL.csv", header=TRUE)

# Basic example with 2 roc objects
roc1 <- roc(AUC1$Pred, AUC1$FRAXnb)
roc2 <- roc(AUC1$Pred, AUC1$Ex)
roc3 <- roc(AUC1$Pred, AUC1$Volume)
roc4 <- roc(AUC1$Pred, AUC1$Median)

plot(roc1, type="l", col="black", lwd=5, legacy.axes = TRUE, xlim=c(0.4,0))
par(new=TRUE)
plot(roc2, type="l", col="grey", lwd=5, legacy.axes = TRUE, xlim=c(0.4,0))
par(new=TRUE)
plot(roc3, type="l", col="steelblue1", lwd=5, legacy.axes = TRUE, xlim=c(0.4,0))
par(new=TRUE)
plot(roc4, type="l", col="steelblue4", lwd=5, legacy.axes = TRUE, xlim=c(0.4,0))

roc4
ci(roc4)
roc.test(roc2, roc4, method="delong")
roc.test(roc1, roc2, method="bootstrap", boot.n=10000)
roc.test(roc1, roc2, method="venkatraman", boot.n=10000)

roc.test(roc1, roc2, method="specificity", specificity=0.9)
roc.test(roc1, roc2, method="sensitivity", sensitivity=0.2)


AUC1 <- read.csv("/Users/quincy/Documents/Research/Hopkins/Paraspine/ANALYSIS/Radiomics/Exam5_Para/Paraspinal_Paper/R/Predictions_ALL_2.csv", header=TRUE)

# Basic example with 2 roc objects
roc1 <- roc(AUC1$Pred, AUC1$FRAXnb)
roc2 <- roc(AUC1$Pred, AUC1$Bone)
roc3 <- roc(AUC1$Pred, AUC1$V)
roc4 <- roc(AUC1$Pred, AUC1$M)
roc5 <- roc(AUC1$Pred, AUC1$BV)
roc6 <- roc(AUC1$Pred, AUC1$BM)

plot(roc1, type="l", col="black", lwd=5, legacy.axes = TRUE, xlim=c(0.4,0))
par(new=TRUE)
plot(roc2, type="l", col="grey", lwd=5, legacy.axes = TRUE, xlim=c(0.4,0))
par(new=TRUE)
plot(roc3, type="l", col="palegreen1", lwd=5, legacy.axes = TRUE, xlim=c(0.4,0))
par(new=TRUE)
plot(roc4, type="l", col="steelblue1", lwd=5, legacy.axes = TRUE, xlim=c(0.4,0))
par(new=TRUE)
plot(roc5, type="l", col="palegreen4", lwd=5, legacy.axes = TRUE, xlim=c(0.4,0))
par(new=TRUE)
plot(roc6, type="l", col="steelblue4", lwd=5, legacy.axes = TRUE, xlim=c(0.4,0))

roc5
ci(roc5)
roc.test(roc2, roc6, method="delong")
roc.test(roc1, roc2, method="bootstrap", boot.n=10000)
roc.test(roc1, roc2, method="venkatraman", boot.n=10000)

roc.test(roc1, roc2, method="specificity", specificity=0.9)
roc.test(roc1, roc2, method="sensitivity", sensitivity=0.2)

library(caret)

###FRAXnb
AUC1[,'Pred']<-factor(AUC1[,'Pred'])
AUC1$FRAXnb <- ifelse(AUC1$FRAXnb > 0.5,1,0)
AUC1[,'FRAXnb']<-factor(AUC1[,'FRAXnb'])
precision <- posPredValue(AUC1$FRAXnb, AUC1$Pred)
recall <- sensitivity(AUC1$FRAXnb, AUC1$Pred)
F1 <- (2 * precision * recall) / (precision + recall)
confusionMatrix(AUC1$FRAXnb, AUC1$Pred, positive ="1")

###Ex
AUC1[,'Pred']<-factor(AUC1[,'Pred'])
AUC1$Ex <- ifelse(AUC1$Ex > 0.5,1,0)
AUC1[,'Ex']<-factor(AUC1[,'Ex'])
precision <- posPredValue(AUC1$Ex, AUC1$Pred)
recall <- sensitivity(AUC1$Ex, AUC1$Pred)
F1 <- (2 * precision * recall) / (precision + recall)
confusionMatrix(AUC1$Ex, AUC1$Pred, positive ="1")

###Volume
AUC1[,'Pred']<-factor(AUC1[,'Pred'])
AUC1$Volume <- ifelse(AUC1$Volume > 0.5,1,0)
AUC1[,'Volume']<-factor(AUC1[,'Volume'])
precision <- posPredValue(AUC1$Volume, AUC1$Pred)
recall <- sensitivity(AUC1$Volume, AUC1$Pred)
F1 <- (2 * precision * recall) / (precision + recall)
confusionMatrix(AUC1$Volume, AUC1$Pred, positive ="1")

###Attenuation
AUC1[,'Pred']<-factor(AUC1[,'Pred'])
AUC1$Median <- ifelse(AUC1$Median > 0.5,1,0)
AUC1[,'Median']<-factor(AUC1[,'Median'])
precision <- posPredValue(AUC1$Median, AUC1$Pred)
recall <- sensitivity(AUC1$Median, AUC1$Pred)
F1 <- (2 * precision * recall) / (precision + recall)
confusionMatrix(AUC1$Median, AUC1$Pred, positive ="1")


####Compare Multiple Classes ROC-AUC####
library(multiROC)

Multi <- data.frame(read.csv(
  file = '/Users/quincy/Documents/Research/HVI/MESA/FINAL/Pred/Multi_M.csv'))

res <- multi_roc(Multi, force_diag=T)
unlist(res$AUC)

multi_roc_auc <- function(true_pred_data, idx) {
  results <- multi_roc(true_pred_data[idx, ])$AUC
  results <- unlist(results)
  return(results)
}

roc_auc_with_ci_res <- roc_auc_with_ci(Multi, conf= 0.95, type='basic', R = 1000)
roc_auc_with_ci_res

#roc_test <- multi_roc(Multi)
#roc_test$Sensitivity
#roc_test$Specificity

n_method <- length(unique(res$Methods))
n_group <- length(unique(res$Groups))
res_df <- data.frame(Specificity= numeric(0), Sensitivity= numeric(0), Group = character(0), AUC = numeric(0), Method = character(0))
for (i in 1:n_method) {
  for (j in 1:n_group) {
    temp_data_1 <- data.frame(Specificity=res$Specificity[[i]][j],
                              Sensitivity=res$Sensitivity[[i]][j],
                              Group=unique(res$Groups)[j],
                              AUC=res$AUC[[i]][j],
                              Method = unique(res$Methods)[i])
    colnames(temp_data_1) <- c("Specificity", "Sensitivity", "Group", "AUC", "Method")
    res_df <- rbind(res_df, temp_data_1)
    
  }
}

ggplot2::ggplot(res_df, ggplot2::aes(x = 1-Specificity, y=Sensitivity)) + ggplot2::geom_path(ggplot2::aes(color = Group, linetype=Method)) + ggplot2::geom_segment(ggplot2::aes(x = 0, y = 0, xend = 1, yend = 1),
      colour='grey', linetype = 'dotdash') + ggplot2::theme_bw() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.justification=c(1, 0),
      legend.position=c(.95, .05), legend.title=ggplot2::element_blank(),
      legend.background = ggplot2::element_rect(fill=NULL, size=0.5,
      linetype="solid", colour ="black"))


####Using Macro and Micro####
#for (i in 1:n_method) {
#  for (j in 1:n_group) {
#    temp_data_1 <- data.frame(Specificity=res$Specificity[[i]][j],
#                              Sensitivity=res$Sensitivity[[i]][j],
#                              Group=unique(res$Groups)[j],
#                              AUC=res$AUC[[i]][j],
#                              Method = unique(res$Methods)[i])
#    colnames(temp_data_1) <- c("Specificity", "Sensitivity", "Group", "AUC", "Method")
#    res_df <- rbind(res_df, temp_data_1)
#    
#  }
#  temp_data_2 <- data.frame(Specificity=res$Specificity[[i]][n_group+1],
#                            Sensitivity=res$Sensitivity[[i]][n_group+1],
#                            Group= "Macro",
#                            AUC=res$AUC[[i]][n_group+1],
#                            Method = unique(res$Methods)[i])
#  temp_data_3 <- data.frame(Specificity=res$Specificity[[i]][n_group+2],
#                            Sensitivity=res$Sensitivity[[i]][n_group+2],
#                            Group= "Micro",
#                            AUC=res$AUC[[i]][n_group+2],
#                            Method = unique(res$Methods)[i])
#  colnames(temp_data_2) <- c("Specificity", "Sensitivity", "Group", "AUC", "Method")
#  colnames(temp_data_3) <- c("Specificity", "Sensitivity", "Group", "AUC", "Method")
#  res_df <- rbind(res_df, temp_data_2)
#  res_df <- rbind(res_df, temp_data_3)
#}