####Compare ROC-AUC####
library(pROC)
library(caret)
sessionInfo()

AUC1 <- read.csv("/path/to/predictions/for/AUC/*.csv", header=TRUE)

# Basic example with 2 roc objects
roc1 <- roc(AUC1$Pred, AUC1$FRAXnb)
roc2 <- roc(AUC1$Pred, AUC1$Expanded)
roc3 <- roc(AUC1$Pred, AUC1$Volume)
roc4 <- roc(AUC1$Pred, AUC1$Median)

plot(roc1, type="l", col="black", lwd=5, legacy.axes = TRUE, xlim=c(0.4,0))
par(new=TRUE)
plot(roc2, type="l", col="grey", lwd=5, legacy.axes = TRUE, xlim=c(0.4,0))
par(new=TRUE)
plot(roc3, type="l", col="steelblue1", lwd=5, legacy.axes = TRUE, xlim=c(0.4,0))
par(new=TRUE)
plot(roc4, type="l", col="steelblue4", lwd=5, legacy.axes = TRUE, xlim=c(0.4,0))

ci(roc1)
roc.test(roc1, roc2, method="delong")

###FRAXnb
AUC1[,'Pred']<-factor(AUC1[,'Pred'])
AUC1$FRAXnb <- ifelse(AUC1$FRAXnb > 0.5,1,0)
AUC1[,'FRAXnb']<-factor(AUC1[,'FRAXnb'])
precision <- posPredValue(AUC1$FRAXnb, AUC1$Pred)
recall <- sensitivity(AUC1$FRAXnb, AUC1$Pred)
F1 <- (2 * precision * recall) / (precision + recall)
confusionMatrix(AUC1$FRAXnb, AUC1$Pred, positive ="1")

###Expanded FRAXnb
AUC1[,'Pred']<-factor(AUC1[,'Pred'])
AUC1$Expanded <- ifelse(AUC1$Expanded > 0.5,1,0)
AUC1[,'Expanded']<-factor(AUC1[,'Expanded'])
precision <- posPredValue(AUC1$Expanded, AUC1$Pred)
recall <- sensitivity(AUC1$Expanded, AUC1$Pred)
F1 <- (2 * precision * recall) / (precision + recall)
confusionMatrix(AUC1$Expanded, AUC1$Pred, positive ="1")

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
#  colnames(temp_data_2) <- c("Specificity", "Sensitivity", "Group", "AUC", "Method")
#  colnames(temp_data_3) <- c("Specificity", "Sensitivity", "Group", "AUC", "Method")
#  res_df <- rbind(res_df, temp_data_2)
#  res_df <- rbind(res_df, temp_data_3)
#}
