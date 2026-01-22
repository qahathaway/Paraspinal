
####Compare 2 ROC-AUC####
library(pROC)
sessionInfo()

AUC1 <- read.csv("/Users/quincy/Documents/Research/Hopkins/Paraspine/ANALYSIS/Radiomics/Exam5_Para/Diabetes/New_Diabetes_Radiomics.csv", header=TRUE)

# Basic example with 2 roc objects
roc1 <- roc(AUC1$Change, AUC1$Baseline)
roc2 <- roc(AUC1$Change, AUC1$Radiomics)
roc3 <- roc(AUC1$Pred, AUC1$IMAT2)
roc4 <- roc(AUC1$Pred, AUC1$All2)

plot(roc1, type="l", col="lightsteelblue1", lwd=5, legacy.axes = TRUE, xlim=c(0.4,0))
par(new=TRUE)
plot(roc2, type="l", col="lightsteelblue3", lwd=5, legacy.axes = TRUE, xlim=c(0.4,0))
par(new=TRUE)
plot(roc3, type="l", col="slategray4", lwd=5, legacy.axes = TRUE, xlim=c(0.4,0))
par(new=TRUE)
plot(roc4, type="l", col="steelblue4", lwd=5, legacy.axes = TRUE, xlim=c(0.4,0))

roc1

library(caret)

###FRAXnb
AUC1[,'Pred']<-factor(AUC1[,'Pred'])
AUC1$FRAXnb2 <- ifelse(AUC1$FRAXnb2 > 0.5,1,0)
AUC1[,'FRAXnb2']<-factor(AUC1[,'FRAXnb2'])
precision <- posPredValue(AUC1$FRAXnb2, AUC1$Pred)
recall <- sensitivity(AUC1$FRAXnb2, AUC1$Pred)
F1 <- (2 * precision * recall) / (precision + recall)
confusionMatrix(AUC1$FRAXnb2, AUC1$Pred, positive ="1")

###Ex
AUC1[,'Pred']<-factor(AUC1[,'Pred'])
AUC1$Ex2 <- ifelse(AUC1$Ex2 > 0.5,1,0)
AUC1[,'Ex2']<-factor(AUC1[,'Ex2'])
precision <- posPredValue(AUC1$Ex2, AUC1$Pred)
recall <- sensitivity(AUC1$Ex2, AUC1$Pred)
F1 <- (2 * precision * recall) / (precision + recall)
confusionMatrix(AUC1$Ex2, AUC1$Pred, positive ="1")

###Volume
AUC1[,'Pred']<-factor(AUC1[,'Pred'])
AUC1$Volume2 <- ifelse(AUC1$Volume2 > 0.5,1,0)
AUC1[,'Volume2']<-factor(AUC1[,'Volume2'])
precision <- posPredValue(AUC1$Volume2, AUC1$Pred)
recall <- sensitivity(AUC1$Volume2, AUC1$Pred)
F1 <- (2 * precision * recall) / (precision + recall)
confusionMatrix(AUC1$Volume2, AUC1$Pred, positive ="1")

###Attenuation
AUC1[,'Pred']<-factor(AUC1[,'Pred'])
AUC1$Median <- ifelse(AUC1$Median > 0.5,1,0)
AUC1[,'Median']<-factor(AUC1[,'Median'])
precision <- posPredValue(AUC1$Median, AUC1$Pred)
recall <- sensitivity(AUC1$Median, AUC1$Pred)
F1 <- (2 * precision * recall) / (precision + recall)
confusionMatrix(AUC1$Median, AUC1$Pred, positive ="1")

roc2
ci(roc2)
roc.test(roc2, roc3, method="delong")
roc.test(roc1, roc2, method="bootstrap", boot.n=10000)
roc.test(roc1, roc2, method="venkatraman", boot.n=10000)

roc.test(roc1, roc2, method="specificity", specificity=0.9)
roc.test(roc1, roc2, method="sensitivity", sensitivity=0.2)


# Install and load the sensemakr package
library(sensemakr)
library(mlr)

my_data <- read.csv("/Users/quincy/Documents/Research/Hopkins/Paraspine/ANALYSIS/Radiomics/Exam5_Para/Paraspinal_Paper/Paraspinal_Manuscript_Sensitivity_Analysis.csv", header=TRUE)

imputed = impute(my_data, target = character(0), classes = list(numeric = imputeMedian(), integer = imputeMedian()))
my_data <- as.data.frame(imputed$data)

model1 <- lm(vBMD ~Age + Gender + Race + Site + BMI + Steroids + Bisphosphonates +
               Alcohol + Alcohol_Consumption + Tobacco + Pack_Years, data = my_data)

model2 <- lm(vBMD ~Age + Gender + Race + Site + BMI + HTN + Diabetes +
               Liver + Kidney + Creatinine + Cancer + Steroids + Bisphosphonates + Warfarin +
               Thyroid + Alcohol + Alcohol_Consumption +
               Tobacco + Pack_Years + Moderate_Phys + Vigorous_Phys + AGATSTON, data = my_data)

model3 <- lm(vBMD ~Age + Gender + Race + Site + BMI + HTN + Diabetes +
               Liver + Kidney + Creatinine + Cancer + Steroids + Bisphosphonates + Warfarin +
               Thyroid + Alcohol + Alcohol_Consumption +
               Tobacco + Pack_Years +
               Moderate_Phys + Vigorous_Phys + AGATSTON + Shape, data = my_data)

model4 <- lm(vBMD ~Age + Gender + Race + Site + BMI + HTN + Diabetes +
               Liver + Kidney + Creatinine + Cancer + Steroids + Bisphosphonates + Warfarin +
               Thyroid + Alcohol + Alcohol_Consumption +
               Tobacco + Pack_Years +
               Moderate_Phys + Vigorous_Phys + AGATSTON + Median_New, data = my_data)

model5 <- lm(vBMD ~Age + Gender + Race + Site + BMI + HTN + Diabetes +
               Liver + Kidney + Creatinine + Cancer + Steroids + Bisphosphonates + Warfarin +
               Thyroid + Alcohol + Alcohol_Consumption +
               Tobacco + Pack_Years +
               Moderate_Phys + Vigorous_Phys + AGATSTON + Shape + Median_New, data = my_data)

model6 <- lm(vBMD ~Median_New + Gender + Race + Site + BMI + HTN + Diabetes +
               Liver + Kidney + Creatinine + Cancer + Steroids + Bisphosphonates + Warfarin +
               Thyroid + Alcohol + Alcohol_Consumption +
               Tobacco + Pack_Years + Manufacturer + Model + Age +
               Moderate_Phys + Vigorous_Phys + AGATSTON + Shape, data = my_data)

anova(model4, model5)
AIC(model1, model2)
BIC(model1, model2)
summary(model2)$r.squared
summary(model6)$adj.r.squared
summary(model6)
plot(model5)

residuals <- residuals(model5)
predictions <- predict(model5)
actual_values <- my_data$vBMD
mse <- mean(residuals^2)
rmse <- sqrt(mse)
mae <- mean(abs(residuals))
print(paste("MSE:", mse))
print(paste("RMSE:", rmse))
print(paste("MAE:", mae))
summary(model5)$adj.r.squared


sensitivity <- sensemakr(model = model5, 
                         treatment = "Median_New",
                         benchmark_covariates = "Bisphosphonates",
                         kd = 1:3)

ovb_minimal_reporting(sensitivity, format = "html")

# Print and plot a summary of the results
summary(sensitivity)
plot(sensitivity)
plot(sensitivity, sensitivity.of = "t-value")
plot(sensitivity, type = "extreme")


library(randtests)
library(mdatools)

my_data <- read.csv("/Users/quincy/Documents/Research/Hopkins/Paraspine/ANALYSIS/Radiomics/Exam5_Para/Paraspinal_Paper/Paraspinal_Manuscript_Randomness.csv", header=TRUE)
my_data <- read.csv("/Users/quincy/Documents/Research/Hopkins/Paraspine/ANALYSIS/Radiomics/Exam5_Para/Paraspinal_Paper/Paraspinal_Manuscript_Randomness_2.csv", header=TRUE)

Random <- runs.test(my_data$original_firstorder_Median, plot = "TRUE")

library(sas7bdat)
x<-read.sas7bdat("/Users/quincy/Documents/Research/Hopkins/Paraspine/ANALYSIS/Radiomics/MESAe16_BisphosMeds_20220502.sas7bdat", debug=TRUE)
write.csv(x,"/Users/quincy/Documents/Research/Hopkins/Paraspine/ANALYSIS/Radiomics/MESAe16_BisphosMeds_20220502.csv", row.names = FALSE) 


####Compare MultipvBMD####Compare Multiple Classes ROC-AUC####
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