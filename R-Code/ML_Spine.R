
####Regression Predictions####
library(tidyverse)
library(broom)

##Load Dataset##
MESA <- data.frame(read.csv(file = '/Users/quincy/Documents/Research/Hopkins/Paraspine/ANALYSIS/Radiomics/Exam5_Para/Paraspinal_Orange_10_R_LowBMI.csv'))

##Impute Median Values##
library(mlr)
imputed = impute(MESA, target = character(0), classes = list(numeric = imputeMedian(), integer = imputeMedian()))
final <- as.data.frame(imputed$data)

library(dplyr)
set.seed(1)
out2 <- final %>%
  group_by(Controled_BMI) %>%
  sample_n(10)

######HF Risk Prediction#####
##Training/Testing Split##
set.seed(100)
index <- createDataPartition(final$New_VCF, p = 0.6, list = FALSE)
train <- final[index, ]
test <- final[-index, ]


###Binary Classifier###
##Univariate Analysis##
#Univariate Analysis - Binary
explanatory_vars <- c(colnames(final)[8:29])

explanatory_vars %>% str_c("Low_BMI ~ ", .)

models <- explanatory_vars %>%       # begin with variables of interest
  str_c("Low_BMI ~ ", .) %>%         # combine each variable into formula ("outcome ~ variable of interest")
  
  # iterate through each univariate formula
  map(                               
    .f = ~glm(                       # pass the formulas one-by-one to glm()
      formula = as.formula(.x),      # within glm(), the string formula is .x
      family = "binomial",           # specify type of glm (logistic)
      data = final)) %>%          # dataset
  
  # tidy up each of the glm regression outputs from above
  map(
    .f = ~tidy(
      .x, 
      exponentiate = TRUE,           # exponentiate 
      conf.int = TRUE)) %>%          # return confidence intervals
  
  # collapse the list of regression outputs in to one data frame
  bind_rows() %>% 
  
  # round all numeric columns
  mutate(across(where(is.numeric), round, digits = 4))

models
summary(models)
write.csv(models, file="/Users/quincy/Documents/Research/Hopkins/Paraspine/ANALYSIS/Radiomics/Exam5_Para/Univariate_Exam5_Para.csv")

## score logistic regression models
glmDemHF = glm(Low_BMI~age5c+gender1+race1c+creatin5c+ostrd5c+warf5c+
                 thry5c+curalc5+cig5c+hrmrep1+pamcm5c+vBMD+original_firstorder_Median,data=final,family=binomial)

summary(glmDemHF)
confint(glmDemHF)

glmFuncHF <- glm(New_VCF~Age+vBMD,data=train,family=binomial)

glmDrugHF <- glm(HF~ATEZOLIZUMAB+CAPECITABINE+CARBOPLATIN+CELECOXIB+
                     CYCLOPHOSPHAMIDE+DOCETAXEL+DOXORUBICIN+ERIBULIN+
                     EVEROLIMUS+GEMCITABINE+METHOTREXATE+NERATINIB+
                     PACLITAXEL+PALBOCICLIB+PEMBROLIZUMAB+PERTUZUMAB+
                     TRASTUZUMAB,data=train,family=binomial)

glmOmicsHF <- glm(HF~IVS_PRE_CONVENTIONAL_max_ED+
                  IVS_PRE_CONVENTIONAL_std_ED+
                  IVS_PRE_CONVENTIONAL_std_ES+
                  IVS_PRE_GLRLM_GLNU_ES+
                  IVS_PRE_GLCM_Entropy_log10_ES+
                  IVS_PRE_NGLDM_Coarseness_ES+
                  PW_PRE_DISCRETIZED_HISTO_Entropy_log10_ED+
                  IVS_PRE_DISCRETIZED_std_ES+
                  PW_PRE_DISCRETIZED_HISTO_Energy.Uniformity_ED+
                  PW_PRE_DISCRETIZED_std_ED,data=train,family=binomial)

###AUC###
##Test Score AUC##
Binary_HF_Test=Score(list(glmFuncHF),formula=New_VCF~1,
         data=test,plots=c("calibration","ROC"))

##Cross-Validation Train Score AUC##
Binary_HF_Train=Score(list(glmDemHF, glmFuncHF, glmDrugHF, glmOmicsHF),formula=HF~1,
                     data=train, method="bootcv", B=10, plots=c("calibration","ROC"))

##Binary AUC##
ROC_Test <- plotROC(Binary_HF_Test, lwd=8, col=c("#4682B4"))
ROC_Train <- plotROC(Binary_HF_Train, lwd=8, col=c("#4682B4", "#AF46B4", "#B47846","#4BB446"))


######Predict Risk#####
risk_score_DemM <- predictRisk(coxDemM,times=c(3939),newdata=final)
risk_score_FuncM <- predictRisk(coxFuncM,times=c(3939),newdata=final)
risk_score_DrugM <- predictRisk(coxDrugM,times=c(3939),newdata=final)
risk_score_OmicsM <- predictRisk(coxOmicsM,times=c(3939),newdata=final)

risk_score_DemHF <- predictRisk(coxDemHF,times=c(3754),newdata=final)
risk_score_FuncHF <- predictRisk(coxFuncHF,times=c(3754),newdata=final)
risk_score_DrugHF <- predictRisk(coxDrugHF,times=c(3754),newdata=final)
risk_score_OmicsHF <- predictRisk(coxOmicsHF,times=c(3754),newdata=final)

df_list <- list(risk_score_DemM, risk_score_FuncM, risk_score_DrugM, risk_score_OmicsM)
df <- data.frame(df_list)
write.csv(df, file = "/Users/quincy/Documents/Research/Brijesh - CardioOnc/020923/NEW_CSV/RiskScore_M.csv")

df_listHF <- list(risk_score_DemHF, risk_score_FuncHF, risk_score_DrugHF, risk_score_OmicsHF)
dfHF <- data.frame(df_listHF)
write.csv(dfHF, file = "/Users/quincy/Documents/Research/Brijesh - CardioOnc/020923/NEW_CSV/RiskScore_HF.csv")




#####Hazard Ratio#####
##Load Dataset##
MESA <- data.frame(read.csv(
  file = '/Users/quincy/Documents/Research/Brijesh - CardioOnc/020923/NEW_CSV/RiskScore_HF_Test.csv'))

##Impute Median Values##
library(mlr)
imputed = impute(MESA, target = character(0), classes = list(numeric = imputeMedian(), integer = imputeMedian()))
final <- as.data.frame(imputed$data)

##Plot the baseline survival function##
fit <- surv_fit(Surv(HF_Duration, HF) ~Omics_Prob,
                data = final)

#ggsurvplot(fit, data = final, pval = TRUE, break.time.by = 500)

ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  fun="event",
  #  conf.int.style = "step",  # customize style of confidence intervals
  xlab = "Time in days",
  ylab = "Incidence",# customize X axis label.
  fontsize = 2,
  break.time.by = 500,     # break X axis in time intervals by 200.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table = "abs_pct",  # absolute number and percentage at risk.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.y.text = TRUE,
  risk.table.font = 3,
  risk.table.title = "Number at Risk (%)",
  risk.table.pos = "out",
  cumevents.title = "Cumulative Events",
  cumevents = FALSE,
  cumevents.font = 3,
  font.x = 12,
  font.y = 12,
  ncensor.plot = FALSE,      # plot the number of censored subjects at time t
  #  surv.median.line = "hv",  # add the median survival pointer.
  ylim = c(0, 1),
  legend.labs = 
    c("No HF", "HF"),    # change legend labels.
  palette = 
    c("black", "red") # custom color palettes.
)




#####NRI/IDI######
##Load Risk Score
MESA <- data.frame(read.csv(
  file = '/Users/quincy/Documents/Research/Brijesh - CardioOnc/020923/NEW_CSV/RiskScore_M_Test.csv'))

##Impute Median Values##
library(mlr)
imputed = impute(MESA, target = character(0), classes = list(numeric = imputeMedian(), integer = imputeMedian()))
final <- as.data.frame(imputed$data)

###Sample 1###
# specify column number of outcome variable
cOutcome <- 2
# specify column numbers of non-genetic predictors
cNonGenPred <- c(4)
# specify column numbers of non-genetic predictors that are categorical
cNonGenPredCat <- c(0)
# specify column numbers of genetic predictors
cGenPred <- c(0)
# specify column numbers of genetic predictors that are categorical
cGenPredCat <- c(0)
# fit logistic regression model
riskmodel1 <- fitLogRegModel(data=final, cOutcome=cOutcome,
                            cNonGenPreds=cNonGenPred, cNonGenPredsCat=cNonGenPredCat,
                            cGenPreds=cGenPred, cGenPredsCat=cGenPredCat)
summary(riskmodel1)

###Sample 2###
# specify column numbers of non-genetic predictors
cNonGenPred <- c(7)

# fit logistic regression model
riskmodel2 <- fitLogRegModel(data=final, cOutcome=cOutcome,
                            cNonGenPreds=cNonGenPred, cNonGenPredsCat=cNonGenPredCat,
                            cGenPreds=cGenPred, cGenPredsCat=cGenPredCat)
summary(riskmodel2)


# obtain multivariate OR(95% CI) for all predictors of the fitted model
ORmultivariate(riskModel=riskmodel1)
ORmultivariate(riskModel=riskmodel2)

# obtain predicted risks
predRisk1 <- predRisk(riskmodel1)
predRisk2 <- predRisk(riskmodel2)

# specify cutoff values for risk categories
cutoff <- c(0.0,0.15,1.0)
#cutoff <- c(0, 0.27, 0.36, 0.78, 1)
#compute reclassification measures
reclassification(data=final, cOutcome=cOutcome,
                 predrisk1=predRisk1, predrisk2=predRisk2, cutoff=cutoff)

# specify cutoff values for risk categories
cutoff <- c(0,.10,.20,.30,.40,.50,.60,.70,.80,.90,1)
# compute reclassification measures
reclassification(data=final, cOutcome=cOutcome,
                 predrisk1=predRisk1, predrisk2=predRisk2, cutoff=cutoff)



#####C-Index#####
##Load Dataset##
MESA <- data.frame(read.csv(file = '/Users/quincy/Documents/Research/Brijesh - CardioOnc/020923/NEW_CSV/test_HF.csv'))

##Impute Median Values##
library(mlr)
imputed = impute(MESA, target = character(0), classes = list(numeric = imputeMedian(), integer = imputeMedian()))
final <- as.data.frame(imputed$data)

##HF##
##Training/Testing Split##
set.seed(100)
index <- createDataPartition(final$HF, p = 0.6, list = FALSE)
train <- final[index, ]
test <- final[-index, ]

coxDemHF <- coxph(Surv(HF_Duration, HF)~Cancer_Diagnosis+Age+BMI+HTN+HLD+COPD+DM+ASA+BB+Statin+
                    ACEi_ARB_Entresto+CCB+Diuretic+insulin+NSAID+smokinghx+currentsmoking+
                    Alcohol+Race,data=final,x=TRUE,y=TRUE)
coxFuncHF <- coxph(Surv(HF_Duration, HF)~Cancer_Diagnosis+Age+BMI+HTN+HLD+COPD+DM+ASA+BB+Statin+
                     ACEi_ARB_Entresto+CCB+Diuretic+insulin+NSAID+smokinghx+currentsmoking+
                     Alcohol+Race+Func_Prob,data=final,x=TRUE,y=TRUE)
coxDrugHF <- coxph(Surv(HF_Duration, HF)~Cancer_Diagnosis+Age+BMI+HTN+HLD+COPD+DM+ASA+BB+Statin+
                     ACEi_ARB_Entresto+CCB+Diuretic+insulin+NSAID+smokinghx+currentsmoking+
                     Alcohol+Race+Drug_Prob,data=final,x=TRUE,y=TRUE)
coxOmicsHF <- coxph(Surv(HF_Duration, HF)~Cancer_Diagnosis+Age+BMI+HTN+HLD+COPD+DM+ASA+BB+Statin+
                      ACEi_ARB_Entresto+CCB+Diuretic+insulin+NSAID+smokinghx+currentsmoking+
                      Alcohol+Race+Omics_Prob,data=final,x=TRUE,y=TRUE)

###HF###
ApparrentCindexHF <- pec::cindex(list("COXPH Dem"=coxDemHF,
                                      "COXPH Function"=coxFuncHF,
                                      "COXPH Drugs"=coxDrugHF,
                                      "COXPH Omics"=coxOmicsHF),
                                 formula=Surv(HF_Duration, HF)~1,data=final,
                                 eval.times=seq(0,3754,1), pred.times=seq(0,3754,1))

col = c("#4682B4", "#AF46B4", "#B47846", "#4BB446")

plot(ApparrentCindexHF, legend = FALSE, xlim=c(0,4000), ylim=c(0.5,1.0), lwd = 5, col = col)
write.csv(ApparrentCindexHF$AppCindex, file = "/Users/quincy/Documents/Research/Brijesh - CardioOnc/020923/Concordance_HF.csv")


##Mortality Risk Prediction##
##Training/Testing Split##
set.seed(500)
index2 <- createDataPartition(final$Mortality, p = 0.6, list = FALSE)
train <- final[index2, ]
test <- final[-index2, ]

###Mortalty###
ApparrentCindexM <- pec::cindex(list("COXPH Dem"=coxDemM,
                                     "COXPH Function"=coxFuncM,
                                     "COXPH Drugs"=coxDrugM,
                                     "COXPH Omics"=coxOmicsM),
                                formula=Surv(Mortality_Duration, Mortality)~1,data=test,
                                eval.times=seq(0,3939,1), pred.times=seq(0,3939,1))

col = c("#4682B4", "#AF46B4", "#B47846", "#4BB446")

plot(ApparrentCindexM, legend = FALSE, xlim=c(0,4000), ylim=c(0.5,1.0), lwd = 5, col = col)
write.csv(ApparrentCindexM$AppCindex, file = "/Users/quincy/Documents/Research/Brijesh - CardioOnc/020923/Concordance_Mortality.csv")
