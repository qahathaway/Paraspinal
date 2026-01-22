
### Linear Regression and Sensitivity Analsyis ###
library(sensemakr)
library(mlr)

my_data <- read.csv("/path/to/data/*.csv", header=TRUE)

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
               Moderate_Phys + Vigorous_Phys + AGATSTON + Median, data = my_data)

model5 <- lm(vBMD ~Age + Gender + Race + Site + BMI + HTN + Diabetes +
               Liver + Kidney + Creatinine + Cancer + Steroids + Bisphosphonates + Warfarin +
               Thyroid + Alcohol + Alcohol_Consumption +
               Tobacco + Pack_Years +
               Moderate_Phys + Vigorous_Phys + AGATSTON + Shape + Median, data = my_data)

model6 <- lm(vBMD ~Age + Gender + Race + Site + BMI + HTN + Diabetes +
               Liver + Kidney + Creatinine + Cancer + Steroids + Bisphosphonates + Warfarin +
               Thyroid + Alcohol + Alcohol_Consumption +
               Tobacco + Pack_Years + Manufacturer + Model +
               Moderate_Phys + Vigorous_Phys + AGATSTON + Shape + Median, data = my_data)

anova(model1, model2)
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
                         treatment = "Median",
                         benchmark_covariates = "Bisphosphonates",
                         kd = 1:3)

ovb_minimal_reporting(sensitivity, format = "html")

# Print and plot a summary of the results
summary(sensitivity)
plot(sensitivity)
plot(sensitivity, sensitivity.of = "t-value")
plot(sensitivity, type = "extreme")
