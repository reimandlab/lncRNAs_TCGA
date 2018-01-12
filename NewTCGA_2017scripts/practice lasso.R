





#practice lasso 
canc_data$time = as.numeric(canc_data$time)
canc_data$status[canc_data$status=="Alive"] <- 0
canc_data$status[canc_data$status=="Dead"] <- 1
canc_data$status <- as.numeric(canc_data$status)

#x is an n*p matrix of covariate values, each row corresponds to a patient and each
#column is a covariate 

#y is an n*2 matrix with a column time and status 

gene_data = log1p(canc_data[,1:5778])

x <- model.matrix( ~ ., gene_data)
y <- Surv(canc_data$time, canc_data$status)

library(glmnet)

fit = glmnet(x, y, family = "cox")

coef(fit, s = 0.05)

#Since the Cox Model is not commonly used for prediction, we do not give an illustrative example on prediction. 
#If needed, users can refer to the help file by typing help(predict.glmnet).
#Also, the function cv.glmnet can be used to compute k
#-fold cross-validation for the Cox model. The usage is similar to that for other families except for two main differences.

cvfit = cv.glmnet(x, y, family = "cox", alpha =1, nfolds=5) #uses cross validation to select
#the best lamda and then use lambda to see which features remain in model 

cvfit$lambda.min #left vertical line
cvfit$lambda.1se #right vertical line 

#active covariates and their coefficients 
coef.min = coef(cvfit, s = "lambda.min") 
active.min = which(coef.min != 0)
index.min = coef.min[active.min]

index.min
coef.min