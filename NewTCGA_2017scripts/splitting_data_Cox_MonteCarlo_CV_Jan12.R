###source_code_Cox_MonteCarlo_CV_Jan12.R

###Purpose--------------------------------------------------------------------

#Using each cancer's dataset, divide it into training and test sets (MonteCarlo, 100 times)
#Each time, run univariate survival analysis to identify individual lncRNAs that are significantly
#associated with survival 

#retained feature number shouldn't exceed number of events (deaths) in training set 

#Then Run Lasso --> feature selection, keep track of this list during each iteration 
#penalty parameter was chosen based on the fivefold cross-validation within the training set 

#apply trained model to the test set for prediction 
#calculate C-index using survcomp package 

###[1.] Load requirements  

###---------------------------------------------------------------
###Load libraries and data 
###---------------------------------------------------------------

source("source_code_Cox_MonteCarlo_CV_Jan12.R")
require(caTools)
set.seed(101) 

###[2.] Data splitting 

###---------------------------------------------------------------
###Split dataset into training and testing 
###---------------------------------------------------------------

#Going to work on each cancer seperatley 
#for now start with one cancer 
cancer = cancers[[1]]
canc_data = rna[canc == cancer]
canc_data = as.data.frame(canc_data)
rownames(canc_data) = canc_data$patient
gene_data = t(canc_data[,1:(ncol(rna)-5)])
	
#1. remove any genes that have 0 counts within cancer
sums = apply(canc_data[,1:(ncol(rna)-5)], 2, sum) #134 with 0 expression in ALL patients 
zeroes = names(sums[which(sums ==0)]) #what are they?
#in pcawg 
z <- which(colnames(canc_data) %in% zeroes)
if(!(length(z)==0)){
	canc_data = canc_data[,-z]
}

smp_size <- floor(0.8 * nrow(canc_data))
## set the seed to make your partition reproductible
set.seed(123)
train_ind <- sample(seq_len(nrow(canc_data)), size = smp_size)

train <- canc_data[train_ind, ]
test <- canc_data[-train_ind, ]

###[2.] Feature selection 

###---------------------------------------------------------------
###Extract significant predictors in training data 
###---------------------------------------------------------------

source("survival_scriptJan12.R")

#surv_test is the loaded survival function which takes in gene 
#name as input and conducts a univariate cox regression model 

genes = as.list(colnames(train)[1:5778])
tests_survival2 = llply(genes, surv_test, .progress="text")
tests_survival3 = ldply(tests_survival2, rbind)
tests_survival3$pval = as.numeric(tests_survival3$pval)

#keep sig, how many ? don't want to be more than # of deaths
#so filter by HR and keep the max 
tests_survival3 = subset(tests_survival3, pval <=0.05)
length(unique(tests_survival3$gene))

#using LASSO 
#with 10fold cross validation to get optimal lamda penalization parameter
z <- which(colnames(canc_data) %in% tests_survival3$gene)
canc_data = canc_data[,c(z, 5779:5783)]
#practice lasso 
canc_data$time = as.numeric(canc_data$time)
canc_data$status[canc_data$status=="Alive"] <- 0
canc_data$status[canc_data$status=="Dead"] <- 1
canc_data$status <- as.numeric(canc_data$status)

#x is an n*p matrix of covariate values, each row corresponds to a patient and each
#column is a covariate 
#y is an n*2 matrix with a column time and status 

gene_data = canc_data[,1:(ncol(canc_data)-5)]
for(i in 1:ncol(gene_data)){
	gene_data[,i] = scale(gene_data[,i])[, 1]
}

x <- model.matrix( ~ ., gene_data)
y <- Surv(canc_data$time, canc_data$status)

library(glmnet)
fit = glmnet(x, y, family = "cox")

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


#limit survival analysis to just those genes in training and test set 
cox_model = coxph(Surv(training_data$Survival,training_data$Status) ~ ., data=training_data) 

#create survival estimates on validation data
pred_validation = predict(cox_model, newdata = test)

# Determine concordance
cindex_validation = concordance.index (pred_validation, surv.time = validation_data$Survival,
                                       surv.event=validation_data$Status, method = "noether")

###---------------------------------------------------------------
###Survival function 
###---------------------------------------------------------------




	
	


