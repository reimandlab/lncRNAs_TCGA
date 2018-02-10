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

#start with only lncRNA_intergenic
lincs = subset(fantom, CAT_geneClass == "lncRNA_intergenic")
z = which(colnames(rna) %in% lincs$gene)
rna = as.data.frame(rna)
rna = rna[,c(z, 5786:5790)]

###[2.] Data splitting 

###---------------------------------------------------------------
###Split dataset into training and testing 
###---------------------------------------------------------------

library(glmnet)
library(survcomp)
library(caret)

#Going to work on each cancer seperatley 
#for now start with one cancer 
cancer = cancers[[1]]
canc_data = rna[which(rna$canc == cancer),]
canc_data = as.data.frame(canc_data)

#clinical variables 
add_clin = fread("OV_clinical_core.txt", data.table=F)
clin = canc_data[,2359:ncol(canc_data)]

#clinical variables 
add_clin = fread("OV_clinical_core.txt", data.table=F)
add_clin = add_clin[,-5]
colnames(add_clin)[1] = "patient"
add_clin = merge(add_clin, clin, by="patient")
rownames(add_clin) = add_clin$patient
z = which(add_clin$stage == "[Not Available]")
add_clin = add_clin[-z,]
z = which(add_clin$grade %in% c("G4", "GB", "GX"))
add_clin = add_clin[-z,]

add_clin$stage[add_clin$stage == "IIA"] = 2
add_clin$stage[add_clin$stage == "IIB"] = 2
add_clin$stage[add_clin$stage == "IIC"] = 2
add_clin$stage[add_clin$stage == "IIIA"] = 3
add_clin$stage[add_clin$stage == "IIIB"] = 3
add_clin$stage[add_clin$stage == "IIIC"] = 3
add_clin$stage[add_clin$stage == "IV"] = 4

add_clin$grade[add_clin$grade == "G1"] = 1
add_clin$grade[add_clin$grade == "G2"] = 2
add_clin$grade[add_clin$grade == "G3"] = 3

rownames(canc_data) = canc_data$patient
canc_data = subset(canc_data, patient %in% add_clin$patient)
add_clin = add_clin[,-c(6,7)]

gene_data = t(canc_data[,1:(ncol(rna)-5)])
	
#1. remove any genes that have 0 counts within cancer
sums = apply(canc_data[,1:(ncol(rna)-5)], 2, sum) #134 with 0 expression in ALL patients 
zeroes = names(sums[which(sums ==0)]) #what are they?
#in pcawg 
z <- which(colnames(canc_data) %in% zeroes)
if(!(length(z)==0)){
	canc_data = canc_data[,-z]
}

## set the seed to make your partition reproductible
cinds = c()
genes_results = list()
set.seed(123)

for(j in 1:100){
smp_size <- floor(0.7 * nrow(canc_data))
train_ind <- sample(seq_len(nrow(canc_data)), size = smp_size)
train <- canc_data[train_ind, ]
test <- canc_data[-train_ind, ]

sums = apply(train[,1:(ncol(train)-5)], 2, sum) #134 with 0 expression in ALL patients 
zeroes = names(sums[which(sums ==0)]) #what are they?
#in pcawg 
z <- which(colnames(train) %in% zeroes)
if(!(length(z)==0)){
	train = train[,-z]
}

###[2.] Feature selection 

###---------------------------------------------------------------
###Extract significant predictors in training data 
###---------------------------------------------------------------

source("survival_scriptJan12.R")

#surv_test is the loaded survival function which takes in gene 
#name as input and conducts a univariate cox regression model 
genes = as.list(colnames(train)[1:(ncol(train)-5)])
#genes = genes[1:100]
tests_survival2 = llply(genes, surv_test, .progress="text")
tests_survival3 = ldply(tests_survival2, rbind)
tests_survival3$pval = as.numeric(tests_survival3$pval)

#keep sig, how many ? don't want to be more than # of deaths
#so filter by HR and keep the max 
tests_survival3 = subset(tests_survival3, pval <=0.05)
length(unique(tests_survival3$gene))
sortme <- function(dt, sort.field) dt[order(-abs(dt[[sort.field]]))]
tests_survival3$HR = as.numeric(tests_survival3$HR)
tests_survival3 = as.data.table(tests_survival3)
sort.field = "HR"
tests_survival3 = sortme(tests_survival3, sort.field)
#z = table(train$status)[2]
tests_survival3 = as.data.frame(tests_survival3)
#tests_survival3 = tests_survival3[1:z, ]

#using LASSO 
#with 10fold cross validation to get optimal lamda penalization parameter
z <- which(colnames(train) %in% tests_survival3$gene)
train = train[,c(z, (ncol(train)-4):ncol(train))]
#practice lasso 
train$time = as.numeric(train$time)
train$status[train$status=="Alive"] <- 0
train$status[train$status=="Dead"] <- 1
train$status <- as.numeric(train$status)

#x is an n*p matrix of covariate values, each row corresponds to a patient and each
#column is a covariate 
#y is an n*2 matrix with a column time and status 

train$patient = rownames(train)
train = merge(train, add_clin, by = c("patient", "canc", "sex"))

gene_data = train[,4:(ncol(train)-5)]
gene_data = log1p(gene_data)

#Add high/low tags to each gene 
for(k in 1:ncol(gene_data)){
	median2 <- quantile(as.numeric(gene_data[,k]), 0.5)
  	if(median2 ==0){
  	median2 = mean(as.numeric(gene_data[,k]))
  	}
  	for(m in 1:nrow(gene_data)){
    genexp <- gene_data[m,k]
    if(genexp >= median2){
      gene_data[m,k] <- 1
      }
    if(genexp < median2){
      gene_data[m,k] <- 0
      }
    } 
}


#gene_data = log1p(gene_data)
x <- model.matrix( ~., gene_data)
y <- Surv(train$time, train$status)

fit = glmnet(x, y, family = "cox")

#Since the Cox Model is not commonly used for prediction, we do not give an illustrative example on prediction. 
#If needed, users can refer to the help file by typing help(predict.glmnet).
#Also, the function cv.glmnet can be used to compute k
#-fold cross-validation for the Cox model. The usage is similar to that for other families except for two main differences.

cvfit = cv.glmnet(x, y, family = "cox", alpha =1) #uses cross validation to select

#the best lamda and then use lambda to see which features remain in model 
cvfit$lambda.min #left vertical line
cvfit$lambda.1se #right vertical line 

#active covariates and their coefficients 
coef.min = coef(cvfit, s = "lambda.min") 
active.min = which(coef.min != 0)
index.min = coef.min[active.min]

genes_keep = rownames(coef.min)[active.min]
print(genes_keep)

if(!(length(genes_keep)==0)){

training_data = train[,c(1, which(colnames(train) %in% genes_keep), (ncol(train)-4):ncol(train))]

#limit survival analysis to just those genes in training and test set 
training_data$status[training_data$status=="Alive"] <- 0
training_data$status[training_data$status=="Dead"] <- 1
training_data$status <- as.numeric(training_data$status)
training_data$time <- as.numeric(training_data$time)
training_data$age <- as.numeric(training_data$age)
training_data$grade <- as.numeric(training_data$grade)
training_data$stage <- as.numeric(training_data$stage)
rownames(training_data) = training_data$patient 
training_data = training_data[,-1]

training_data[,1:(ncol(training_data)-5)] = log1p(training_data[,1:(ncol(training_data)-5)] )
#Add high/low tags to each gene 
for(k in 1:(ncol(training_data)-5)){
	median2 <- quantile(as.numeric(training_data[,k]), 0.5)
  	if(median2 ==0){
  	median2 = mean(as.numeric(training_data[,k]))
  	}
  	for(m in 1:nrow(training_data)){
    genexp <- training_data[m,k]
    if(genexp >= median2){
      training_data[m,k] <- 1
      }
    if(genexp < median2){
      training_data[m,k] <- 0
      }
    } 
}

cox_model = coxph(Surv(time, status)  ~ ., data = training_data)

#create survival estimates on validation data
test$patient = rownames(test)
test = merge(test, add_clin, by = c("patient", "canc", "sex"))

test = test[,which(colnames(test) %in% colnames(training_data))]

test$status[test$status=="Alive"] <- 0
test$status[test$status=="Dead"] <- 1
test$status <- as.numeric(test$status)
test$time <- as.numeric(test$time)
test$age <- as.numeric(test$age)
test$grade <- as.numeric(test$grade)
test$stage <- as.numeric(test$stage)

test[,1:(ncol(test)-5)] = log1p(test[,1:(ncol(test)-5)])
#Add high/low tags to each gene 
for(k in 1:(ncol(test)-5)){
	median2 <- quantile(as.numeric(test[,k]), 0.5)
  	if(median2 ==0){
  	median2 = mean(as.numeric(test[,k]))
  	}
  	for(m in 1:nrow(test)){
    genexp <- test[m,k]
    if(genexp >= median2){
      test[m,k] <- 1
      }
    if(genexp < median2){
      test[m,k] <- 0
      }
    } 
}	
pred_validation = predict(cox_model, newdata = test)

# Determine concordance
cindex_validation = concordance.index(pred_validation, surv.time = test$time,
                                       surv.event=test$status, method = "noether")

cindex_validation = cindex_validation$c.index
genes_results[j] = list(as.list(genes_keep))
cinds = c(cinds, cindex_validation)
print(cindex_validation)
}
}

print("done")

###---------------------------------------------------------------
###Survival function 
###---------------------------------------------------------------

saveRDS(genes_results, file="lncRNAS_AND_wclinical_updated_code_binary_predictors_list_of_sig_genes_CVJan16.rds")
saveRDS(cinds, file="lncRNAS_AND_wclinical_updated_code_binary_predictors_cindes_CV_100timesJan16.RDS")


	
	


