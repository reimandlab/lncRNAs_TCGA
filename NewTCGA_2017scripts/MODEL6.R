###MODEL6.R

#multivariate model with 3 lncRNAs and clinical variables 

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

#check if this person is in my analysis: TCGA-61-2095

library(glmnet)
library(survcomp)
library(caret)

#Going to work on each cancer seperatley 
#for now start with one cancer 
cancer = cancers[[1]] #ovarian 
canc_data = rna[which(rna$canc == cancer),]
canc_data = as.data.frame(canc_data)

clin = canc_data[,2359:ncol(canc_data)]

#clinical variables 
add_clin = fread("OV_clinical_core.txt", data.table=F)
add_clin = add_clin[,-5]
colnames(add_clin)[1] = "patient"
#add_clin = merge(add_clin, clin, by="patient")
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
#add_clin = add_clin[,-c(6,7)]

top3genes = c("ENSG00000227544", "ENSG00000235823", "ENSG00000258082")
canc_data = canc_data[,c(which(colnames(canc_data) %in% top3genes), 2359:ncol(canc_data))]
canc_data = canc_data[,-c(4,7)]

add_clin = merge(add_clin, canc_data, by="patient")

## set the seed to make your partition reproductible
cinds = c()
set.seed(123)

############
###START####
############

for(j in 1:100){
smp_size <- floor(0.7 * nrow(add_clin))
train_ind <- sample(seq_len(nrow(add_clin)), size = smp_size)
train <- add_clin[train_ind, ]
test <- add_clin[-train_ind, ]

###[2.] Feature selection 

###---------------------------------------------------------------
###Extract significant predictors in training data 
###---------------------------------------------------------------

###---------------------------------------------------------------
###Extract significant predictors in training data 
###---------------------------------------------------------------

source("survival_scriptJan12.R")
train[,5:7] = log1p(train[,5:7])

for(k in 5:7){
    median2 <- quantile(as.numeric(train[,k]), 0.5)
    if(median2 ==0){
    median2 = mean(as.numeric(train[,k]))
    }
    for(m in 1:nrow(train)){
    genexp <- train[m,k]
    if(genexp >= median2){
      train[m,k] <- 1
      }
    if(genexp < median2){
      train[m,k] <- 0
      }
    } 
}

train$time = as.numeric(train$time)
train$status[train$status=="Alive"] <- 0
train$status[train$status=="Dead"] <- 1
train$status <- as.numeric(train$status)
train$age = as.numeric(train$age)
train$stage = as.numeric(train$stage)
train$grade = as.numeric(train$grade)
train = train[,-c(1)]

cox_model = coxph(Surv(time, status)  ~ ., data = train)

#create survival estimates on validation data
test = test[,which(colnames(test) %in% colnames(train))]

test$status[test$status=="Alive"] <- 0
test$status[test$status=="Dead"] <- 1
test$status <- as.numeric(test$status)
test$time <- as.numeric(test$time)
test$age = as.numeric(test$age)
test$stage = as.numeric(test$stage)
test$grade = as.numeric(test$grade)

#Add high/low tags to each gene 
test[,4:6] = log1p(test[,4:6])
for(k in 4:6){
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
cinds = c(cinds, cindex_validation)
print(cindex_validation)
}

print("done")


###---------------------------------------------------------------
###Survival function 
###---------------------------------------------------------------

saveRDS(cinds, file="MODEL6_202_OV_pats_3bestGenesANDclinicalVariables_CV_100timesJan19.RDS")


	
	


