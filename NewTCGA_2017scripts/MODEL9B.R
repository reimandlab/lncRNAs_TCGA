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

rownames(canc_data) = canc_data$patient

top3genes = c("ENSG00000227544", "ENSG00000235823", "ENSG00000258082")
top3genes = top3genes[2]
canc_data = canc_data[,c(which(colnames(canc_data) %in% top3genes), 2359:ncol(canc_data))]
canc_data = canc_data[,-c(2,5)]

## set the seed to make your partition reproductible
cinds = c()
set.seed(123)

############
###START####
############

for(j in 1:100){
smp_size <- floor(0.7 * nrow(canc_data))
train_ind <- sample(seq_len(nrow(canc_data)), size = smp_size)
train <- canc_data[train_ind, ]
test <- canc_data[-train_ind, ]

###[2.] Feature selection 

###---------------------------------------------------------------
###Extract significant predictors in training data 
###---------------------------------------------------------------

###---------------------------------------------------------------
###Extract significant predictors in training data 
###---------------------------------------------------------------

source("survival_scriptJan12.R")
train[,1] = log1p(train[,1])
k = 1
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

train$time = as.numeric(train$time)
train$status[train$status=="Alive"] <- 0
train$status[train$status=="Dead"] <- 1
train$status <- as.numeric(train$status)

train = train[,-c(4)]

cox_model = coxph(Surv(time, status)  ~ ., data = train)

#create survival estimates on validation data
test = test[,which(colnames(test) %in% colnames(train))]

test$status[test$status=="Alive"] <- 0
test$status[test$status=="Dead"] <- 1
test$status <- as.numeric(test$status)
test$time <- as.numeric(test$time)

#Add high/low tags to each gene 
test[,1] = log1p(test[,1])
k = 1    
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

saveRDS(cinds, file="MODEL9B_305_OV_pats_bestGenes_CV_100timesJan19.RDS")


	
	


