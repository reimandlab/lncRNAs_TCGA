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

source("source_code_Cox_MonteCarlo_CV_PAAD_Feb5.R")
require(caTools)

#start with only lncRNA_intergenic
lincs = subset(fantom, (CAT_geneClass == "lncRNA_intergenic") & (CAT_geneCategory %in% c("e_lncRNA", "p_lncRNA_intergenic")))
z = which(colnames(rna) %in% lincs$gene)
rna = as.data.frame(rna)
rna = rna[,c(z, 5749:5753)]

###[2.] Data splitting 

###---------------------------------------------------------------
###Split dataset into training and testing 
###---------------------------------------------------------------

#check if this person is in my analysis: TCGA-61-2095

library(glmnet)
library(survcomp)
library(caret)
library(stringr)

#Going to work on each cancer seperatley 
#for now start with one cancer 
cancer = cancers[[1]] #PAAD 
canc_data = rna[which(rna$canc == cancer),]
canc_data = as.data.frame(canc_data)

clin = canc_data[,1513:ncol(canc_data)]

#New clinical file from Firehose 
newclin = readRDS("185_PAAD_pats_clinical_data_Feb5_firehose.rds")
z = which(newclin$patient.bcr_patient_barcode %in% clin$patient)
newclin = newclin[z,]
colss = colnames(newclin)
colsskeep = c(which(str_detect(colss, "stage")), which(str_detect(colss, "grade")), which(str_detect(colss, "days")), which(str_detect(colss, "followup")))
newclin = newclin[,c(1:25, colsskeep)]

#see if time and status mathch old information 
clin$newtime = ""
clin$newstatus =""
for(i in 1:nrow(newclin)){
  pat = newclin$patient.bcr_patient_barcode[i]
  z = which(clin$patient %in% pat)
  t = newclin$patient.days_to_death[i]
  if(!(is.na(t))){
    clin$newtime[z] = t
    clin$newstatus[z] = "Dead"
  }
  if(is.na(t)){
    clin$newtime[z] = newclin$patient.days_to_last_followup[i]
    clin$newstatus[z] = "Alive"
  }
}

z = which(is.na(clin$time))
#clin = clin[-z,]
ind <- clin$time == clin$newtime
clin[ind, "newtime"] <- "EQUAL"
ind <- clin$status == clin$newstatus
clin[ind, "newstatus"] <- "EQUAL"

clin = subset(clin, time >0)

rownames(canc_data) = canc_data$patient
canc_data = subset(canc_data, patient %in% clin$patient)

clin$grade = ""
clin$stage = ""
for(i in 1:nrow(clin)){
  pat = clin$patient[i]
  z =which(newclin$patient.bcr_patient_barcode ==pat)
  clin$stage[i] = newclin$patient.stage_event.pathologic_stage[z]
  clin$grade[i] = newclin$patient.neoplasm_histologic_grade[z]
}

#173 PAAD patients with clinical data 

#change stage and grade to numeric values
clin$grade[clin$grade == "g1"] = 1
clin$grade[clin$grade == "g2"] = 2
clin$grade[clin$grade == "g3"] = 3
clin$grade[clin$grade == "g4"] = 4
z = which(clin$grade == "gx")
clin = clin[-z,]

clin$stage[clin$stage == "stage i"] = 1
clin$stage[clin$stage == "stage ia"] = 1
clin$stage[clin$stage == "stage ib"] = 1

clin$stage[clin$stage == "stage iia"] = 2
clin$stage[clin$stage == "stage iib"] = 2

clin$stage[clin$stage == "stage iii"] = 3
clin$stage[clin$stage == "stage iv"] = 4

z = which(is.na(clin$stage))
clin = clin[-z,]
z = which(is.na(clin$grade))
#clin = clin[-z,]

clin$grade = as.numeric(clin$grade)
clin$stage = as.numeric(clin$stage) #169 PAAD patients remain 

## set the seed to make your partition reproductible
set.seed(914)

############
###START####
############

genes_results = list()
cinds = c()
for(j in 1:100){
smp_size <- floor(0.7 * nrow(canc_data))
train_ind <- sample(seq_len(nrow(canc_data)), size = smp_size)
train <- canc_data[train_ind, c(1:(ncol(canc_data)-5), ncol(canc_data))]
test <- canc_data[-train_ind, c(1:(ncol(canc_data)-5), ncol(canc_data))]

###[2.] Feature selection 

###---------------------------------------------------------------
###Extract significant predictors in training data 
###---------------------------------------------------------------

source("survival_scriptJan23_newClindata_OV.R")
means = apply(train[,1:1512], 2, mean) #134 with 0 expression in ALL patients 
zeroes = names(means[which(means ==0)]) #what are they?
z <- which(colnames(train) %in% zeroes)
if(!(length(z)==0)){
  train = train[,-z]
}
train[,1:(ncol(train)-1)] = log1p(train[,1:(ncol(train)-1)])
sums = apply(train[,1:(ncol(train)-1)], 2, sum) #134 with 0 expression in ALL patients 
zeroes = names(sums[which(sums ==0)]) #what are they?
z <- which(colnames(train) %in% zeroes)
if(!(length(z)==0)){
  train = train[,-z]
}
means = apply(train[,1:(ncol(train)-1)], 2, mean) #134 with 0 expression in ALL patients 
zeroes = names(means[which(means <1)]) #what are they?
z <- which(colnames(train) %in% zeroes)
if(!(length(z)==0)){
  train = train[,-z]
}

train = merge(train, clin, by=c("patient"))
rownames(train) = train$patient
train = train[,-1]

test$patient = rownames(test)
test = merge(test, clin, by = c("patient"))
rownames(test) = test$patient
test = test[,-1]

#surv_test is the loaded survival function which takes in gene 
#name as input and conducts a univariate cox regression model 
genes = as.list(colnames(train)[1:(ncol(train)-8)])
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
train = train[,c(z, (ncol(train)-6), (ncol(train)-5))]
#practice lasso 
train$time = as.numeric(train$time)
train$status[train$status=="Alive"] <- 0
train$status[train$status=="Dead"] <- 1
train$status <- as.numeric(train$status)

#x is an n*p matrix of covariate values, each row corresponds to a patient and each
#column is a covariate 
#y is an n*2 matrix with a column time and status 

gene_data = as.data.frame(train[,1:(ncol(train)-2)])
colnames(gene_data) = colnames(train)[1:(ncol(train)-2)]
means = apply(gene_data, 2, mean)

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

#there might be a imbalance between how many patients in each group after dividing patients by the mean 
genes = as.list(colnames(gene_data))

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

genes_results[j] = list(as.list(genes_keep))
train = train[,c(which(colnames(train) %in% c(genes_keep,"time", "status")))]
cox_model = coxph(Surv(time, status)  ~ ., data = train)

#create survival estimates on validation data
test = test[,which(colnames(test) %in% colnames(train))]
test$status[test$status=="Alive"] <- 0
test$status[test$status=="Dead"] <- 1
test$status <- as.numeric(test$status)
test$time <- as.numeric(test$time)

#Add high/low tags to each gene 
test[,1:(ncol(test)-2)] = log1p(test[,1:(ncol(test)-2)])
for(k in 1:(ncol(test)-2)){
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

#if predicotrs just lncRNAs 
pred_validation = predict(cox_model, newdata = test)
# Determine concordance
cindex_validation = concordance.index(pred_validation, surv.time = test$time,
                                       surv.event=test$status, method = "noether")

cindex_validation = cindex_validation$c.index
cinds = c(cinds, cindex_validation)
print(cindex_validation)
print(paste("run is", j))
}



saveRDS(genes_results, file="ALL_PAAD_just_intergenic_pats_binary_predictors_list_of_sig_genes_CV_5folds_March8.rds")




































	
	


