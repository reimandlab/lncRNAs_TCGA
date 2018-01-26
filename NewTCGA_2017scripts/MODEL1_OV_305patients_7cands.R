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
###Load Libraries
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
library(stringr)

#Going to work on each cancer seperatley 
#for now start with one cancer 
cancer = cancers[[1]] #ovarian 
canc_data = rna[which(rna$canc == cancer),]
canc_data = as.data.frame(canc_data)

clin = canc_data[,2359:ncol(canc_data)]

#New clinical file from Firehose 
newclin = readRDS("591_OV_pats_clinical_data_Jan23_firehose.rds")
z = which(newclin$patient.bcr_patient_barcode %in% clin$patient)
newclin = newclin[z,]
colss = colnames(newclin)
colsskeep = which(str_detect(colss, "age"))[1]
colsskeep = c(colsskeep, which(str_detect(colss, "stage")), which(str_detect(colss, "grade")))
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

ind <- clin$time == clin$newtime
clin[ind, "newtime"] <- "EQUAL"
ind <- clin$status == clin$newstatus
clin[ind, "newstatus"] <- "EQUAL"

rownames(canc_data) = canc_data$patient
canc_data = subset(canc_data, patient %in% clin$patient)
colnames(newclin)[12] = "patient"

clin$grade = ""
clin$stage = ""
clin$age = ""
for(i in 1:nrow(clin)){
  pat = clin$patient[i]
  z =which(newclin$patient ==pat)
  clin$stage[i] = newclin$patient.stage_event.clinical_stage[z]
  clin$grade[i] = newclin$patient.neoplasm_histologic_grade[z]
  clin$age[i] = newclin$patient.age_at_initial_pathologic_diagnosis.1[z]
}

z = which(clin$grade %in% c("gb", "gx"))
clin = clin[-z,]

z = which(is.na(clin$stage))
clin = clin[-z,]

z = which(is.na(clin$grade))
clin = clin[-z,]

clin$stage[clin$stage == "stage ic"] = 1
clin$stage[clin$stage %in% c("stage iia", "stage iib", "stage iic")] = 2
clin$stage[clin$stage %in% c("stage iiia", "stage iiib", "stage iiic")] = 3
clin$stage[clin$stage == "stage iv"] = 4

clin$grade[clin$grade == "g1"] = 1
clin$grade[clin$grade == "g2"] = 2
clin$grade[clin$grade == "g3"] = 3
clin$grade[clin$grade == "g4"] = 4

###---------------------------------------------------------------
###Cands - 10 batches of 1000 CV 
###---------------------------------------------------------------

genes = readRDS("final_candidates_10batches_of1000CV_305ovarian_cancer_patientsJan23.RDS")
cands = as.data.frame(table(genes$gene))
colnames(cands)[1] = "gene"
cands$name = ""
for(i in 1:nrow(cands)){
  n = fantom$CAT_geneName[which(fantom$gene %in% cands$gene[i])]
  cands$name[i] = n
}


###set up for validation 

top3genes = as.character(cands$gene)
canc_data = canc_data[,c(which(colnames(canc_data) %in% top3genes), 2359:ncol(canc_data))]
canc_data = canc_data[,-c(8,11)] #dont need cancer and sex columns

clin = merge(clin, canc_data, by=c("patient", "time", "status"))

## set the seed to make your partition reproductible
cinds = c()
set.seed(123)

############
###START####
############

for(j in 1:1000){
smp_size <- floor(0.7 * nrow(clin))
train_ind <- sample(seq_len(nrow(clin)), size = smp_size)
train <- clin[train_ind, ]
test <- clin[-train_ind, ]

###[2.] Feature selection 

###---------------------------------------------------------------
###Extract significant predictors in training data 
###---------------------------------------------------------------

###---------------------------------------------------------------
###Extract significant predictors in training data 
###---------------------------------------------------------------

source("survival_scriptJan23_newClindata_OV.R")
train[,11:ncol(train)] = log1p(train[,11:ncol(train)])

for(k in 11:ncol(train)){
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
train = train[,-c(1, 4, 5, 6, 7)] #remove cancer, patient and sex columns, newtime and newstatus 

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
test[,6:ncol(test)] = log1p(test[,6:ncol(test)])
for(k in 6:ncol(test)){
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

saveRDS(cinds, file="MODEL1_305_OV_pats_CAND_genes_AND_STAGE_GRADE_AGE_CV_1000timesJan26.RDS")


	
	


