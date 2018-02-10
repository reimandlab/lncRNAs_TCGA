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
library(stringr)

#Going to work on each cancer seperatley 
#for now start with one cancer 
cancer = cancers[[2]] #KIRC 
canc_data = rna[which(rna$canc == cancer),]
canc_data = as.data.frame(canc_data)

clin = canc_data[,2359:ncol(canc_data)]

#New clinical file from Firehose 
newclin = readRDS("537_KIRC_pats_clinical_data_Jan23_firehose.rds")
z = which(newclin$patient.bcr_patient_barcode %in% clin$patient)
newclin = newclin[z,]
colss = colnames(newclin)
colsskeep = c(which(str_detect(colss, "stage")), which(str_detect(colss, "grade")))
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

clin = subset(clin, time >0)

rownames(canc_data) = canc_data$patient
canc_data = subset(canc_data, patient %in% clin$patient)
colnames(newclin)[11] = "patient"

clin$grade = ""
clin$stage = ""
for(i in 1:nrow(clin)){
  pat = clin$patient[i]
  z =which(newclin$patient ==pat)
  clin$stage[i] = newclin$patient.stage_event.pathologic_stage[z]
  clin$grade[i] = newclin$patient.neoplasm_histologic_grade[z]
}

## set the seed to make your partition reproductible
set.seed(123)

############
###START####
############

genes_results = list()

for(j in 1:1000){
smp_size <- floor(0.7 * nrow(canc_data))
train_ind <- sample(seq_len(nrow(canc_data)), size = smp_size)
train <- canc_data[train_ind, c(1:(ncol(canc_data)-5), ncol(canc_data))]
test <- canc_data[-train_ind, c(1:(ncol(canc_data)-5), ncol(canc_data))]

###[2.] Feature selection 

###---------------------------------------------------------------
###Extract significant predictors in training data 
###---------------------------------------------------------------

###---------------------------------------------------------------
###Extract significant predictors in training data 
###---------------------------------------------------------------

source("survival_scriptJan23_newClindata_OV.R")
train[,1:2358] = log1p(train[,1:2358])
sums = apply(train[,1:2358], 2, sum) #134 with 0 expression in ALL patients 
zeroes = names(sums[which(sums ==0)]) #what are they?
z <- which(colnames(train) %in% zeroes)
if(!(length(z)==0)){
  train = train[,-z]
}
train = merge(train, clin, by=c("patient"))
rownames(train) = train$patient
train = train[,-1]

#surv_test is the loaded survival function which takes in gene 
#name as input and conducts a univariate cox regression model 
genes = as.list(colnames(train)[1:(ncol(train)-8)])
genes = genes[1:100]

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

gene_data = train[,1:(ncol(train)-2)]

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
genes_results[j] = list(as.list(genes_keep))

}

saveRDS(genes_results, file="ALL_KIRC_pats473_binary_predictors_list_of_sig_genes_CV1000Jan23_batch1.rds")




































	
	


