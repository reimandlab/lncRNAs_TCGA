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

source("source_code_Cox_MonteCarlo_CV_Mar13.R")
require(caTools)
#check if this person is in my analysis: TCGA-61-2095
library(glmnet)
library(survcomp)
library(caret)
library(stringr)

#start with only lncRNA_intergenic
#lincs = subset(fantom, (CAT_geneClass == "lncRNA_intergenic") & (CAT_geneCategory %in% c("e_lncRNA", "p_lncRNA_intergenic")))

lnc_info = read.csv("fantom_genebased_evidence_supp_table_17.csv")
lnc_info = lnc_info[which(str_detect(lnc_info$CAT_geneID, "ENSG")),]
lnc_info = subset(lnc_info, lnc_info$CAT_geneCategory %in% c("e_lncRNA", "p_lncRNA_divergent", "p_lncRNA_intergenic"))

#shorten the gene names 
extract3 <- function(row){
  gene <- as.character(row[[1]])
  ens <- gsub("\\..*","",gene)
  return(ens)
}
lnc_info[,1] <- apply(lnc_info[,1:2], 1, extract3) #5049 lncRNAs 
lnc_info = subset(lnc_info, (CAT_geneClass == "lncRNA_intergenic") & (CAT_geneCategory %in% c("e_lncRNA", "p_lncRNA_intergenic")))
#2091 

#how many conserved 
#1207 have conserved exons 
z1 = which(!(lnc_info$exon_RS_score == "__na"))

#924 have conserved TIRs 
z2 = which(!(lnc_info$TIR_RS_score == "__na"))

conserved = unique(c(z1, z2))
lnc_info = lnc_info[conserved, ]
colnames(lnc_info)[1] = "gene"

#z = which(colnames(rna) %in% lnc_info$gene) <-------- 
rna = as.data.frame(rna)
#rna = rna[,c(z, (ncol(rna)-5):ncol(rna))] <---------

###[2.] Data splitting 

###---------------------------------------------------------------
###Split dataset into training and testing 
###---------------------------------------------------------------

rna = subset(rna, canc %in% c("Kidney renal clear cell carcinoma", "Liver hepatocellular carcinoma", 
  "Ovarian serous cystadenocarcinoma", "Pancreatic adenocarcinoma"))

#remove 0 sums
sums = apply(rna[,1:(ncol(rna)-5)], 2, sum)
z = which(sums==0)
#rna = rna[,-(which(colnames(rna) == names(z)))]

#Going to work on each cancer seperatley 
#for now start with one cancer
cancers = cancers[which(cancers %in% rna$canc)] 
cancer = cancers[[3]] #KIRC 
canc_data = rna[which(rna$canc == cancer),]
canc_data = as.data.frame(canc_data)

clin = canc_data[,(ncol(canc_data)-4):ncol(canc_data)]

#New clinical file from Firehose 
newclin = readRDS("537_KIRC_pats_clinical_data_Jan23_firehose.rds")
z = which(newclin$patient.bcr_patient_barcode %in% clin$patient)
newclin = newclin[z,]
colss = colnames(newclin)
colsskeep = c(which(str_detect(colss, "stage")), which(str_detect(colss, "age")), which(str_detect(colss, "grade")), which(str_detect(colss, "days")), which(str_detect(colss, "followup")))
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
clin$age = ""
for(i in 1:nrow(clin)){
  pat = clin$patient[i]
  z =which(newclin$patient.bcr_patient_barcode ==pat)
  clin$stage[i] = newclin$patient.stage_event.pathologic_stage[z]
  clin$grade[i] = newclin$patient.neoplasm_histologic_grade[z]
  clin$age[i] = newclin$patient.age_at_initial_pathologic_diagnosis[z]
}

#change stage and grade to numeric values
clin$grade[clin$grade == "g1"] = 1
clin$grade[clin$grade == "g2"] = 2
clin$grade[clin$grade == "g3"] = 3
clin$grade[clin$grade == "g4"] = 4
z = which(clin$grade == "gx")
clin = clin[-z,]

clin$stage[clin$stage == "stage i"] = 1
clin$stage[clin$stage == "stage ii"] = 2
clin$stage[clin$stage == "stage iii"] = 3
clin$stage[clin$stage == "stage iv"] = 4

z = which(is.na(clin$stage))
clin = clin[-z,]
z = which(is.na(clin$grade))
clin = clin[-z,]

clin$grade = as.numeric(clin$grade)
clin$stage = as.numeric(clin$stage)
clin$age = as.numeric(clin$age)
clin$time = as.numeric(clin$time)
clin$status[clin$status=="Alive"] <- 0
clin$status[clin$status=="Dead"] <- 1

#how many lncRNAs have median < 1? 
canc_data[,1:(ncol(canc_data)-5)] = log1p(canc_data[,1:(ncol(canc_data)-5)])

detectable = readRDS("PCAWG_detectable_genes_4cancers_March20.rds")
detectable = subset(detectable, canc == "Kidney-RCC")
detectable$canc = "Kidney renal clear cell carcinoma"

z = which(colnames(canc_data) %in% c(detectable$gene, "canc", "time", "status", "sex", "patient"))
canc_data = canc_data[,z]

## set the seed to make your partition reproductible
set.seed(911)

############
###START####
############

genes_results = list()
cinds_clin = c()
cinds_justlncs = c()
cinds_combined = c()

for(j in 1:100){
smp_size <- floor(0.7 * nrow(canc_data))
train_ind <- sample(seq_len(nrow(canc_data)), size = smp_size)
train <- canc_data[train_ind, c(1:(ncol(canc_data)-6), ncol(canc_data))]
test <- canc_data[-train_ind, c(1:(ncol(canc_data)-6), ncol(canc_data))]

###[2.] Feature selection 

###---------------------------------------------------------------
###Extract significant predictors in training data 
###---------------------------------------------------------------

#use training medians to split the test set 
medians = apply(train[,1:(ncol(train)-1)], 2, median)
for(k in 1:length(medians)){
    med = medians[k]
    if(med ==0){
    #if median = 0 then anyone greater than zero is 1 
    l1 = which(train[,k] > 0)
    l2 = which(train[,k] ==0)
    train[l1,k] = 1
    train[l2, k] = 0
    }

    if(!(med ==0)){
    l1 = which(train[,k] >= med)
    l2 = which(train[,k] < med)
    train[l1,k] = 1
    train[l2, k] = 0
    }
}  

source("survival_script_march14_already_labelled_highlow.R")
#means = apply(train[,1:(ncol(train)-1)], 2, mean) #134 with 0 expression in ALL patients 
#zeroes = names(means[which(means ==0)]) #what are they?
#z <- which(colnames(train) %in% zeroes)
#if(!(length(z)==0)){
#  train = train[,-z]
#}
#train[,1:(ncol(train)-1)] = log1p(train[,1:(ncol(train)-1)])
sums = apply(train[,1:(ncol(train)-1)], 2, sum) #134 with 0 expression in ALL patients 
zeroes = names(sums[which(sums <50)]) #if less than 5, then unbalanced 
z <- which(colnames(train) %in% zeroes)
if(!(length(z)==0)){
  train = train[,-z]
}
#means = apply(train[,1:(ncol(train)-1)], 2, mean) #134 with 0 expression in ALL patients 
#zeroes = names(means[which(means <1)]) #what are they?
#z <- which(colnames(train) %in% zeroes)
#if(!(length(z)==0)){
 # train = train[,-z]
#}

train = merge(train, clin, by=c("patient"))
rownames(train) = train$patient
train = train[,-1]
train$status = as.numeric(train$status)

test = merge(test, clin, by=c("patient"))
rownames(test) = test$patient
test = test[,-1]
test$status = as.numeric(test$status)

#surv_test is the loaded survival function which takes in gene 
#name as input and conducts a univariate cox regression model 
genes = as.list(colnames(train)[1:(ncol(train)-9)])
#genes = genes[1:100]
tests_survival2 = llply(genes, surv_test, .progress="text")
tests_survival3 = ldply(tests_survival2, rbind)
tests_survival3$pval = as.numeric(tests_survival3$pval)
#keep sig, how many ? don't want to be more than # of deaths
#so filter by HR and keep the max 
tests_survival3$fdr = p.adjust(tests_survival3$pval, method="fdr")
tests_survival3 = subset(tests_survival3, fdr <=0.05)
if(!(dim(tests_survival3)[1])==0){
length(unique(tests_survival3$gene))
sortme <- function(dt, sort.field) dt[order(-abs(dt[[sort.field]]))]
tests_survival3$HR = as.numeric(tests_survival3$HR)
tests_survival3 = as.data.table(tests_survival3)
sort.field = "HR"
tests_survival3 = sortme(tests_survival3, sort.field)

tests_survival3 = as.data.table(tests_survival3)
tests_survival3$coef = as.numeric(unlist(tests_survival3$coef))
tests_survival3 = tests_survival3[order(-abs(coef))]

z = table(train$status)[2]
tests_survival3 = as.data.frame(tests_survival3)
if(z < dim(tests_survival3)[1]){
tests_survival3 = tests_survival3[1:z, ]
}

#using LASSO 
#with 5fold cross validation to get optimal lamda penalization parameter
z <- which(colnames(train) %in% tests_survival3$gene)
train = train[,c(z, (ncol(train)-8):ncol(train))]

###[1]. just clinical model
clin_cox_model = coxph(Surv(time, status)  ~ grade +  stage + age, data = train)
test = test[,which(colnames(test) %in% colnames(train))]
pred_validation = predict(clin_cox_model, newdata = test)
# Determine concordance
cindex_validation = concordance.index(pred_validation, surv.time = test$time,
                                       surv.event=test$status, method = "noether")

cindex_validation = cindex_validation$c.index
cinds_clin = c(cinds_clin, cindex_validation)

if(!(dim(tests_survival3)[1] <2)){

###[2]. clinical and lncRNA expression model
gene_data = as.data.frame(train[,1:(ncol(train)-9)])
#colnames(gene_data) = colnames(train)[1:(ncol(train)-2)]

#there might be a imbalance between how many patients in each group after dividing patients by the mean 
genes = as.list(colnames(gene_data))
#gene_data = log1p(gene_data)
#gene_data = cbind(gene_data, train[,(ncol(train)-2):ncol(train)])
x <- model.matrix( ~., gene_data)
y <- Surv(train$time, train$status)

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
trainlncs = train[,c(which(colnames(train) %in% c(genes_keep,"time", "status")))]
justlncs = coxph(Surv(time, status)  ~ ., data = trainlncs)
#which ones actually have significant p-value in multivariate model
#keep = names(which(summary(justlncs)$coefficients[,5] <=0.05))
#if(!(length(keep)==0)){
#trainlncs = trainlncs[,c(which(colnames(trainlncs) %in% c(keep,"time", "status")))]
#}

justlncs_cox_model = justlncs

#create survival estimates on validation data
#label test set as 0/1s for each lncRNA based on median expression in the training set 
z = which(names(medians) %in% colnames(trainlncs))
medians = medians[z]
testlncs = test[,which(colnames(test) %in% colnames(trainlncs))]
#Add high/low tags to each gene 
for(k in 1:(ncol(testlncs)-2)){
    med = medians[which(names(medians) == colnames(testlncs)[k])]
    if(med ==0){
    #if median = 0 then anyone greater than zero is 1 
    for(m in 1:nrow(testlncs)){
    genexp <- testlncs[m,k]
    if(genexp > 0){
      testlncs[m,k] <- 1
      }
    if(genexp == 0 ){
      testlncs[m,k] <- 0
      }
    } 
    }
    if(!(med ==0)){
    for(m in 1:nrow(testlncs)){
    genexp <- testlncs[m,k]
    if(genexp >= med){
      testlncs[m,k] <- 1
      }
    if(genexp < med){
      testlncs[m,k] <- 0
      }
    } 
    }
}  

colss = colnames(testlncs)
colsskeep = which(str_detect(colss, "ENSG"))
#test[,colsskeep] = log1p(test[,colsskeep])

pred_validation = predict(justlncs_cox_model, newdata = testlncs)
# Determine concordance
cindex_validation = concordance.index(pred_validation, surv.time = testlncs$time,
                                       surv.event=testlncs$status, method = "noether")

cindex_validation = cindex_validation$c.index
cinds_justlncs = c(cinds_justlncs, cindex_validation)

#combined lncRNAs plus the chosen lncRNAs from above
train = train[,c(which(colnames(train) %in% colnames(trainlncs)), which(colnames(train) %in% c("stage", "age", "grade")))]
combined = coxph(Surv(time, status)  ~ ., data = train)
test = test[,which(colnames(test) %in% colnames(train))]
test = test[,c("age", "stage", "grade")]
testlncs = cbind(test, testlncs)

pred_validation = predict(combined, newdata = testlncs)
# Determine concordance
cindex_validation = concordance.index(pred_validation, surv.time = testlncs$time,
                                       surv.event=testlncs$status, method = "noether")

cindex_validation = cindex_validation$c.index
cinds_combined = c(cinds_combined, cindex_validation)
print(paste("run is", j))
}
}
}

saveRDS(cinds_clin, file="ALL_KIRC_pats_prebin_predictors_clin_cindeces_CV_100_March19nokeep.rds")
saveRDS(cinds_combined, file="ALL_KIRC_prebin_predictors_combined_cindeces_CV_100_March19nokeep.rds")
saveRDS(cinds_justlncs, file="ALL_KIRC_prebin_predictors_justlncs_cindeces_CV_100_March19nokeep.rds")
saveRDS(genes_results, file="ALL_KIRC_pats_prebin_predictors_list_of_sig_genes_CV_100March19nokeep.rds")





































	
	


