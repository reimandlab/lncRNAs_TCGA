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

###canc data
source("lihc_source_canc_dataMar21.R")
canc_data = readRDS("LIHC_tcga_RNA_data_only_detectable_iPCAWG_lncs_mar21.rds")
corlncs = readRDS("LIHC_tcga_RNA_data_only_detectable_iPCAWG_lncs_mar21_mostcorrelated_lncs.rds")
#canc_data = corlncs

## set the seed to make your partition reproductible
set.seed(911)

############
###START####
############

genes_results = list()
cinds_clin = c()
cinds_justlncs = c()
cinds_combined = c()

for(j in 1:1000){
smp_size <- floor(0.8 * nrow(canc_data))
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

###########################################################################################
######DATA---------------------------------------------------------------------------------
###########################################################################################

train = merge(train, clin, by=c("patient"))
rownames(train) = train$patient
train = train[,-1]
train$status = as.numeric(train$status)

test = merge(test, clin, by=c("patient"))
rownames(test) = test$patient
test = test[,-1]
test$status = as.numeric(test$status)


###########################################################################################
######PRE-FILTERING LNCRNAS----------------------------------------------------------------
###########################################################################################

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
tests_survival3 = subset(tests_survival3, pval <=0.05)
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

###########################################################################################
######CLINICAL VARIABLES ONLY--------------------------------------------------------------
###########################################################################################

z <- which(colnames(train) %in% tests_survival3$gene)
train = train[,c(z, (ncol(train)-8):ncol(train))]

clin_cox_model = coxph(Surv(time, status)  ~ grade +  stage + age, data = train)
test = test[,which(colnames(test) %in% colnames(train))]
pred_validation = predict(clin_cox_model, newdata = test)
# Determine concordance
cindex_validation = concordance.index(pred_validation, surv.time = test$time,
                                       surv.event=test$status, method = "noether")

cindex_validation = cindex_validation$c.index
#SAVE C-INDEX
cinds_clin = c(cinds_clin, cindex_validation)

###########################################################################################
######LASSO LNCRNA EXPRESSION--------------------------------------------------------------
###########################################################################################

if(!(dim(tests_survival3)[1] <2)){
gene_data = as.data.frame(train[,1:(ncol(train)-9)])

genes = as.list(colnames(gene_data))
#gene_data = log1p(gene_data)
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
if(!(length(genes_keep)==0)){
genes_results[j] = list(as.list(genes_keep))
trainlncs = train[,c(which(colnames(train) %in% c(genes_keep,"time", "status")))]

###########################################################################################
######LNCRNA VARIABLES ONLY----------------------------------------------------------------
###########################################################################################

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
#SAVE C-INDEX
cinds_justlncs = c(cinds_justlncs, cindex_validation)

###########################################################################################
######LNCRNA AND CLINICAL VARIABLES--------------------------------------------------------
###########################################################################################
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
#SAVE C-INDEX
cinds_combined = c(cinds_combined, cindex_validation)
print(paste("run is", j))
}
}
}
}

saveRDS(cinds_clin, file="8020_LIHC_100CV_justclin_detectable_correlated_lncs_PCAWGtcga_mar21.rds")
saveRDS(cinds_combined, file="8020_LIHC_100CV_combined_detectable_correlated_lncs_PCAWGtcga_mar21.rds")
saveRDS(cinds_justlncs, file="8020_LIHC_100CV_justlncs_detectable_correlated_lncs_PCAWGtcga_mar21.rds")
saveRDS(genes_results, file="8020_LIHC_100CV_SIG_genes_detectable_correlated_lncs_PCAWGtcga_mar21.rds")





































	
	


