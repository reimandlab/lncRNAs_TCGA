###---------------------------------------------------------------
###Load libraries and data - April 16th 
###---------------------------------------------------------------

source("source_code_Cox_MonteCarlo_CV_April12.R")
require(caTools)

#start with only lncRNA_intergenic
#lincs = subset(fantom, CAT_geneClass == "lncRNA_intergenic")
#z = which(colnames(rna) %in% lincs$gene)
#rna = as.data.frame(rna)
#rna = rna[,c(z, (ncol(rna)-5):ncol(rna))]

###[2.] Data splitting 

###---------------------------------------------------------------
###PCA using lncRNA expression 
#can then compare how using all genes compared to just using
#the ones chosen by LASSO at the end 
###---------------------------------------------------------------

#check if this person is in my analysis: TCGA-61-2095
library(glmnet)
library(survcomp)
library(caret)
library(stringr)
library(factoextra)

rna = as.data.frame(rna)
dim(rna)
dim(pcg)
dim(norm)
dim(met)

###---------------------------------------------------------------

#function that tests each lncRNA's survival 

#1. remove discrepancy 
z = which(rna$vital_status == "[Discrepancy]")
rna = rna[-z,]

#2. list of cancers to apply function to 
cancers = as.list(unique(rna$Cancer))

#3. function that splits data into cancers 
get_canc = function(canc){
	canc_data = rna[which(rna$Cancer == canc),]
	return(canc_data)
}

canc_datas = llply(cancers, get_canc)

#4. perform 1000CV survival LASSO on each cancer 
set.seed(911)

lasso_cv = function(dat){
	genes_results = list()
	cinds_clin = c()
	cinds_justlncs = c()
	cinds_combined = c()
	set.seed(101) 

	stage = as.data.table(table(dat$clinical_stage))
	rm = stage$V1[which(stage$N <10)]
	z = which(dat$clinical_stage %in% rm)
	if(!(length(z)==0)){
	dat = dat[-z,]}

	grade = as.data.table(table(dat$histological_grade))
	rm = grade$V1[which(grade$N <10)]
	z = which(dat$histological_grade %in% rm)
	if(!(length(z)==0)){
	dat = dat[-z,]}

	for(j in 1:5){

		smp_size <- floor(0.8 * nrow(dat))
		train_ind <- sample(seq_len(nrow(dat)), size = smp_size)
		train <- dat[train_ind, ]
		test <- dat[-train_ind, ]

		train$OS.time = as.numeric(train$OS.time)
		z = which(is.na(train$OS.time))
		if(!(length(z)==0)){
		train = train[-z,]} 

		test$OS.time = as.numeric(test$OS.time)
		z = which(is.na(test$OS.time))
		if(!(length(z)==0)){
		test = test[-z,]} 

		###---------------------------------------------------------------
		###Extract significant predictors in training data 	
		###---------------------------------------------------------------

		#use training medians to split the test set 
		medians = apply(train[,2:(ncol(train)-34)], 2, median)
		for(k in 1:length(medians)){
		    med = medians[k]
		    k = k +1
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

		#survival script
		source("survival_script_march14_already_labelled_highlow.R")
		sums = apply(train[,2:(ncol(train)-34)], 2, sum) #134 with 0 expression in ALL patients 
		zeroes = names(sums[which(sums <10)]) #if less than 50, then unbalanced 
		z <- which(colnames(train) %in% zeroes)
		if(!(length(z)==0)){
		  train = train[,-z]
		}

		genes = as.list(colnames(train)[2:(ncol(train)-34)])
		#genes = genes[1:100]
		tests_survival2 = llply(genes, surv_test, .progress="text")
		print("p1")
		tests_survival3 = ldply(tests_survival2, rbind)
		print("p2")
		tests_survival3$pval = as.numeric(tests_survival3$pval)
		tests_survival3$fdr = p.adjust(tests_survival3$pval, method="fdr")
		tests_survival3 = subset(tests_survival3, fdr <=0.2)
		if(!(dim(tests_survival3)[1])==0){
		length(unique(tests_survival3$gene))
		sortme <- function(dt, sort.field) dt[order(-abs(dt[[sort.field]]))]
		tests_survival3$HR = as.numeric(tests_survival3$HR)
		tests_survival3 = as.data.table(tests_survival3)
		sort.field = "HR"
		tests_survival3 = sortme(tests_survival3, sort.field)

		print("pass")

		tests_survival3 = as.data.table(tests_survival3)
		tests_survival3$coef = as.numeric(unlist(tests_survival3$coef))
		tests_survival3 = tests_survival3[order(-abs(coef))]

		z = table(train$OS)[2]
		tests_survival3 = as.data.frame(tests_survival3)
		if(z < dim(tests_survival3)[1]){
		tests_survival3 = tests_survival3[1:z, ]
		}

		print("pass2")

		z <- which(colnames(train) %in% tests_survival3$gene)
		train = train[,c(z, 1, (ncol(train)-33):ncol(train))]

		if(!(dim(tests_survival3)[1] <2)){
			#-----------------------------------
			#TRAINING using just lncRNAs 
			#-----------------------------------
			gene_data = as.data.frame(train[,1:(ncol(train)-35)])
			genes = as.list(colnames(gene_data))
			x <- model.matrix( ~., gene_data)
			train$OS = as.numeric(train$OS)
			y <- Surv(train$OS.time, train$OS)
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
			trainlncs = train[,c(which(colnames(train) %in% c(genes_keep,"OS", "OS.time")))]
			justlncs = coxph(Surv(OS.time, OS)  ~ ., data = trainlncs)
			justlncs_cox_model = justlncs
			
			print("pass3")

			#-----------------------------------
			#TESTING using just lncRNAs 
			#-----------------------------------
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
			} #end loop assigning high/low to test set 
			
			print("pass4")

			colss = colnames(testlncs)
			colsskeep = which(str_detect(colss, "ENSG"))
			testlncs$OS = as.numeric(testlncs$OS)    
			pred_validation = predict(justlncs_cox_model, newdata = testlncs)
			# Determine concordance
			cindex_validation = concordance.index(pred_validation, surv.time = testlncs$OS.time,
			                                       surv.event=testlncs$OS, method = "noether")

			cindex_validation = cindex_validation$c.index
			cinds_justlncs = c(cinds_justlncs, cindex_validation)

			#-------------end lncRNAs as predictors-------------------------------------------
			#-----------------------------------
			#TRAINING using just clinical variables  
			#-----------------------------------
			clin_train = as.data.frame(train[,(ncol(train)-34):ncol(train)])
			rownames(clin_train) = clin_train$patient
			clin_train = clin_train[,c(4:6, 10)]
			clin_train$age_at_initial_pathologic_diagnosis = as.numeric(clin_train$age_at_initial_pathologic_diagnosis)
			test_clin = test[,which(colnames(test) %in% colnames(clin_train))]

			#remove columns with less than 2 contrasts 
			check_contrasts = function(col){
				check = dim(table(col))
				if(check >1){
					return("keep")
				}
			}
			keep_train = unlist(apply(clin_train, 2, check_contrasts))
			keep_test = unlist(apply(test_clin, 2, check_contrasts))
			keep = keep_test[which(keep_train %in% keep_test)]
			z = which(colnames(clin_train) %in% names(keep))
			clin_train = clin_train[,z]
			trainclin = train[,c(which(colnames(train) %in% c(colnames(clin_train),"OS", "OS.time")))]
			trainclin$age_at_initial_pathologic_diagnosis = as.numeric(trainclin$age_at_initial_pathologic_diagnosis)
			justclin = coxph(Surv(OS.time, OS)  ~ ., data = trainclin)
			testclin = test[,which(colnames(test) %in% colnames(trainclin))]
			testclin$OS = as.numeric(testclin$OS)
			testclin$OS.time = as.numeric(testclin$OS.time)
			testclin$age_at_initial_pathologic_diagnosis = as.numeric(testclin$age_at_initial_pathologic_diagnosis)
			pred_validation = predict(justclin, newdata = testclin)
			# Determine concordance
			cindex_validation = concordance.index(pred_validation, surv.time = testlncs$OS.time,
			                                       surv.event=testlncs$OS, method = "noether")
			cindex_validation = cindex_validation$c.index
			cinds_clin = c(cinds_clin, cindex_validation)

			#-------------end clinical variables as predictors------------------------------------
			#-----------------------------------
			#TRAINING using lncRNAs and clinical variables  
			#-----------------------------------

			variables = unique(c(colnames(trainclin), colnames(trainlncs)))
			trainboth = train[,which(colnames(train) %in% variables)]
			trainboth$age_at_initial_pathologic_diagnosis = as.numeric(trainboth$age_at_initial_pathologic_diagnosis)
			both = coxph(Surv(OS.time, OS)  ~ ., data = trainboth)
			testboth = test[,which(colnames(test) %in% colnames(trainboth))]
			testboth$OS = as.numeric(testboth$OS)
			testboth$OS.time = as.numeric(testboth$OS.time)
			testboth$age_at_initial_pathologic_diagnosis = as.numeric(testboth$age_at_initial_pathologic_diagnosis)
			pred_validation = predict(both, newdata = testboth)
			# Determine concordance
			cindex_validation = concordance.index(pred_validation, surv.time = testboth$OS.time,
			                                       surv.event=testboth$OS, method = "noether")
			cindex_validation = cindex_validation$c.index
			cinds_combined = c(cinds_combined, cindex_validation)

		}#dim(pvalues)[1] >2
	}

	}#end 1:1000 loop
	return(cinds_combined)
	return(cinds_clin)
	return(cinds_justlncs)
	return(genes_results)
}

canc_datas_survival_results = llply(canc_datas, lasso_cv)


























