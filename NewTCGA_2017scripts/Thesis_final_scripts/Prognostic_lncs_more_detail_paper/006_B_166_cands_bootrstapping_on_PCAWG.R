source("source_script_006_166_cands_bootrstapping_on_PCAWG.R")

#PCAWG datafile --> pcawg_data_all
#TCGA datafile --> tcga_data_setup

#For now just looking at lncRNAs 
#Need to add clinical variables later... 
#well asap 

#-----------------------------------------------------------------------------------------------------------------------------------
#3 - train survival models using all of TCGA data for one lncRNA at a time and test on samples of PCAWG data via bootstrapping 
#-----------------------------------------------------------------------------------------------------------------------------------

library(boot)

#testing
#dtt = tcga_data_setup[[3]]

run_cv = function(dtt){
	
	set.seed(101) 

	all_runs_results = as.data.frame(matrix(ncol=5))
	colnames(all_runs_results) = c("lncRNA", "Cancer", "C-index", "type", "round")

	#Can just train model outside of the loop 

		#train set ------------------------------
		train <- dtt #ALL of TGCA cohort 
		train$OS.time = as.numeric(train$OS.time)
		z = which(is.na(train$OS.time)) 
		if(!(length(z)==0)){
		train = train[-z,]} 

		#survival script
		source("survival_script_march14_already_labelled_highlow.R")
		z = which(str_detect(colnames(train), "ENSG"))
		genes = colnames(train)[z]

		#test set --------------------------------
	    #1. Get PCAWG data
		z = which(unlist(cancers_tests) %in% dtt$Cancer[[1]])
		test = pcawg_data_all[[z]]

		#test set -------------------------------
		colnames(test)[colnames(test)=="status"] = "OS"
		colnames(test)[colnames(test)=="time"] = "OS.time"
		test$OS = as.numeric(test$OS)
		test$OS.time = as.numeric(test$OS.time)
		z = which(is.na(test$OS.time))
		if(!(length(z)==0)){
		test = test[-z,]} 

		#make sure all lncRNAs are present 
		z = which(str_detect(colnames(test), "ENSG"))
		genes = colnames(test)[z]
		#keep all clinical columns
		z1 = which(!(str_detect(colnames(train), "ENSG")))
		z2 = which(colnames(train) %in% genes)
		train = train[,c(z2,z1)]

		set.seed(103) 

		#set-up indices for boostrapping
		samp_size <- 1 * (nrow(test))
		iter <- 100

		# here are 15 blocks of 5 numbers, which will index rows of your matrix x
		samp_mat <- matrix(sample(1:nrow(test), samp_size*iter, replace=T),
                   ncol=samp_size, byrow=T)

		#set-up clinical data --------------------
		#train set clinical 
		z = which(str_detect(colnames(train), "ENSG"))
		clin_train = train[,-z]
		clin_train = clin_train[,which(colnames(clin_train) %in% c("age_at_initial_pathologic_diagnosis", "gender", "race", "clinical_stage", 
			"histological_grade", "OS", "OS.time"))]
		clin_train$age_at_initial_pathologic_diagnosis = as.numeric(clin_train$age_at_initial_pathologic_diagnosis)
		
		#remove columns with less than 2 contrasts 
		check_contrasts = function(col){
			check = dim(table(col))
				if(check >1){
					return("keep")
				}
			}
		
		keep_train = unlist(apply(clin_train, 2, check_contrasts))
		keep = keep_train
		
		z = which(names(keep) %in% "OS.time")
		if(length(z)==0){
			keep = c(keep, "OS.time")
			names(keep)[length(keep)] = "OS.time"
		}

		z = which(colnames(clin_train) %in% names(keep))
		clin_train = clin_train[,z]


	for(j in 1:100){
		
		print(j)
		
		#Bootstrapping PCAWG data

		#2. Sample with replacement of PCAWG data
		#70%? 
		test_set = test[samp_mat[j,],]

		#--------------------------------------------------------------------------------------
		#TRAINING using just lncRNAs 
		#--------------------------------------------------------------------------------------

		#train and test model on each lncRNA candidate  
		train_test_lnc = function(lnc){
			#TRAIN 
			trainlncs = train[,c(which(colnames(train) %in% c(lnc,"OS", "OS.time")))]
			print(lnc)
			colnames(trainlncs)[1] = "lncRNA"
			trainlncs$OS = as.numeric(trainlncs$OS)
			trainlncs$OS.time = as.numeric(trainlncs$OS.time)
			#model trained on 70% of the data for just one lncRNA 
			justlncs = coxph(Surv(OS.time, OS)  ~ lncRNA, data = trainlncs)
			tidy(justlncs)
			#TEST
			testlncs = test_set[,c(which(colnames(test_set) %in% c(lnc,"OS", "OS.time")))]
			colnames(testlncs)[1] = "lncRNA"
			testlncs$OS = as.numeric(testlncs$OS)   
			testlncs$OS.time = as.numeric(testlncs$OS.time)
			pred_validation = predict(justlncs, newdata = testlncs)
			# Determine concordance
			cindex_validation = concordance.index(pred_validation, surv.time = testlncs$OS.time,
			                                       surv.event=testlncs$OS, method = "noether")

			cindex_validation = cindex_validation$c.index
			row = c(lnc, train$Cancer[1],cindex_validation)
			return(row)

		} #end train_test_lnc function 

		lnc_cindices = llply(genes, train_test_lnc, .progress="text")
		lnc_cindices = do.call(rbind.data.frame, lnc_cindices)
		lnc_cindices$type = "lncRNAonly"
		colnames(lnc_cindices) = c("lncRNA", "Cancer", "C-index", "type")
		print(lnc_cindices)


		#--------------------------------------------------------------------------------------
		#TRAINING using just clinical variables 
		#--------------------------------------------------------------------------------------
		#test set clinical 
		colnames(test_set)[colnames(test_set)=="sex"] = "gender"
		colnames(test_set)[colnames(test_set)=="donor_age_at_diagnosis"] = "age_at_initial_pathologic_diagnosis"

		test_clin = test_set[,which(colnames(test_set) %in% colnames(clin_train))]
		test_clin$age_at_initial_pathologic_diagnosis = as.numeric(test_clin$age_at_initial_pathologic_diagnosis)
		keep_test = unlist(apply(test_clin, 2, check_contrasts))
		z = which(colnames(test_clin) %in% c(names(keep_test), "OS"))
		test_clin = test_clin[,z]

			clin_train_round = clin_train[,which(colnames(clin_train) %in% colnames(test_clin))]
			if(length(which(colnames(clin_train_round) %in% "gender"))==1){
			clin_train_round$gender[clin_train_round$gender == "MALE"] = "male"
			clin_train_round$gender[clin_train_round$gender == "FEMALE"] = "female"
			order = c("female", "male")
			clin_train_round$gender <- factor(clin_train_round$gender, levels = order)
			test_clin$gender <- factor(test_clin$gender, levels = order)
			}

			#make sure there aren't any columns with just one value
      		rm = which(sapply(clin_train_round, function(x) { length(unique(x)) }) == 1)
      		if(!(length(rm))==0){
        		clin_train_round = clin_train_round[,-rm]
     		}
			clin_train_round$OS = as.numeric(clin_train_round$OS)
			clin_train_round$OS.time = as.numeric(clin_train_round$OS.time)

			justclin = coxph(Surv(OS.time, OS)  ~ ., data = clin_train_round)
			#test_set 
			test_clin$OS = as.numeric(test_clin$OS)
			test_clin$OS.time = as.numeric(test_clin$OS.time)
			test_clin$age_at_initial_pathologic_diagnosis = as.numeric(test_clin$age_at_initial_pathologic_diagnosis)
			pred_validation = predict(justclin, newdata = test_clin)
			# Determine concordance
			cindex_validation = concordance.index(pred_validation, surv.time = test_clin$OS.time,
			                                       surv.event=test_clin$OS, method = "noether", na.rm=TRUE)
			cindex_validation = cindex_validation$c.index
			#which clinical variables were used
			variables_used = colnames(test_clin)
			variables_used = paste(variables_used, collapse="--")
			clinical_cindices = c(variables_used, train$Cancer[1], cindex_validation, "ClinicalVariables")
			names(clinical_cindices) = c("clinVars", "Cancer", "C-index", "type")


		#--------------------------------------------------------------------------------------
		#TRAINING using lncRNAs and clinical variables 
		#--------------------------------------------------------------------------------------

		variables = colnames(clin_train_round)
		train$patient = rownames(train)
		test_set$patient = rownames(test_set)
		clin_train_round$patient = rownames(clin_train_round)
		test_clin$patient = rownames(test_clin)

     	#train and test model on each lncRNA candidate  
		train_test_lnc_plus_clin = function(lnc){
			#TRAIN 
			trainlncs = train[,c(which(colnames(train) %in% c(lnc, "patient")))]
			trainlncs = merge(trainlncs, clin_train_round, by = "patient")
			print(lnc)
			print(colnames(trainlncs))
			rownames(trainlncs) = trainlncs$patient
			trainlncs$patient = NULL
			colnames(trainlncs)[1] = "lncRNA"
			trainlncs$OS = as.numeric(trainlncs$OS)
			trainlncs$OS.time = as.numeric(trainlncs$OS.time)
			trainlncs$age_at_initial_pathologic_diagnosis = as.numeric(trainlncs$age_at_initial_pathologic_diagnosis)
			#model trained on 70% of the data for one lncRNA + clinical data 
			lnc_clin = coxph(Surv(OS.time, OS)  ~., data = trainlncs)
	
			#TEST
			testlncs = test_set[,c(which(colnames(test_set) %in% c(lnc, "patient")))]
			testlncs = merge(testlncs, test_clin, by = "patient")
			rownames(testlncs) = testlncs$patient
			testlncs$patient = NULL
			colnames(testlncs)[1] = "lncRNA"
			testlncs$OS = as.numeric(testlncs$OS)   
			testlncs$OS.time = as.numeric(testlncs$OS.time)
			testlncs$age_at_initial_pathologic_diagnosis = as.numeric(testlncs$age_at_initial_pathologic_diagnosis)

			pred_validation = predict(lnc_clin, newdata = testlncs)
			z = which(is.na(pred_validation))
			if(!length(z)==0){
			pred_validation = pred_validation[-z]}
			pats = names(pred_validation)
			testlncs = testlncs[which(rownames(testlncs) %in% pats),]

			# Determine concordance
			cindex_validation = concordance.index(pred_validation, surv.time = testlncs$OS.time,
			                                       surv.event=testlncs$OS, method = "noether")

			cindex_validation = cindex_validation$c.index
			row = c(lnc, train$Cancer[1],cindex_validation)
			return(row)
		} #end train_test_lnc function 

		lnc_clin_cindices = llply(genes, train_test_lnc_plus_clin, .progress="text")
		lnc_clin_cindices = do.call(rbind.data.frame, lnc_clin_cindices)
		lnc_clin_cindices$type = "lncRNA&clin"
		colnames(lnc_clin_cindices) = c("lncRNA", "Cancer", "C-index", "type")

		print(lnc_clin_cindices)
		print(clinical_cindices)
		print(lnc_cindices)

		all_round = (rbind(lnc_clin_cindices, clinical_cindices, lnc_cindices))
		all_round$round = j 
		all_runs_results = rbind(all_runs_results, all_round)	


	}#end loop j in 1:100

	all_runs_results = all_runs_results[-1,]
	return(all_runs_results)

}#end function run_cv


all_cancers_results = llply(setup_data, run_cv, .progress="text")

#for each cancer type --> summarize number of c-indices that were labelled as NA and
#number of patients selected in each iteration 

saveRDS(all_cancers_results, file="lncRNAs_100_external_PCAWG_CVs_individual_cands_june19.rds")





















