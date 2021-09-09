set.seed(911)

#load all libraries and functions 
source("check_lnc_exp_cancers.R")
source("survival_script_march14_already_labelled_highlow.R")

#setWD
setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNAs_2019_manuscript")

#get candidates files
#Data--------------------------------------------
allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = filter(allCands, data=="TCGA") #179 unique lncRNA-cancer combos, #166 unique lncRNAs 

#which cancer types are the non-unique lncRNAs from?
allCands$Combo = NULL
allCands = allCands[,c("gene", "coef", "HR", "pval", "cancer", "gene_name")]
allCands = allCands[!duplicated(allCands), ]
cands_dups = unique(allCands$gene[which(duplicated(allCands$gene))])

#-----------------------------------------------------------------------------------------------------------------------------------
#1 - get cancer data 
#-----------------------------------------------------------------------------------------------------------------------------------

cancers = unique(allCands$cancer)

get_canc = function(canc){
	canc_dat = subset(rna, Cancer == canc)
	return(canc_dat)
}	

cancer_data = llply(cancers, get_canc)

get_canc_data = function(dtt){
  #get cancer specific candidates 
  z = which(colnames(dtt) %in% c(as.character(allCands$gene[allCands$cancer == dtt$Cancer[1]]), "age_at_initial_pathologic_diagnosis", 
    "OS.time", "OS", "gender", "race", "patient", "clinical_stage", "histological_grade", "treatment_outcome_first_course", 
    "new_tumor_event_type", "Cancer", "PFI", "PFI.time"))
  dtt = dtt[,z]
  #check that all cands have exp >100 in at least 15 patients 
  z = which(str_detect(colnames(dtt), "ENSG"))
  canc = dtt$Cancer[1]
  exp_cut = 100
  for(y in 1:length(z)){
  	lnc = colnames(dtt)[z[y]]
  	#check 15 patients yes/no
  	res = (get_num_pats(lnc, canc, exp_cut))
  	print(res)
  }
  return(dtt)
}

filtered_data = llply(cancer_data, get_canc_data)

#-----------------------------------------------------------------------------------------------------------------------------------
#2 - set up for cross validation 
#-----------------------------------------------------------------------------------------------------------------------------------

set_up_cv = function(dtt){
	dat = dtt
	genes_results = list()
	cinds_clin = c()
	cinds_justlncs = c()
	cinds_combined = c()
	#set.seed(101) 

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

	return(dat)
}

setup_data = llply(filtered_data, set_up_cv)

#-----------------------------------------------------------------------------------------------------------------------------------
#2 - run cross validation looking at one lncRNA at a time? 
#-----------------------------------------------------------------------------------------------------------------------------------

run_cv = function(dtt){
	
	set.seed(101) 

	all_runs_results = as.data.frame(matrix(ncol=5))
	colnames(all_runs_results) = c("lncRNA", "Cancer", "C-index", "type", "round")

	for(j in 1:1000){
		
		print(j)
		smp_size <- floor(0.7 * nrow(dtt))
		train_ind <- sample(seq_len(nrow(dtt)), size = smp_size)
		train <- dtt[train_ind, ]
		test <- dtt[-train_ind, ]

		#train set ------------------------------
		train$OS.time = as.numeric(train$OS.time)
		z = which(is.na(train$OS.time))
		if(!(length(z)==0)){
		train = train[-z,]} 
		rownames(train) = train$patient
		train$patient = NULL

		#test set -------------------------------

		test$OS.time = as.numeric(test$OS.time)
		z = which(is.na(test$OS.time))
		if(!(length(z)==0)){
		test = test[-z,]} 
		rownames(test) = test$patient
		test$patient = NULL

		#------------------------------------------
		#use training medians to split the test set 
		#------------------------------------------
		z = which(str_detect(colnames(train), "ENSG"))
		
		if(length(z)==1){
			medians = median(train[,z])
			names(medians) = colnames(train)[z]
		}

		if(length(z) > 1){
		medians = apply(train[,z], 2, median)
		}

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
			} #end for k in 1:length(medians) 

		#survival script
		z = which(str_detect(colnames(train), "ENSG"))
		genes = colnames(train)[z]

		#Add high/low tags to each gene 
			z = which(str_detect(colnames(test), "ENSG"))
			for(k in 1:length(z)){
			    med = medians[which(names(medians) == colnames(test)[k])]
			    if(med ==0){
			    #if median = 0 then anyone greater than zero is 1 
			    for(m in 1:nrow(test)){
			    genexp <- test[m,k]
			    if(genexp > 0){
			      test[m,k] <- 1
			      }
			    if(genexp == 0 ){
			      test[m,k] <- 0
			      }
			    } 
			    }
			    if(!(med ==0)){
			    for(m in 1:nrow(test)){
			    genexp <- test[m,k]
			    if(genexp >= med){
			      test[m,k] <- 1
			      }
			    if(genexp < med){
			      test[m,k] <- 0
			      }
			    } 
			    }
				} #end loop assigning high/low to test set 


		#--------------------------------------------------------------------------------------
		#TRAINING using just ALL lncRNAs 
		#--------------------------------------------------------------------------------------

		#train and test model on each lncRNA candidate  
		
			#TRAIN 
			#trainlncs = train[,c(which(colnames(train) %in% c(lnc,"OS", "OS.time")))]
			#print(lnc)
			#colnames(trainlncs)[1] = "lncRNA"
			
			z1 = which(str_detect(colnames(train), "ENSG"))
			z2 = which(colnames(train) %in% c("OS", "OS.time"))
			trainlncs = train[,c(z1,z2)]

			trainlncs$OS = as.numeric(trainlncs$OS)
			trainlncs$OS.time = as.numeric(trainlncs$OS.time)
			#model trained on 70% of the data for just one lncRNA 
			justlncs = coxph(Surv(OS.time, OS)  ~ ., data = trainlncs)
	
			#TEST

			z1 = which(str_detect(colnames(test), "ENSG"))
			z2 = which(colnames(test) %in% c("OS", "OS.time"))
			testlncs = test[,c(z1,z2)]

			#testlncs = test[,c(which(colnames(test) %in% c(lnc,"OS", "OS.time")))]
			#colnames(testlncs)[1] = "lncRNA"
			testlncs$OS = as.numeric(testlncs$OS)   
			testlncs$OS.time = as.numeric(testlncs$OS.time)
			pred_validation = predict(justlncs, newdata = testlncs)
			# Determine concordance
			cindex_validation = concordance.index(pred_validation, surv.time = testlncs$OS.time,
			                                       surv.event=testlncs$OS, method = "noether")

			cindex_validation = cindex_validation$c.index

		lnc_cindices = c("all_lncs", train$Cancer[1], cindex_validation, "lncRNAonly")
		names(lnc_cindices) = c("lncRNA", "Cancer", "C-index", "type")

		#--------------------------------------------------------------------------------------
		#TRAINING using just clinical variables 
		#--------------------------------------------------------------------------------------
		#train set clinical 
		z = which(str_detect(colnames(train), "ENSG"))
		clin_train = train[,-z]
		clin_train = clin_train[,which(colnames(clin_train) %in% c("age_at_initial_pathologic_diagnosis", "gender", "race", "clinical_stage", 
			"histological_grade", "OS", "OS.time"))]
		clin_train$age_at_initial_pathologic_diagnosis = as.numeric(clin_train$age_at_initial_pathologic_diagnosis)
		#test set clinical 
		test_clin = test[,which(colnames(test) %in% colnames(clin_train))]
		test_clin$age_at_initial_pathologic_diagnosis = as.numeric(test_clin$age_at_initial_pathologic_diagnosis)

		#remove columns with less than 2 contrasts 
		check_contrasts = function(col){
			check = dim(table(col))
				if(check >1){
					return("keep")
				}
			}
		
		keep_train = unlist(apply(clin_train, 2, check_contrasts))
		keep_test = unlist(apply(test_clin, 2, check_contrasts))
	
		keep = keep_train[which(keep_train %in% keep_test)]
		
		z = which(names(keep) %in% "OS.time")
		if(length(z)==0){
			keep = c(keep, "OS.time")
			names(keep)[length(keep)] = "OS.time"
		}

		z = which(colnames(clin_train) %in% names(keep))
		clin_train = clin_train[,z]
		z = which(colnames(test_clin) %in% names(keep))
		test_clin = test_clin[,z]

			#make sure the levels are the same 
			for(l in 1:(ncol(clin_train)-2)){
				if(!(colnames(clin_train)[l] == "age_at_initial_pathologic_diagnosis")){
					clin_train[,l] = as.factor(clin_train[,l])
				}
			}

			for(l in 1:(ncol(test_clin)-2)){
				if(!(colnames(test_clin)[l] == "age_at_initial_pathologic_diagnosis")){
					test_clin[,l] <- factor(test_clin[,l], levels=levels(clin_train[,which(colnames(clin_train) %in% colnames(clin_train)[l])]))
				}
			}

			#make sure there aren't any columns with just one value
      		rm = which(sapply(clin_train, function(x) { length(unique(x)) }) == 1)
      		if(!(length(rm))==0){
        		clin_train = clin_train[,-rm]
     		}
			clin_train$OS = as.numeric(clin_train$OS)
			clin_train$OS.time = as.numeric(clin_train$OS.time)

			justclin = coxph(Surv(OS.time, OS)  ~ ., data = clin_train)
			#test 
			test_clin$OS = as.numeric(test_clin$OS)
			test_clin$OS.time = as.numeric(test_clin$OS.time)
			test_clin$age_at_initial_pathologic_diagnosis = as.numeric(test_clin$age_at_initial_pathologic_diagnosis)
			pred_validation = predict(justclin, newdata = test_clin)
			# Determine concordance
			cindex_validation = concordance.index(pred_validation, surv.time = test_clin$OS.time,
			                                       surv.event=test_clin$OS, method = "noether", na.rm=TRUE)
			cindex_validation = cindex_validation$c.index
			#which clinical variables were used
			variables_used = colnames(test_clin)[1:(ncol(test_clin)-2)]
			variables_used = paste(variables_used, collapse="--")
			clinical_cindices = c(variables_used, train$Cancer[1], cindex_validation, "ClinicalVariables")
			names(clinical_cindices) = c("clinVars", "Cancer", "C-index", "type")

			clin_train$patient = rownames(clin_train)
			test_clin$patient = rownames(test_clin)

		#--------------------------------------------------------------------------------------
		#TRAINING using lncRNAs and clinical variables 
		#--------------------------------------------------------------------------------------

		variables = colnames(clin_train)
		train$patient = rownames(train)
		test$patient = rownames(test)

     	#train and test model on each lncRNA candidate  
     	z1=which(str_detect(colnames(train), "ENSG"))
     	z2=which(colnames(train) %in% "patient")

			#TRAIN 
			trainlncs = train[,c(z1,z2)]
			trainlncs = merge(trainlncs, clin_train, by = "patient")
			#print(lnc)
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
			z1=which(str_detect(colnames(test), "ENSG"))
     		z2=which(colnames(test) %in% "patient")
			testlncs = test[,c(z1,z2)]
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
		
		lnc_clin_cindices = c("all_lncs", train$Cancer[1],cindex_validation, "lncRNA&clin")
		names(lnc_clin_cindices) = c("lncRNA", "Cancer", "C-index", "type")

		print(lnc_clin_cindices)
		print(clinical_cindices)
		print(lnc_cindices)

		all_round = as.data.frame((rbind(lnc_clin_cindices, clinical_cindices, lnc_cindices)))
		all_round$round = j 
		all_runs_results = rbind(all_runs_results, all_round)

	}#end loop j in 1:100

	all_runs_results = all_runs_results[-1,]
	return(all_runs_results)

}#end function run_cv


all_cancers_results = llply(setup_data, run_cv, .progress="text")
#saveRDS(all_cancers_results, file="lncRNAs_100_internal_CVs_individual_cands_june19.rds")
saveRDS(all_cancers_results, file="lncRNAs_1000_internal_CVs_individual_cands_multivariate_feb25.rds")

#get average and median c-index for each lncRNA clincial variabels per cancer type

r = readRDS("lncRNAs_1000_internal_CVs_individual_cands_multivariate_feb25.rds")
r = ldply(r)
r = as.data.table(r)
z = which(r$type == "lncRNA&clin")
r = r[-z,]
list_r = split(r, by = "Cancer")

get_sum = function(canc){
	#get summ for c-indices for each lnc and clin
	colnames(canc)[3] = "cindex"
	z = which(is.na(canc$cindex))
	if(!(length(z)==0)){
	canc = canc[-z,]}
	canc$cindex = as.numeric(canc$cindex)
	lncs <- group_by(canc, lncRNA)
	s = as.data.table(summarise(lncs, mean_cindex = mean(cindex), median_cindex = median(cindex)))
	clin = s[nrow(s),]
	s$imp = s$median_cindex - clin$median_cindex
	return(s)
}

lists = llply(list_r, get_sum)
lists = ldply(lists)
lists = as.data.table(lists)
z= which(str_detect(lists$lncRNA, "ENSG"))
lists = lists[z,]





