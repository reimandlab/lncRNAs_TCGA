set.seed(911)

#setWD
setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ")

#load all libraries and functions 
source("check_lnc_exp_cancers.R")

#get candidates files

#------DATA---------------------------------------------------------

#loaded from check_lnc_exp_cancers.R
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
#z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
#ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,8]))
ucsc <- ucsc[-z,]
colnames(ucsc)[8] = "HGNC.symbol"

#------FEATURES-----------------------------------------------------

#cands -- ion channels 
cands = read.csv("ION_CHANNELS_targets_and_families.csv")
cands = merge(cands, ucsc, by = "HGNC.symbol")

#-----------------------------------------------------------------------------------------------------------------------------------
#1 - get cancer data 
#-----------------------------------------------------------------------------------------------------------------------------------


#get dataset for each cancer

cancers = unique(all$type)

get_canc = function(canc){
	z = which(all$type == canc)
	canc_dat = all[z,]
	print('cool')
	return(canc_dat)
}

cancer_data  = llply(cancers, get_canc)

get_canc_data = function(dtt){
  #get cancer specific candidates 
  z = which(colnames(dtt) %in% c(as.character(cands$hg19.ensGene.name2), "age_at_initial_pathologic_diagnosis.y", 
    "OS.time", "OS", "gender.y", "race.y", "patient", "clinical_stage.y", "histological_grade.y", "treatment_outcome_first_course.y", 
    "new_tumor_event_type.y", "Cancer"))
  dtt = dtt[,z]
  #check that all cands have exp >100 in at least 15 patients 
  z = which(str_detect(colnames(dtt), "ENSG"))
  canc = dtt$Cancer[1]
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
	colnames(all_runs_results) = c("IC", "Cancer", "Cindex", "type", "round")

	for(j in 1:100){
		
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
		    l1 = which(train[,z[k]] > 0)
		    l2 = which(train[,z[k]] ==0)
		    train[l1,z[k]] = 1
		    train[l2,z[k]] = 0
		    }

		    if(!(med ==0)){
		    l1 = which(train[,z[k]] >= med)
		    l2 = which(train[,z[k]] < med)
		    train[l1,z[k]] = 1
		    train[l2,z[k]] = 0
		    }
			} #end for k in 1:length(medians) 

		#survival script
		source("survival_script_march14_already_labelled_highlow.R")
		z = which(str_detect(colnames(train), "ENSG"))
		genes = colnames(train)[z]

		#Add high/low tags to each gene 
			z = which(str_detect(colnames(test), "ENSG"))
			for(k in 1:length(z)){
			    print(k)
			    med = medians[which(names(medians) == colnames(test)[z[k]])]
			    if(med ==0){
			    #if median = 0 then anyone greater than zero is 1 
			    for(m in 1:nrow(test)){
			    genexp <- test[m,z[k]]
			    if(genexp > 0){
			      test[m,z[k]] <- 1
			      }
			    if(genexp == 0 ){
			      test[m,z[k]] <- 0
			      }
			    } 
			    }
			    if(!(med ==0)){
			    for(m in 1:nrow(test)){
			    genexp <- test[m,z[k]]
			    if(genexp >= med){
			      test[m,z[k]] <- 1
			      }
			    if(genexp < med){
			      test[m,z[k]] <- 0
			      }
			    } 
			    }
				} #end loop assigning high/low to test set 


		#--------------------------------------------------------------------------------------
		#TRAINING using just lncRNAs 
		#--------------------------------------------------------------------------------------

		#train and test model on each lncRNA candidate  
		train_test_lnc = function(ic){
			#TRAIN 
			trainlncs = train[,c(which(colnames(train) %in% c(ic,"OS", "OS.time")))]
			print(ic)
			colnames(trainlncs)[3] = "IonChannel"
			trainlncs$OS = as.numeric(trainlncs$OS)
			trainlncs$OS.time = as.numeric(trainlncs$OS.time)
			#model trained on 70% of the data for just one lncRNA 
			justlncs = coxph(Surv(OS.time, OS)  ~ IonChannel, data = trainlncs)
	
			#TEST
			testlncs = test[,c(which(colnames(test) %in% c(ic,"OS", "OS.time")))]
			colnames(testlncs)[3] = "IonChannel"
			testlncs$OS = as.numeric(testlncs$OS)   
			testlncs$OS.time = as.numeric(testlncs$OS.time)
			pred_validation = predict(justlncs, newdata = testlncs)
			# Determine concordance
			cindex_validation = concordance.index(pred_validation, surv.time = testlncs$OS.time,
			                                       surv.event=testlncs$OS, method = "noether")

			cindex_validation = cindex_validation$c.index
			row = c(ic, train$Cancer[1],cindex_validation)
			return(row)
		} #end train_test_lnc function 

		lnc_cindices = llply(genes, train_test_lnc, .progress="text")
		lnc_cindices = do.call(rbind.data.frame, lnc_cindices)
		lnc_cindices$type = "IConly"
		colnames(lnc_cindices) = c("IC", "Cancer", "Cindex", "type")

		#--------------------------------------------------------------------------------------
		#TRAINING using just clinical variables 
		#--------------------------------------------------------------------------------------
		#train set clinical 
		z = which(str_detect(colnames(train), "ENSG"))
		clin_train = train[,-z]
		clin_train = clin_train[,which(colnames(clin_train) %in% c("age_at_initial_pathologic_diagnosis.y", "gender.y", "race.y", "clinical_stage.y", 
			"histological_grade.y", "OS", "OS.time"))]
		clin_train$age_at_initial_pathologic_diagnosis.y = as.numeric(clin_train$age_at_initial_pathologic_diagnosis.y)
		#test set clinical 
		test_clin = test[,which(colnames(test) %in% colnames(clin_train))]
		test_clin$age_at_initial_pathologic_diagnosis.y = as.numeric(test_clin$age_at_initial_pathologic_diagnosis.y)

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
			for(l in 3:(ncol(clin_train))){
				if(!(colnames(clin_train)[l] == "age_at_initial_pathologic_diagnosis.y")){
					clin_train[,l] = as.factor(clin_train[,l])
				}
			}

			for(l in 3:(ncol(test_clin))){
				if(!(colnames(test_clin)[l] == "age_at_initial_pathologic_diagnosis.y")){
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
			test_clin$age_at_initial_pathologic_diagnosis.y = as.numeric(test_clin$age_at_initial_pathologic_diagnosis.y)
			pred_validation = predict(justclin, newdata = test_clin)
			# Determine concordance
			cindex_validation = concordance.index(pred_validation, surv.time = test_clin$OS.time,
			                                       surv.event=test_clin$OS, method = "noether", na.rm=TRUE)
			cindex_validation = cindex_validation$c.index
			#which clinical variables were used
			variables_used = colnames(test_clin)
			variables_used = paste(variables_used, collapse="--")
			clinical_cindices = c(variables_used, train$Cancer[1], cindex_validation, "ClinicalVariables")
			names(clinical_cindices) = c("ClinVars", "Cancer", "Cindex", "type")

			clin_train$patient = rownames(clin_train)
			test_clin$patient = rownames(test_clin)

		#--------------------------------------------------------------------------------------
		#TRAINING using lncRNAs and clinical variables 
		#--------------------------------------------------------------------------------------

		variables = colnames(clin_train)
		train$patient = rownames(train)
		test$patient = rownames(test)

     	#train and test model on each lncRNA candidate  
		train_test_lnc_plus_clin = function(ic){
			#TRAIN 
			trainlncs = train[,c(which(colnames(train) %in% c(ic, "patient")))]
			trainlncs = merge(trainlncs, clin_train, by = "patient")
			print(ic)
			print(colnames(trainlncs))
			rownames(trainlncs) = trainlncs$patient
			trainlncs$patient = NULL
			colnames(trainlncs)[1] = "IonChannel"
			trainlncs$OS = as.numeric(trainlncs$OS)
			trainlncs$OS.time = as.numeric(trainlncs$OS.time)
			trainlncs$age_at_initial_pathologic_diagnosis.y = as.numeric(trainlncs$age_at_initial_pathologic_diagnosis.y)
			#model trained on 70% of the data for one lncRNA + clinical data 
			lnc_clin = coxph(Surv(OS.time, OS)  ~., data = trainlncs)
	
			#TEST
			testlncs = test[,c(which(colnames(test) %in% c(ic, "patient")))]
			testlncs = merge(testlncs, test_clin, by = "patient")
			rownames(testlncs) = testlncs$patient
			testlncs$patient = NULL
			colnames(testlncs)[1] = "IonChannel"
			testlncs$OS = as.numeric(testlncs$OS)   
			testlncs$OS.time = as.numeric(testlncs$OS.time)
			testlncs$age_at_initial_pathologic_diagnosis.y = as.numeric(testlncs$age_at_initial_pathologic_diagnosis.y)

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
			row = c(ic, train$Cancer[1],cindex_validation)
			return(row)
		} #end train_test_lnc function 

		ic_clin_cindices = llply(genes, train_test_lnc_plus_clin, .progress="text")
		ic_clin_cindices = do.call(rbind.data.frame, ic_clin_cindices)
		ic_clin_cindices$type = "IonChannel&clin"
		colnames(ic_clin_cindices) = c("IonChannel", "Cancer", "Cindex", "type")

		print(ic_clin_cindices)
		print(clinical_cindices)
		print(lnc_cindices)

		names(ic_clin_cindices) = c("IC", "Cancer", "Cindex", "type")
		names(clinical_cindices) = c("IC", "Cancer", "Cindex", "type")
		names(lnc_cindices) = c("IC", "Cancer", "Cindex", "type")

		all_round = (rbind(ic_clin_cindices, clinical_cindices, lnc_cindices))
		all_round$round = j 
		all_runs_results = rbind(all_runs_results, all_round)

	}#end loop j in 1:100

	all_runs_results = all_runs_results[-1,]
	return(all_runs_results)

}#end function run_cv


all_cancers_results = llply(setup_data, run_cv, .progress="text")

saveRDS(all_cancers_results, file="ION_CHANNELS_28_Cancers_100_internal_CVs_individual_cands_Dec13_v2.rds")

print("DONE and SAVED")


#get average and median c-index for each IonChannel clincial variabels per cancer type

#r = readRDS("ION_CHANNELS_28_Cancers_100_internal_CVs_individual_cands_Dec13.rds")
#r = ldply(r)
#r = as.data.table(r)
#z = which(r$type == "IonChannel&clin")
#r = r[-z,]
#list_r = split(r, by = "Cancer")

#get_sum = function(canc){
	#get summ for c-indices for each lnc and clin
#	colnames(canc)[3] = "cindex"
#	z = which(is.na(canc$cindex))
#	if(!(length(z)==0)){
#	canc = canc[-z,]}
#	canc$cindex = as.numeric(canc$cindex)
#	lncs <- group_by(canc, IonChannel)
#	s = as.data.table(summarise(lncs, mean_cindex = mean(cindex), median_cindex = median(cindex)))
#	clin = s[nrow(s),]
#	s$imp = s$median_cindex - clin$median_cindex
#	return(s)
#}

#lists = llply(list_r, get_sum)
#lists = ldply(lists)
#lists = as.data.table(lists)
#z= which(str_detect(lists$IonChannel, "ENSG"))
#lists = lists[z,]





