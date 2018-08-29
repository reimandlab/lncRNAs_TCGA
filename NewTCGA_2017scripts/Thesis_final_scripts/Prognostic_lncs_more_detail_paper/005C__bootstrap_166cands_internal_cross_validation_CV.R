#4. perform 1000CV survival LASSO on each cancer 

source("universal_LASSO_survival_script.R")

set.seed(911)

setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ")

#get candidates files
#Data--------------------------------------------
allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = filter(allCands, data=="TCGA") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 

#which cancer types are the non-unique lncRNAs from?
allCands$Combo = NULL
allCands = allCands[,c("gene", "coef", "HR", "pval", "cancer", "CAT_geneName")]
allCands = allCands[!duplicated(allCands), ]
cands_dups = unique(allCands$gene[which(duplicated(allCands$gene))])

#UCSC gene info
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
#z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
#ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,8]))
ucsc <- ucsc[-z,]

#fantom 
fantom <- fread("lncs_wENSGids.txt", data.table=F) #6088 lncRNAs 
extract3 <- function(row){
	gene <- as.character(row[[1]])
	ens <- gsub("\\..*","",gene)
	return(ens)
}
fantom[,1] <- apply(fantom[,1:2], 1, extract3)
#remove duplicate gene names (gene names with multiple ensembl ids)
z <- which(duplicated(fantom$CAT_geneName))
z <- which(fantom$CAT_geneName %in% rm)
rm <- fantom$CAT_geneName[z]
fantom <- fantom[-z,]

#-----------------------------------------------------------------------------------------------------------------------------------
#1 - get cancer data 
#-----------------------------------------------------------------------------------------------------------------------------------

z = which(cancers %in% allCands$cancer)
cancer_data = canc_datas[z] #cancers list and canc_datas list should be the same 

get_canc_data = function(dtt){
  #get cancer specific candidates 
  z = which(colnames(dtt) %in% c(as.character(allCands$gene[allCands$cancer == dtt$Cancer[1]]), "age_at_initial_pathologic_diagnosis", 
    "OS.time", "OS", "gender", "race", "patient", "clinical_stage", "histological_grade", "treatment_outcome_first_course", 
    "new_tumor_event_type", "Cancer"))
  dtt = dtt[,z]
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
	
	rownames(dtt) = dtt$patient
	set.seed(101) 

	all_runs_results = as.data.frame(matrix(ncol=9))
	colnames(all_runs_results) = c("lncRNA", "Cancer", "C-index", "wald_pval", "se_concordance", "aic", "hr", "type", "round")

	#set-up indices for boostrapping
	samp_size <- 1 * (nrow(dtt))
	iter <- 100

	# here are 15 blocks of 5 numbers, which will index rows of your matrix x
	samp_mat <- matrix(sample(1:nrow(dtt), samp_size*iter, replace=T),
        ncol=samp_size, byrow=T)

	for(j in 1:100){
		
		#fit Cox-PH model using this set of bootstrap 
		test_set = dtt[samp_mat[j,],]

		#fit model and extract c-index, HR and p-value 
		z = which(str_detect(colnames(test_set), "ENSG"))
		
		if(length(z)==1){
			medians = median(test_set[,z])
			names(medians) = colnames(test_set)[z]
		}

		if(length(z) > 1){
		medians = apply(test_set[,z], 2, median)
		}

		for(k in 1:length(medians)){
		    med = medians[k]
		    if(med ==0){
		    #if median = 0 then anyone greater than zero is 1 
		    l1 = which(test_set[,k] > 0)
		    l2 = which(test_set[,k] ==0)
		    test_set[l1,k] = 1
		    test_set[l2, k] = 0
		    }

		    if(!(med ==0)){
		    l1 = which(test_set[,k] >= med)
		    l2 = which(test_set[,k] < med)
		    test_set[l1,k] = 1
		    test_set[l2, k] = 0
		    }
		} #end for k in 1:length(medians) 

		#survival script
		source("survival_script_march14_already_labelled_highlow.R")
		z = which(str_detect(colnames(test_set), "ENSG"))
		genes = colnames(test_set)[z]

		#--------------------------------------------------------------------------------------
		#TRAINING using just lncRNAs 
		#--------------------------------------------------------------------------------------

		#train and test model on each lncRNA candidate  
		train_test_lnc = function(lnc){
			#TRAIN 
			trainlncs = test_set[,c(which(colnames(test_set) %in% c(lnc,"OS", "OS.time")))]
			print(lnc)
			colnames(trainlncs)[1] = "lncRNA"
			trainlncs$OS = as.numeric(trainlncs$OS)
			trainlncs$OS.time = as.numeric(trainlncs$OS.time)
			#model trained on 70% of the data for just one lncRNA 
			justlncs = coxph(Surv(OS.time, OS)  ~ lncRNA, data = trainlncs)
			g = glance(justlncs)
			cindex_validation = g$concordance
			pval = g$p.value.wald
			se_concordance = g$std.error.concordance
			aic = g$AIC	
			hr = summary(justlncs)$coefficients[2]
			row = c(lnc, test_set$Cancer[1],cindex_validation, pval, se_concordance, aic, hr)
			return(row)
		} #end train_test_lnc function 

		lnc_cindices = llply(genes, train_test_lnc, .progress="text")
		lnc_cindices = do.call(rbind.data.frame, lnc_cindices)
		lnc_cindices$type = "lncRNAonly"
		colnames(lnc_cindices) = c("lncRNA", "Cancer", "C-index", "wald_pval", "se_concordance", "aic", "hr", "type")

		#--------------------------------------------------------------------------------------
		#TRAINING using just clinical variables 
		#--------------------------------------------------------------------------------------
		#train set clinical 
		z = which(str_detect(colnames(test_set), "ENSG"))
		clin_train = test_set[,-z]
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

			#make sure the levels are the same 
			for(l in 1:(ncol(clin_train)-2)){
				if(!(colnames(clin_train)[l] == "age_at_initial_pathologic_diagnosis")){
					clin_train[,l] = as.factor(clin_train[,l])
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
			g = glance(justclin)
			cindex_validation = g$concordance
			pval = g$p.value.wald
			se_concordance = g$std.error.concordance
			aic = g$AIC	
			hr = "na"
#			row = c(lnc, test_set$Cancer[1],cindex_validation, pval, se_concordance, aic, hr)

			#which clinical variables were used
			variables_used = colnames(clin_train)[1:(ncol(clin_train)-2)]
			variables_used = paste(variables_used, collapse="--")
			clinical_cindices = c(variables_used, test_set$Cancer[1], cindex_validation, pval, se_concordance, aic, hr, "ClinicalVariables")
			names(clinical_cindices) = c("clinVars", "Cancer", "C-index", "wald_pval", "se_concordance", "aic", "hr", "type")

			clin_train$patient = rownames(clin_train)

		#--------------------------------------------------------------------------------------
		#TRAINING using lncRNAs and clinical variables 
		#--------------------------------------------------------------------------------------

		variables = colnames(clin_train)

     	#train and test model on each lncRNA candidate  
		train_test_lnc_plus_clin = function(lnc){
			#TRAIN 
			trainlncs = test_set[,c(which(colnames(test_set) %in% c(lnc,"OS", "OS.time")))]
			print(lnc)
			colnames(trainlncs)[1] = "lncRNA"
			trainlncs$OS = as.numeric(trainlncs$OS)
			trainlncs$OS.time = as.numeric(trainlncs$OS.time)
			trainlncs$patient = rownames(trainlncs)
			clin_train$patient = rownames(clin_train)
			trainlncs = merge(trainlncs, clin_train, by=c("patient", "OS", "OS.time"))
			#model trained on 70% of the data for just one lncRNA 
			trainlncs$patient = NULL
			justlncs = coxph(Surv(OS.time, OS)  ~ ., data = trainlncs)
			g = glance(justlncs)
			cindex_validation = g$concordance
			pval = g$p.value.wald
			se_concordance = g$std.error.concordance
			aic = g$AIC	
			hr = summary(justlncs)$coefficients[2]
			row = c(lnc, test_set$Cancer[1],cindex_validation, pval, se_concordance, aic, hr)
			return(row)

		} #end train_test_lnc function 

		lnc_clin_cindices = llply(genes, train_test_lnc_plus_clin, .progress="text")
		lnc_clin_cindices = do.call(rbind.data.frame, lnc_clin_cindices)
		lnc_clin_cindices$type = "combined"
		colnames(lnc_clin_cindices) = c("lncRNA", "Cancer", "C-index", "wald_pval", "se_concordance", "aic", "hr", "type")

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

#------DO NOT RUN-----------------------------------------------------------------------------------------
##all_cancers_results = llply(setup_data, run_cv, .progress="text")

##saveRDS(all_cancers_results, file="lncRNAs_100_bootstrapping_individual_cands_Cox_PH_models_Aug28.rds")


#------SUMMARIZE RESULTS----------------------------------------------------------------------------------

res = ldply(all_cancers_results)

#Within each cancer type plot values
cancers = unique(res$Cancer)

get_sum_boot = function(canc){
	resul = subset(res, Cancer == canc)

	values = colnames(resul)[3:7]
	#summarize each metric

	get_val_plot = function(val){
		z = which(colnames(resul) == val)
		resul[,z] = as.numeric(resul[,z])
		colnames(resul)[z] = "val"
		resul = as.data.table(resul)
		resul = as.data.table(filter(resul, type %in% c("lncRNAonly", "ClinicalVariables")))
		resul = resul[order(val)]
		z = which(str_detect(resul$lncRNA, "ENSG"))
		resul$lncRNA[-z] = "Clinical"
		order = c("Clinical", unique(resul$lncRNA[z]))
		resul$lncRNA = factor(resul$lncRNA, levels=order)
		ggboxplot(resul, x="lncRNA", y="val", color="black", fill="type", palette=c("lightpink3", "lightsteelblue4", "grey")) + 
		stat_compare_means(label = "p.signif", method = "wilcox.test", ref.group = "Clinical") + 
		xlab("Predictor")+ theme_classic() 


	}


}



























