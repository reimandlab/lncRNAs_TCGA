set.seed(911)

source("/u/kisaev/lncRNAs_TCGA/NewTCGA_2017scripts/Thesis_final_scripts/Prognostic_lncs_more_detail_paper/final_scripts_2019/revisions_2020/load_data.R")
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
	z = which(canc_dat$race %in% c("[Not Available]", "[Not Evaluated]", "[Unknown]"))
	if(!(length(z)==0)){
		canc_dat$race[z] = "NA"
	}
	return(canc_dat)
}	

cancer_data = llply(cancers, get_canc)

get_canc_data = function(dtt){
  #get cancer specific candidates 
  z = which(colnames(dtt) %in% c(as.character(allCands$gene[allCands$cancer == dtt$Cancer[1]]), "age_at_initial_pathologic_diagnosis", 
    "OS.time", "OS", "gender", "race", "patient", "clinical_stage", "histological_grade", "treatment_outcome_first_course", 
    "new_tumor_event_type", "Cancer", "PFI", "PFI.time"))
  dtt = dtt[,..z]
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
#test all combinations of models via risk score get HR, pvalue and c-index
#-----------------------------------------------------------------------------------------------------------------------------------

run_cv = function(dtt){
	
	set.seed(101) 

	all_runs_results = as.data.frame(matrix(ncol=5))
	colnames(all_runs_results) = c("lncRNA", "Cancer", "C-index", "type", "round")

	dtt$OS.time = as.numeric(dtt$OS.time)
	z = which(is.na(dtt$OS.time))
	if(!(length(z)==0)){
		dtt = dtt[-z,]} 

	dtt = as.data.frame(dtt)	
	rownames(dtt) = dtt$patient

	dtt$patient = NULL

	z = which(str_detect(colnames(dtt), "ENSG"))
		
	if(length(z)==1){
			medians = median(dtt[,z])
			names(medians) = colnames(dtt)[z]
	}

	if(length(z) > 1){
		medians = apply(dtt[,z], 2, median)
	}

	for(k in 1:length(medians)){
		    med = medians[k]
		    if(med ==0){
		    #if median = 0 then anyone greater than zero is 1 
		    l1 = which(dtt[,k] > 0)
		    l2 = which(dtt[,k] ==0)
		    dtt[l1,k] = 1
		    dtt[l2,k] = 0
		    }

		    if(!(med ==0)){
		    l1 = which(dtt[,k] >= med)
		    l2 = which(dtt[,k] < med)
		    dtt[l1,k] = 1
		    dtt[l2, k] = 0
		    }
		} #end for k in 1:length(medians) 

	#survival script
	z = which(str_detect(colnames(dtt), "ENSG"))
	genes = colnames(dtt)[z]

	#--------------------------------------------------------------------------------------
	#TRAINING using just ALL lncRNAs 
	#--------------------------------------------------------------------------------------
			
	gene_surv = function(gene){

	z1 = which(str_detect(colnames(dtt), gene))
	z2 = which(colnames(dtt) %in% c("OS", "OS.time"))
	all_z = c(z1,z2)
	dtt_model = dtt[,all_z]

	dtt_model$OS = as.numeric(dtt_model$OS)
	dtt_model$OS.time = as.numeric(dtt_model$OS.time)
	justlncs = coxph(Surv(OS.time, OS)  ~ ., data = dtt_model)
	
	risk_scores = predict(justlncs, newdata = dtt_model, type="risk")

	#use risk scores to fit model
	dtt_model$risk_score = risk_scores
	median= median(dtt_model$risk_score)
	dtt_model$risk[dtt_model$risk_score == max(dtt_model$risk_score)] = 1
	dtt_model$risk[dtt_model$risk_score == min(dtt_model$risk_score)] = 0
	justlncs = coxph(Surv(OS.time, OS)  ~ risk, data = dtt_model)

	# Determine concordance
	cindex_validation = glance(justlncs)$concordance
   	hr = summary(justlncs)$coefficients[2]
   	pval=summary(justlncs)$coefficients[5]
   
   	ci05=summary(justlncs)$conf.int[3]
   	ci95=summary(justlncs)$conf.int[4]

	lnc_cindices = c(gene, dtt$Cancer[1], cindex_validation, "lncRNAonly", hr, pval, ci05, ci95)
	names(lnc_cindices) = c("lncRNA", "Cancer", "Cindex", "type", "HR", "pval", "ci05", "ci95")

	dtt_model$risk = NULL
	dtt_model$risk_score = NULL

	#--------------------------------------------------------------------------------------
	#TRAINING using just clinical variables 
	#--------------------------------------------------------------------------------------
		
	#train set clinical 
	z = which(str_detect(colnames(dtt), "ENSG"))
	clin_model = dtt[,-z]
	z = (which(colnames(clin_model) %in% c("age_at_initial_pathologic_diagnosis", "gender", "race", "clinical_stage", 
			"histological_grade", "OS", "OS.time")))
	clin_model = clin_model[,z]
	clin_model$age_at_initial_pathologic_diagnosis = as.numeric(clin_model$age_at_initial_pathologic_diagnosis)
	
	#remove columns with less than 2 contrasts 
	check_contrasts = function(col){
			check = dim(table(col))
				if(check >1){
					return("keep")
				}
			}
		
	keep = unlist(apply(clin_model, 2, check_contrasts))
		
	z = which(names(keep) %in% "OS.time")
		if(length(z)==0){
			keep = c(keep, "OS.time")
			names(keep)[length(keep)] = "OS.time"
		}

	z = which(colnames(clin_model) %in% names(keep))
	clin_model = clin_model[,z]
			
	#make sure there aren't any columns with just one value
    rm = which(sapply(clin_model, function(x) { length(unique(x)) }) == 1)
     	if(!(length(rm))==0){
        		clin_model = clin_model[,-rm]
     		}
	clin_model$OS = as.numeric(clin_model$OS)
	clin_model$OS.time = as.numeric(clin_model$OS.time)
	clin_model$patient = rownames(clin_model)
	
	#--------------------------------------------------------------------------------------
	#TRAINING using lncRNAs and clinical variables 
	#--------------------------------------------------------------------------------------

	variables = colnames(clin_model)
	dtt_model$patient = rownames(dtt_model)

    z1=which(str_detect(colnames(dtt_model), "ENSG"))
    z2=which(colnames(dtt_model) %in% "patient")
    all_z = c(z1,z2)

	#TRAIN 
	dtt_model = dtt_model[,all_z]
	dtt_model_full = merge(dtt_model, clin_model, by = "patient")
	print(colnames(dtt_model_full))
	rownames(dtt_model_full) = dtt_model_full$patient
	dtt_model_full$patient = NULL
	dtt_model_full$OS = as.numeric(dtt_model_full$OS)
	dtt_model_full$OS.time = as.numeric(dtt_model_full$OS.time)
	dtt_model_full$age_at_initial_pathologic_diagnosis = as.numeric(dtt_model_full$age_at_initial_pathologic_diagnosis)

	#model using all variables
	lnc_clin = coxph(Surv(OS.time, OS)  ~ ., data = dtt_model_full)
	risk_scores = predict(lnc_clin, newdata = dtt_model_full, type="risk")

	#use risk scores to fit model
	dtt_model_full$risk_score = risk_scores
	median= median(dtt_model_full$risk_score)
	dtt_model_full$risk[dtt_model_full$risk_score >=median] = 1
	dtt_model_full$risk[dtt_model_full$risk_score <median] = 0
	lnc_clin = coxph(Surv(OS.time, OS)  ~ risk, data = dtt_model_full)

	# Determine concordance
	cindex_validation = glance(lnc_clin)$concordance
	hr = summary(lnc_clin)$coefficients[2]
	pval=summary(lnc_clin)$coefficients[5]
   
   	ci05=summary(lnc_clin)$conf.int[3]
   	ci95=summary(lnc_clin)$conf.int[4]

	#z = which(is.na(pred_validation))
	#if(!length(z)==0){
	#	pred_validation = pred_validation[-z]}
	
	#pats = names(pred_validation)
	#testlncs = testlncs[which(rownames(testlncs) %in% pats),]
	lnc_clin_cindices = c(gene, dtt$Cancer[1],cindex_validation, "lncRNA&clin", hr, pval, ci05, ci95)
	names(lnc_clin_cindices) = c("lncRNA", "Cancer", "C-index", "type", "HR", "pval", "ci05", "ci95")

	print(lnc_clin_cindices)
	print(lnc_cindices)

	all_round = as.data.frame((rbind(lnc_clin_cindices, lnc_cindices)))
	return(all_round)
	}

	all_genes = as.data.frame(ldply(llply(genes, gene_surv)))
	return(all_genes)

}#end function 


all_cancers_results = llply(setup_data, run_cv, .progress="text")
saveRDS(all_cancers_results, file="lncRNAs_clinical_variables_risk_score_models_univariate.rds")




