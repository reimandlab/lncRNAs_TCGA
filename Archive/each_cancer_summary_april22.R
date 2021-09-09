source("universal_LASSO_survival_script.R")

set.seed(911)

#for each cancer type 
#how many patients in cancer with all the required data, clinical, gene expression...
#how many alive versus dead 
#how many had recurrence 
#how manny detectable lncRNAs in cancer 
#how variable is the lncRNA expression within cancer type 
#how is lncRNA expression associated with age, sex, stage, grade, recurrence, TP53 mutation... 


get_patients = function(dat){
	print(dim(dat))
	print(dat$Cancer[1])
	pats = length(unique(dat$patient))
	
	clin = dat[,c(1,((ncol(dat)-32):ncol(dat)))]
	clin = clin[,c("patient", "type", "gender", "race", "clinical_stage", "histological_type", "histological_grade",
		"vital_status", "tumor_status", "new_tumor_event_type", "treatment_outcome_first_course", "OS", "DSS",
		"DFI", "PFI", "Cancer")]

	rna = dat[,1:(ncol(dat)-34)]
	rna_clin = merge(rna, clin, by="patient")

	surv = as.data.frame(table(dat$OS))
	treat_outcome = as.data.frame(table(dat$treatment_outcome_first_course))
	tum_status = as.data.frame(table(dat$tumor_status))

	print(surv)
	print(treat_outcome)
	print(tum_status)

	return(rna_clin)
}


rna_clin_dats = llply(canc_datas, get_patients

#look at how clinical variables vary between high and low lncRNA expression patients 

#maybe here first at this step should not just label high/low but look which lncRNAs
#are actually expressed 

 lnc_cors = function(dat){
	
	genes = colnames(dat)[2:(ncol(dat)-15)]	
	medians = apply(dat[,2:(ncol(dat)-15)], 2, median)
		for(k in 1:length(medians)){
		    med = medians[k]
		    k = k +1
		    if(med ==0){
		    #if median = 0 then anyone greater than zero is 1 
		    l1 = which(dat[,k] > 0)
		    l2 = which(dat[,k] ==0)
		    dat[l1,k] = 1
		    dat[l2, k] = 0
		    }

		    if(!(med ==0)){
		    l1 = which(dat[,k] >= med)
		    l2 = which(dat[,k] < med)
		    dat[l1,k] = 1
		    dat[l2, k] = 0
		    }
		}  
	
	return(dat)
	}


labelled_lncs = llply(rna_clin_dats, lnc_cors, .progress="text")

#how many lncs are just 0? --> not detectable at all 

get_gene_clin_results = function(dat){

	#remove those genes with all 0s 
	sums = apply(dat[,2:(ncol(dat)-15)], 2, sum)
	z = which(colnames(dat) %in% names(sums[which(sums ==0)]))
	dat = dat[,-z]
	genes = colnames(dat)[2:(ncol(dat)-15)]
	gene_clin_asoc = function(gene){
		newdat = dat[,c(1, which(colnames(dat) %in% gene), (ncol(dat)-14):ncol(dat))]
	}

}


