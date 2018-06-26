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
#dtt = tcga_data_setup[[1]]

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
		rownames(train) = train$patient
		train$patient = NULL

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

		set.seed(103) 

		#set-up indices for boostrapping
		samp_size <- 0.7 * (nrow(test))
		iter <- 100

		# here are 15 blocks of 5 numbers, which will index rows of your matrix x
		samp_mat <- matrix(sample(1:nrow(test), samp_size*iter, replace=T),
                   ncol=samp_size, byrow=T)


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
			if(!(is.na(cindex_validation))){
				#Visualize KM plot for PCAWG data
				colnames(testlncs)[1] = "gene"
  				testlncs$OS.time = testlncs$OS.time/365
  				fit <- survfit(Surv(OS.time, OS) ~ gene, data = testlncs)
         		 s <- ggsurvplot(
         		 title = paste(lnc, dtt$Cancer[1]),
         		 fit, 
         		 xlab = "Time (Years)", 
         		 surv.median.line = "hv",
         		 font.main = c(16, "bold", "black"),
         		 font.x = c(14, "plain", "black"),
         		 font.y = c(14, "plain", "black"),
         		 font.tickslab = c(14, "plain", "black"),
         		 font.legend = 12,
         		 risk.table.fontsize = 5, 
         		 legend.labs = c("Low Expression", "High Expression"),             # survfit object with calculated statistics.
         		 data = testlncs,      # data used to fit survival curves. 
         		 risk.table = TRUE,       # show risk table.
         		 legend = "right", 
         		 pval = TRUE,             # show p-value of log-rank test.
         		 conf.int = FALSE,        # show confidence intervals for 
                            # point estimaes of survival curves.
         		 xlim = c(0,5),        # present narrower X axis, but not affect
                            # survival estimates.
         		 break.time.by = 1,     # break X axis in time intervals by 500.
         		 #palette = colorRampPalette(mypal)(14), 
         		 palette = mypal[c(4,1)],
         		 ggtheme = theme_bw(), # customize plot and risk table with a theme.
         		 risk.table.y.text.col = T, # colour risk table text annotations.
         		 risk.table.y.text = FALSE # show bars instead of names in text annotations
                            # in legend of risk table
          		)
          		print(s)

			}#end if(!(is.na(cindex_validation)))

			row = c(lnc, train$Cancer[1],cindex_validation)
			return(row)

		} #end train_test_lnc function 

		lnc_cindices = llply(genes, train_test_lnc, .progress="text")
		lnc_cindices = do.call(rbind.data.frame, lnc_cindices)
		lnc_cindices$type = "lncRNAonly"
		colnames(lnc_cindices) = c("lncRNA", "Cancer", "C-index", "type")
		print(lnc_cindices)
		all_round = lnc_cindices
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





















