#4. perform 1000CV survival LASSO on each cancer 

source("universal_LASSO_survival_script.R")

set.seed(911)

setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ")


#all RNA, PCG and clinical data is lwoded in already? 
#clinical columns attached to rna file **

#Four things we want to do here for now:

#1. check correlation of each lncRNA candidate with available clinical variables for cancer type

#2. calculate c-index distribution for each lncRNA prognostic marker comapred to clinical variables available 

#3. visualize results of each lncRNA candidate via forest plot to show differences in HR and p-value between all predictors 

#4. train model using all lncRNA predictors obtained from TCGA analysis and test performance of model on PCAWG, 50% of PCAWG at a time randomly selected --> train model and get c-index distribution , compare to age? whatever clinical data available? 


#get candidates files
#Data--------------------------------------------
allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_May4.rds")
#save only the ones that came from the noFDR appraoch 
allCands = filter(allCands, AnalysisType == "noFDR", data=="TCGA") #173 unique lncRNA-cancer combos, #166 unique lncRNAs 
#23 unique cancer types 

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
rm <- fantom$CAT_geneName[z]
z <- which(fantom$CAT_geneName %in% rm)
fantom <- fantom[-z,]


#Analysis--------------------------------------------

#Function 1
#input: tissue 
#output: list of dataframes by tissue

tissues <- unique(allCands$cancer)

get_tissue_specific <- function(tissue){
	tis <- rna[rna$Cancer==tissue,]
	return(tis)
}
tissues_data <- llply(tissues, get_tissue_specific, .progress="text")

#Function 2 
#input: cancer specific RNA data
#output: lncRNA specific RNA and clinical data 

#Function 2
#for each cancer check correlation of lncRNA with clinical variables 

get_lnc_canc = function(dat){
	canc = dat$Cancer[1]
	lncs = as.character(unique(subset(allCands, cancer == canc)$gene))
	
	#for each lncRNA candidate, look at relationship wtih clinical variables 
	#subset to lncRNA expression and clinical variables 
	
	get_lnc_cor = function(lnc){
		canc = dat$Cancer[1]
		z = which(colnames(dat) == lnc)
		lnc_data = dat[,c(z, 5788:5820)]
		#add lnc risk group 
		median2 = median(as.numeric(lnc_data[,1]))
		lnc_data$lncRNA_exp = ""
		if(median2 == 0){
			lnc_data$lncRNA_exp[lnc_data[,1] == 0] = "Low"
			lnc_data$lncRNA_exp[lnc_data[,1] > 0] = "High"
		}
		if(!(median2 == 0)){
			lnc_data$lncRNA_exp[lnc_data[,1] < median2] = "Low"
			lnc_data$lncRNA_exp[lnc_data[,1] >= median2] = "High"
		}

		lnc_data$OS.time = as.numeric(lnc_data$OS.time)
		lnc_data$OS = as.numeric(lnc_data$OS)

		#change levels low high
		lnc_data$lncRNA_exp <- factor(lnc_data$lncRNA_exp, levels = c("Low", "High"))

		#risk type 
		justlncs = coxph(Surv(OS.time, OS)  ~ lncRNA_exp, data = lnc_data)
		HR = summary(justlncs)$coefficients[2]

		if(HR > 1){
			risk = "High"
		}

		if(HR < 1){
			risk = "Low"
		}

		lnc_data$risk_type = risk

		#which clinical variables have enough filled in patients so that can look at relationships
		check_contrasts = function(col){
				check = dim(table(lnc_data[,col]))
				if(check >1){
					return(col)
				}
			}
		cols = unlist(llply(colnames(lnc_data), check_contrasts))
		lnc_data = lnc_data[,which(colnames(lnc_data) %in% c(cols, "risk_type"))]


		#first model OS versus PFI - progression free interval - 1 if patient has new tumour event 
		#overall survival 

		lnc_data$OS.time = lnc_data$OS.time/365
  		lnc_data$PFI = as.numeric(lnc_data$PFI)
  		lnc_data$PFI.time = as.numeric(lnc_data$PFI.time)/365

		pfi <- survfit( Surv(PFI.time, PFI) ~ lncRNA_exp, data = lnc_data)
		os <- survfit( Surv(OS.time, OS) ~ lncRNA_exp, data = lnc_data)


		os_pval = summary(coxph(Surv(OS.time, OS)  ~ lncRNA_exp, data = lnc_data))$coefficients[5]
		pfi_pval = summary(coxph(Surv(PFI.time, PFI)  ~ lncRNA_exp, data = lnc_data))$coefficients[5]


		fit <- list(PFI = pfi, OS = os)

		#plot curves :           
          s <- ggsurvplot(
          title = paste(lnc, canc, "\nOS HR =", round(as.numeric(HR), digits=4), 
          	"\nOS CoxPH pval=", as.numeric(os_pval),
          	"\nPFI CoxPH pval=", as.numeric(pfi_pval)), 
          fit, 
          xlab = "Time (Years)", combine = TRUE,
          surv.median.line = "hv",
          font.main = c(16, "bold", "black"),
          font.x = c(14, "plain", "black"),
          font.y = c(14, "plain", "black"),
          font.tickslab = c(14, "plain", "black"),
          font.legend = 12,
          risk.table.fontsize = 5, 
          #legend.labs = c("Low Expression", "High Expression"),             # survfit object with calculated statistics.
          data = lnc_data,      # data used to fit survival curves. 
          risk.table = TRUE,       # show risk table.
          legend = "right", 
          pval = TRUE,             # show p-value of log-rank test.
          conf.int = FALSE,        # show confidence intervals for 
                            # point estimaes of survival curves.
          xlim = c(0,5),        # present narrower X axis, but not affect
                            # survival estimates.
          break.time.by = 1,     # break X axis in time intervals by 500.
          #palette = colorRampPalette(mypal)(14), 
          palette = "jco",
          ggtheme = theme_minimal(), # customize plot and risk table with a theme.
          risk.table.y.text.col = T, # colour risk table text annotations.
          risk.table.y.text = FALSE # show bars instead of names in text annotations
                            # in legend of risk table
          )
          print(s)

          row = c(lnc, canc, os_pval, HR, pfi_pval)
          names(row) = c("lncRNA", "Cancer", "OS_pval", "HR", "PFI_pval")
          return(row)

	}#end get_lnc_cor

	lnc_pfi_results = llply(lncs, get_lnc_cor, .progress="text")
	results_lncs1 = as.data.frame(do.call("rbind", lnc_pfi_results))
}


pdf("TCGA_lncRNA_candidates_OS_PFI_survival_plots.pdf", width=10)
all_canc_lnc_data = llply(tissues_data, get_lnc_canc, .progress="text")
dev.off()










































































































