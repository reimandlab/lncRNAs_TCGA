#42candidate_lncs_limma_diffCoexpression .R

#Karina Isaev
#September 5th, 2017

#Purpose: using the top 5 cancers selected for analysis, 
#run co-expression analysis with cancer type and sex? as confounders
#to get list of co-expressed PCGs with NEAT1 to run through m:Explorer  

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#RNA files used here obtained from script: 
#pcawg_analysis_July2017/top5_cancers_extraction_script3.R 
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#Preamble#-------------------------------------------------
options(stringsAsFactors=F)
source("universal_LASSO_survival_script.R")

#Later on - incorporate FANTOM5 and CRISPRi lncRNAs 

#Libraries#------------------------------------------------
library(data.table)
library(survival)
library(ggplot2)
library(ggfortify)
library(cluster)
library(qqman)
library(dplyr)
library(rafalib)
library(RColorBrewer) 
#library(genefilter)
library(gplots) ##Available from CRAN
library(survminer)
library(MASS)
library(Hmisc)
library(gProfileR)
library(wesanderson)
library(ggsci)
library(gridExtra)
library(ggpubr)
library(ggthemes)
library(plyr)

mypal = pal_npg("nrc", alpha = 0.7)(10)


#Data#-------------------------------------------------

#UCSC gene info
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,6]))
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

#split by cancer type 

z = which(str_detect(colnames(pcg), "ENSG"))	
pcg = as.data.frame(pcg)
pcg[,z] <- log1p(pcg[,z])

#2. Get lncRNA - median within each tissue type
tissues <- unique(pcg$type)

#---------------------------------------------------------
#Pre-Processing - set up lnc/PCG matrix for LM
#---------------------------------------------------------
#Function 1
#input: tissue 
#output: list of dataframes by tissue
get_tissue_specific <- function(tissue){
	tis <- pcg[pcg$type==tissue,]
	return(tis)
}
tissues_data <- llply(tissues, get_tissue_specific, .progress="text")

#FUNCTION2 - get PCG data, complete dataset necessary for regression analysis  
getClinical <- function(df){
	#add RCC1 state column

	z = which(colnames(df) %in% c("patient", "ENSG00000180198", "OS", "OS.time", "Cancer"))
	rcc1 = df[,z]
	med = median(as.numeric(rcc1[,2]))
	z = which(rcc1[,2] >= med)
	rcc1$rcc1_tag = ""
	rcc1$rcc1_tag[z] = 1
	rcc1$rcc1_tag[-z] = 0

	order = c(0, 1)
	rcc1$rcc1_tag = as.numeric(rcc1$rcc1_tag)

	rcc1$rcc1_tag <- factor(rcc1$rcc1_tag, levels = order)

	rcc1$OS = as.numeric(rcc1$OS)
 	rcc1$OS.time = as.numeric(rcc1$OS.time)
  
  	lncs = coxph(Surv(OS.time, OS)  ~ rcc1_tag, data = rcc1)
  	row <- c("RCC1", summary(lncs)$coefficients[1,c(1,2,5)],  summary(lncs)$conf.int[1,c(3,4)], rcc1$Cancer[1])
    
    colnames(rcc1)[2] = "RCC1"

  	rcc1$OS.time = rcc1$OS.time/365
  	fit <- survfit(Surv(OS.time, OS) ~ rcc1_tag, data = rcc1)
         
          s <- ggsurvplot(
          title = paste("RCC1", rcc1$Cancer[1], "HR =", round(as.numeric(row[3]), digits=4)),
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
          data = rcc1,      # data used to fit survival curves. 
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
          ggtheme = theme_minimal(), # customize plot and risk table with a theme.
          risk.table.y.text.col = T, # colour risk table text annotations.
          risk.table.y.text = FALSE # show bars instead of names in text annotations
                            # in legend of risk table
          )
          print(s)

	return(row)
}

pdf("RCC1_KM_plots_TCGA_cancers_June14.pdf")
dividedWpcgs <- llply(tissues_data, getClinical, .progress = "text")
dev.off()
print("pass2")

cox_results = do.call(rbind.data.frame, dividedWpcgs)
colnames(cox_results) = c("Predictor", "coef", "HR", "Cox_Pval", "CIlower", "CIhigher", "Cancer")
cox_results$Cox_Pval = as.numeric(cox_results$Cox_Pval)
cox_results = as.data.table(cox_results)
cox_results = cox_results[order(Cox_Pval)]
cox_results$cox_fdr = p.adjust(cox_results$Cox_Pval, method="fdr")

#get in order of decreasing p-value to increasing p-value 

tissues = unique(cox_results$Cancer)

get_tissue_specific <- function(tissue){
	tis <- pcg[pcg$Cancer==tissue,]
	return(tis)
}
tissues_data <- llply(tissues, get_tissue_specific, .progress="text")

#FUNCTION2 - get PCG data, complete dataset necessary for regression analysis  
getClinical <- function(df){
	#add RCC1 state column

	z = which(colnames(df) %in% c("patient", "ENSG00000180198", "OS", "OS.time", "Cancer"))
	rcc1 = df[,z]
	med = median(as.numeric(rcc1[,2]))
	z = which(rcc1[,2] >= med)
	rcc1$rcc1_tag = ""
	rcc1$rcc1_tag[z] = 1
	rcc1$rcc1_tag[-z] = 0

	order = c(0, 1)
	rcc1$rcc1_tag = as.numeric(rcc1$rcc1_tag)

	rcc1$rcc1_tag <- factor(rcc1$rcc1_tag, levels = order)

	rcc1$OS = as.numeric(rcc1$OS)
 	rcc1$OS.time = as.numeric(rcc1$OS.time)
      
    colnames(rcc1)[2] = "RCC1"

  	rcc1$OS.time = rcc1$OS.time/365
  	fit <- survfit(Surv(OS.time, OS) ~ rcc1_tag, data = rcc1)
         
          s <- ggsurvplot(
          title = paste("RCC1", rcc1$Cancer[1], "HR =", round(as.numeric(row[3]), digits=4)),
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
          data = rcc1,      # data used to fit survival curves. 
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
          ggtheme = theme_minimal(), # customize plot and risk table with a theme.
          risk.table.y.text.col = T, # colour risk table text annotations.
          risk.table.y.text = FALSE # show bars instead of names in text annotations
                            # in legend of risk table
          )
          print(s)
}

pdf("RCC1_KM_plots_TCGA_cancers_June14.pdf")
dividedWpcgs <- llply(tissues_data, getClinical, .progress = "text")
dev.off()
print("pass2")
















