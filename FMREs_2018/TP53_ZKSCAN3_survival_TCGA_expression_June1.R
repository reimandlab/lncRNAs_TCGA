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
library(patchwork)

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
	
  #---------------------------------------------------------
  #ZSCKAN3 + TP53 
  #---------------------------------------------------------

	z = which(colnames(df) %in% c("patient", "ENSG00000189298", "ENSG00000141510", "OS", "OS.time", "Cancer"))
	ZKSCAN3 = df[,z]
	med = median(as.numeric(ZKSCAN3[,2]))
	z = which(ZKSCAN3[,2] >= med)
	ZKSCAN3$ZKSCAN3_tag = ""
	ZKSCAN3$ZKSCAN3_tag[z] = 1
	ZKSCAN3$ZKSCAN3_tag[-z] = 0

  #TP53
  med = median(as.numeric(ZKSCAN3[,3]))
  z = which(ZKSCAN3[,3] >= med)
  ZKSCAN3$TP53_tag = ""
  ZKSCAN3$TP53_tag[z] = 1
  ZKSCAN3$TP53_tag[-z] = 0

  order = c(0, 1)
	ZKSCAN3$ZKSCAN3_tag = as.numeric(ZKSCAN3$ZKSCAN3_tag)
	ZKSCAN3$ZKSCAN3_tag <- factor(ZKSCAN3$ZKSCAN3_tag, levels = order)

  #TP53
  ZKSCAN3$TP53_tag = as.numeric(ZKSCAN3$TP53_tag)
  ZKSCAN3$TP53_tag <- factor(ZKSCAN3$TP53_tag, levels = order)

	#New categorial variable --> TP53/Z high high, high low, low high, low, low 
  ZKSCAN3$combo = ""
  get_combo = function(r){
    z = as.numeric(as.character(r[[7]]))
    t = as.numeric(as.character(r[[8]]))
    if((z == 1) & (t == 1)){
      new = "High_Both"
    }
    if((z==1) & (t==0)){
      new = "HighZ_LowT"
    } 
    if((z==0) & (t==1)){
      new = "LowZ_HighT"
    }
    if((z==0) & (t==0)){
      new = "Low_Both"
    }
    return(new)
  }
  ZKSCAN3$combo = apply(ZKSCAN3, 1, get_combo)
  ZKSCAN3$combo = as.factor(ZKSCAN3$combo)

  ZKSCAN3$OS = as.numeric(ZKSCAN3$OS)
 	ZKSCAN3$OS.time = as.numeric(ZKSCAN3$OS.time)
  
  	fit <- coxph(Surv(OS.time, OS)  ~ 1, data = ZKSCAN3)
    fit2 <- coxph(Surv(OS.time, OS)  ~ combo, data = ZKSCAN3)
    lncs = anova(fit2,fit)
    
  	row <- c("ZKSCAN_TP53", lncs[2,4], ZKSCAN3$Cancer[1])
    
    colnames(ZKSCAN3)[2] = "ZKSCAN3"
    colnames(ZKSCAN3)[3] = "TP53"

  	ZKSCAN3$OS.time = ZKSCAN3$OS.time/365
  	fit <- survfit(Surv(OS.time, OS) ~ combo, data = ZKSCAN3)
         
          s <- ggsurvplot(
          title = paste("ZKSCAN3 & TP53", ZKSCAN3$Cancer[1]),
          fit, 
          xlab = "Time (Years)", 
          #surv.median.line = "hv",
          font.main = c(16, "bold", "black"),
          font.x = c(14, "plain", "black"),
          font.y = c(14, "plain", "black"),
          font.tickslab = c(14, "plain", "black"),
          font.legend = 12,
          risk.table.fontsize = 5, 
          #legend.labs = c("Low Expression", "High Expression"),             # survfit object with calculated statistics.
          data = ZKSCAN3,      # data used to fit survival curves. 
          risk.table = TRUE,       # show risk table.
          legend = "right", 
          pval = TRUE,             # show p-value of log-rank test.
          conf.int = FALSE,        # show confidence intervals for 
                            # point estimaes of survival curves.
          xlim = c(0,5),        # present narrower X axis, but not affect
                            # survival estimates.
          break.time.by = 1,     # break X axis in time intervals by 500.
          #palette = colorRampPalette(mypal)(14), 
          palette = mypal,
          ggtheme = theme_minimal(), # customize plot and risk table with a theme.
          risk.table.y.text.col = T, # colour risk table text annotations.
          risk.table.y.text = FALSE # show bars instead of names in text annotations
                            # in legend of risk table
          )
          print(s)


          #boxplot showing difference in TP53 expression between high/low Z groups
          p <- ggboxplot(ZKSCAN3, x = "ZKSCAN3_tag", y = "TP53",
          color = "ZKSCAN3_tag",
          palette = mypal[c(4,1)], title = paste("TP53 expression", "\nExpression", ZKSCAN3$Cancer[1] , sep=" "), 
          add = "jitter", ylab = "TP53 log1p(FPKM)",  ggtheme = theme_bw())
          # Change method
          p1 = p + stat_compare_means(method = "wilcox.test")


          #boxplot showing difference in Z expresion bewteen high low TP53 groups 
          p <- ggboxplot(ZKSCAN3, x = "TP53_tag", y = "ZKSCAN3",
          color = "TP53_tag",
          palette = mypal[c(4,1)], title = paste("ZKSCAN3 expression", "\nExpression", ZKSCAN3$Cancer[1] , sep=" "), 
          add = "jitter", ylab = "ZKSCAN3 log1p(FPKM)",  ggtheme = theme_bw())
          # Change method
          p2 = p + stat_compare_means(method = "wilcox.test")

          p = p1 + p2

          print(p)

	return(row)
}

pdf("ZKSCAN3_TP53_KM_plots_TCGA_cancers_June21.pdf", width=10)
dividedWpcgs <- llply(tissues_data, getClinical, .progress = "text")
dev.off()
print("pass2")

cox_results = do.call(rbind.data.frame, dividedWpcgs)
colnames(cox_results) = c("Predictor", "LogLik_Pval", "Cancer")
cox_results$LogLik_Pval = as.numeric(cox_results$LogLik_Pval)
cox_results = as.data.table(cox_results)
cox_results = cox_results[order(LogLik_Pval)]
cox_results$cox_fdr = p.adjust(cox_results$LogLik_Pval, method="fdr")

#get in order of decreasing p-value to increasing p-value 

tissues = unique(cox_results$Cancer)

get_tissue_specific <- function(tissue){
	tis <- pcg[pcg$Cancer==tissue,]
	return(tis)
}
tissues_data <- llply(tissues, get_tissue_specific, .progress="text")

pdf("ZKSCAN3_TP53_KM_plots_TCGA_cancers_June21.pdf", width=10)
dividedWpcgs <- llply(tissues_data, getClinical, .progress = "text")
dev.off()
print("pass2")

















