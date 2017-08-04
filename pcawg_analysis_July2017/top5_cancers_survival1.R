#top5_cancers_survival1.R

#Karina Isaev
#August 1st, 2017

#Purpose: Using PCAWG, extract the top 5 cancer types with the 
#most patients and good distribution of samples across
#multiple histological subtypes 

#For each cancer type, identify list of 100-200 candidate 
#lncRNAs that wiill be used for further survival and co-expression
#analysis 

#Survival analysis using the 9 lncRNAs highly expressed in all cancer types
#and 50 lncRNAs specific to one cancer type 

#Preamble#-------------------------------------------------
options(stringsAsFactors=F)

#Later on - incorporate FANTOM5 and CRISPRi lncRNAs 

#Libraries#------------------------------------------------
library(colorout)
library(data.table)
library(survival)
library(ggplot2)
library(ggfortify)
library(cluster)
library(qqman)
library(dplyr)
library(rafalib)
library(RColorBrewer) 
library(genefilter)
library(gplots) ##Available from CRAN
library(survminer)
library(MASS)
library(Hmisc)
library(gProfileR)
library(wesanderson)
library(ggsci)
library(gridExtra)
library(ggpubr)
library(factoextra)

mypal = pal_npg("nrc", alpha = 0.7)(10)

#---------------------------------------------------------
#Data
#---------------------------------------------------------

#lncRNA Data 
lncs <- fread("logged_geneExpression_50unique_specific_lncRNAs.txt", data.table=F, sep=";")
pats <- unique(lncs$Patient) #497 unique patients 

lncs_all <- fread("logged_geneExpression_9unique_common_lncRNAs.txt", data.table=F, sep=";")
pats <- unique(lncs_all$Patient) #497 unique patients 

#Clinical file 
clin <- readRDS("Jan26_PCAWG_clinical")
z <- which(clin$icgc_donor_id %in% pats)
clin <- clin[z,]

#Subset lncRNA data to include patients with survival data
lncs <- subset(lncs, lncs$Patient %in% clin$icgc_donor_id) #485 patients left 
lncs_all <- subset(lncs_all, lncs_all$Patient %in% clin$icgc_donor_id) #485 patients left 

#---------------------------------------------------------
#Processing 
#---------------------------------------------------------

#1. For lncs file, subset such that the Patient column 
#matched the ref column 
lncs <- subset(lncs, lncs$Cancer == lncs$ref)

#2. For each patient add survival status and days since last seen 
lncs$status <- ""
lncs$time <- ""

for(i in 1:nrow(lncs)){
	pat <- lncs$Patient[i]
	z <- which(clin$icgc_donor_id %in% pat)
	lncs$status[i] <- clin$donor_vital_status[z]
	t <- clin$donor_survival_time[z]
	if(is.na(t)){
        t <- clin$donor_interval_of_last_followup[z]
        }
        lncs$time[i] <- t
}

#3. Add Median cutoff tag High or Low to each patient per each gene 
lncs$median <- ""
for(i in 1:50){
	gene <- unique(lncs$Gene)[i]
	gene_data <- subset(lncs, lncs$Gene %in% gene)
	median <- median(as.numeric(gene_data$GeneE))
	#median <- quantile(as.numeric(gene_data$GeneE), 0.4)
	for(y in 1:nrow(gene_data)){
		genexp <- gene_data$GeneE[y]
		if(genexp >= median){
			gene_data$median[y] <- 1
		}
		if(genexp < median){
			gene_data$median[y] <- 0
		}
	#add new info into the original lnc dataset
	z <- which(lncs$Patient == gene_data$Patient[y] & lncs$Gene == gene_data$Gene[y])
	lncs$median[z] <- gene_data$median[y]
	}	
}

#---------------------------------------------------------
#Analysis - univariate - cancer specific lncRNAs
#---------------------------------------------------------

#For each unique gene-reference cancer combintation 
#dichotomize patients and make plot showing E difference 
#between high and low groups 
#conduct survival analysis using univariate cox 
#make KM plot 

results_cox <- as.data.frame(matrix(ncol=5)) ; colnames(results_cox) <- c("gene", "coef", "HR", "pval","cancer")

#make list of plots, add title to each plot indicating which gene and cancer is studied 
#change legend and colours 

pdf("survival_results_usingMedian_50unique_lncRNAs.pdf", pointsize=6, width=10, height=10)
require(gridExtra)

for(i in 1:50){
	gene <- unique(lncs$Gene)[i]
	gene_data <- subset(lncs, lncs$Gene %in% gene)

    #plot boxplot showing difference between the two groups 
	title <- paste(unique(gene_data$Gene), unique(gene_data$Cancer)) 
	g <- ggboxplot(gene_data, x= "median", y="GeneE", add = "jitter", fill="median", palette=mypal, order=c("0", "1"))
	g <- ggpar(g, legend = "none") 
	g <- g + labs(title = title, y="log1p(FPKM)", x="Median") + 
     	theme(plot.title = element_text(hjust = 0.5))
     	print(g)

     	#cox
	  	gene_data$status[gene_data$status=="alive"] <- 0
      	gene_data$status[gene_data$status=="deceased"] <- 1
      	gene_data$status <- as.numeric(gene_data$status)
      	gene_data$time <- as.numeric(gene_data$time)
      
      	#cox regression 
      	res.cox <- coxph(Surv(time, status) ~median, data = gene_data)
      	row <- c(gene, summary(res.cox)$coefficients[c(1,2,5)], unique(gene_data$ref))
      	names(row) <- names(results_cox)
      	results_cox <- rbind(results_cox, row)
        
         #plot survival plot
          fit <- survfit(Surv(time, status) ~ median, data = gene_data)
          
          s <- ggsurvplot(
          fit,                     # survfit object with calculated statistics.
          data = gene_data,      # data used to fit survival curves. 
          risk.table = TRUE,       # show risk table.
          legend = "right", 
          pval = TRUE,             # show p-value of log-rank test.
          conf.int = FALSE,        # show confidence intervals for 
                            # point estimaes of survival curves.
          xlim = c(0,2000),        # present narrower X axis, but not affect
                            # survival estimates.
          break.time.by = 500,     # break X axis in time intervals by 500.
          palette = mypal, 
          ggtheme = theme_minimal(), # customize plot and risk table with a theme.
          risk.table.y.text.col = T, # colour risk table text annotations.
          risk.table.y.text = FALSE # show bars instead of names in text annotations
                            # in legend of risk table
          )
          

		  print(s)

}

dev.off()

results_cox <- results_cox[-1,]
results_cox$pval <- as.numeric(results_cox$pval)
results_cox <- as.data.table(results_cox)
results_cox <- results_cox[order(pval)]

results_cox[which(results_cox$pval <=0.1)]


#---------------------------------------------------------
#Analysis - univariate - cancer common lncRNAs
#---------------------------------------------------------



