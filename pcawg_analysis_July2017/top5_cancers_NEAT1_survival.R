#top5_cancers_NEAT1_survival.R

#Karina Isaev
#August 11th, 2017

#Purpose: using the top 5 cancers selected for analysis, 
#run survival analysis in a pancancer approach with cancer 
#type as covariate as Neat1 is highly expressed in all cancers

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#RNA files used here obtained from script: 
#pcawg_analysis_July2017/top5_cancers_extraction_script3.R 
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Preamble#-------------------------------------------------
options(stringsAsFactors=F)

#Later on - incorporate FANTOM5 and CRISPRi lncRNAs 

#Libraries#------------------------------------------------
library("colorout")
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
library(ggthemes)

mypal = pal_npg("nrc", alpha = 0.7)(10)

#Data-------------------------------------------------------

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

#Clinical file 
#clin <- fread("pcawg_specimen_histology_August2016_v6.tsv", data.table=F)
#conversion <- fread("pcawgConversion.tsv", data.table=F)

#RNA-Seq file 
#lncRNA
lnc_rna <- readRDS("5607_pcawg_lncRNAs_RNASeq_data.rds")
lnc_rna <- as.data.frame(lnc_rna)
lnc_rna$patient <- rownames(lnc_rna)

#PCGs
pcg_rna <- readRDS("20166_pcawg_PCGs_RNASeq_data.rds")
pcg_rna <- as.data.frame(pcg_rna)
pcg_rna$patient <- rownames(pcg_rna)

#remove duplicated column names 
dups <- colnames(pcg_rna)[which(duplicated(colnames(pcg_rna)))]   
#save them in a list for future reference 
pcg_rna <- pcg_rna[,-(which(colnames(pcg_rna) %in% dups))]

#Clinical file - available only for 485/497 patients 
clin <- readRDS("Jan26_PCAWG_clinical")
z <- which(clin$icgc_donor_id %in% rownames(lnc_rna))
clin <- clin[z,]

lnc_rna <- lnc_rna[which(rownames(lnc_rna) %in% clin$icgc_donor_id),] #485 patients remain
pcg_rna <- pcg_rna[which(rownames(pcg_rna) %in% clin$icgc_donor_id),] #485 patients remain 

#---------------------------------------------------------
#Processing 
#---------------------------------------------------------

#1. For each patient add survival status and days since last seen 
lnc_rna$status <- ""
lnc_rna$time <- ""
lnc_rna$sex <- ""

#lncs
for(i in 1:nrow(lnc_rna)){
	pat <- rownames(lnc_rna)[i]
	z <- which(clin$icgc_donor_id %in% pat)
	lnc_rna$status[i] <- clin$donor_vital_status[z]
	lnc_rna$sex[i] <- clin$donor_sex[z]
	t <- clin$donor_survival_time[z]
	if(is.na(t)){
        t <- clin$donor_interval_of_last_followup[z]
        }
        lnc_rna$time[i] <- t
}

pcg_rna$status <- ""
pcg_rna$time <- ""
pcg_rna$sex <- ""

#pcgs
for(i in 1:nrow(pcg_rna)){
	pat <- rownames(pcg_rna)[i]
	z <- which(clin$icgc_donor_id %in% pat)
	pcg_rna$status[i] <- clin$donor_vital_status[z]
	lnc_rna$sex[i] <- clin$donor_sex[z]
	t <- clin$donor_survival_time[z]
	if(is.na(t)){
        t <- clin$donor_interval_of_last_followup[z]
        }
        pcg_rna$time[i] <- t
}

#Looking only at NEAT1 
neat1 <- lnc_rna[, c((which(colnames(lnc_rna)=="NEAT1")), 5608:5612)]
neat1[,1] <- log1p(neat1[,1])

#3. Add Median cutoff tag High or Low to each patient per each gene 
neat1$median <- ""
median2 <- median(neat1[,1])
#median <- quantile(as.numeric(gene_data$GeneE), 0.4)
for(y in 1:nrow(neat1)){
	genexp <- neat1[y,1]
	if(genexp >= median2){
		neat1$median[y] <- 1
		}
	if(genexp < median2){
		neat1$median[y] <- 0
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

results_cox <- as.data.frame(matrix(ncol=4)) ; colnames(results_cox) <- c("gene", "coef", "HR", "pval")

#make list of plots, add title to each plot indicating which gene and cancer is studied 
#change legend and colours 

pdf("survival_results_usingMedian_NEAT1_Pancan_lncRNAs.pdf", pointsize=6, width=15, height=14)
require(gridExtra)

#Just looking at Neat1
	gene <- colnames(neat1)[1]
    #plot boxplot showing difference between the two groups and sex
	title <- "Neat1 Expression"
	g <- ggboxplot(neat1, x= "median", y="NEAT1", palette=mypal, order=c("0", "1"), fill = "canc")
	g <- ggpar(g, font.legend = c(8, "plain", "black")) 
	g <- g + labs(title = title, y="log1p(FPKM)", x="Median") + 
     	theme(plot.title = element_text(hjust = 0.5))
     	print(g)

     	#cox
	  	neat1$status[neat1$status=="alive"] <- 0
      	neat1$status[neat1$status=="deceased"] <- 1
      	neat1$status <- as.numeric(neat1$status)
      	neat1$time <- as.numeric(neat1$time)
      
      	#cox regression 
      	res.cox <- coxph(Surv(time, status) ~ median + canc, data = neat1)
      	row <- c(gene, summary(res.cox)$coefficients[1,c(1,2,5)])
      	names(row) <- names(results_cox)
      	results_cox <- rbind(results_cox, row)
        
         #plot survival plot
          fit <- survfit(Surv(time, status) ~ median + canc, data = neat1)
          s <- ggsurvplot(
          fit,                     # survfit object with calculated statistics.
          data = neat1,      # data used to fit survival curves. 
          risk.table = TRUE,       # show risk table.
          legend = "right", 
          pval = TRUE,             # show p-value of log-rank test.
          conf.int = FALSE,        # show confidence intervals for 
                            # point estimaes of survival curves.
          xlim = c(0,2000),        # present narrower X axis, but not affect
                            # survival estimates.
          break.time.by = 500,     # break X axis in time intervals by 500.
          palette = colorRampPalette(mypal)(14), 
          ggtheme = theme_minimal(), # customize plot and risk table with a theme.
          risk.table.y.text.col = T, # colour risk table text annotations.
          risk.table.y.text = FALSE # show bars instead of names in text annotations
                            # in legend of risk table
          )
          

		  print(s)

dev.off()





















