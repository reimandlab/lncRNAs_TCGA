#top5_cancers_median5fpkm_specificFind.R

#Karina Isaev
#August 11th, 2017

#Purpose: using the top 5 cancers selected for analysis, 
#run survival analysis in a pancancer approach with cancer 
#type as covariate as Neat1 is highly expressed in all cancers

#last updated August 28th

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
library(plyr)

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
colnames(fantom)[1] = "gene"

#Clinical file 
#clin <- fread("pcawg_specimen_histology_August2016_v6.tsv", data.table=F)
#conversion <- fread("pcawgConversion.tsv", data.table=F)

#RNA-Seq file 
#lncRNA
lnc_rna <- readRDS("6028_pcawg_lncRNAs_RNASeq_data.rds")
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
#Find cancer sepcific lncRNAs with median E >= 5 FPKM 
#---------------------------------------------------------

#Write function that takes a dataframe, calculates medians and 
#output list of genes with median greater than that 

#lnc_rna[,1:6013] = log1p(lnc_rna[,1:6013])

#---------------------------------------------------------
#Subset lncRNA Expression dataset to those lncRNAs with 
#high expression in at leat one canc 215 total lncRNAs
#---------------------------------------------------------
genes = readRDS("final_candidates_10batches_of1000CV_KIRC_cancer_patientsJan29.RDS")
cands = as.data.frame(table(genes$gene))
colnames(cands)[1] = "gene"
cands$name = ""
for(i in 1:nrow(cands)){
  n = fantom$CAT_geneName[which(fantom$gene %in% cands$gene[i])]
  cands$name[i] = n
}

#only keep those that appreared in 9-10 of 10 batches 
cands = subset(cands, Freq >=9)

#top3genes = c("ENSG00000227486", "ENSG00000227544", "ENSG00000232124", "ENSG00000235572", "ENSG00000249662", 
 # "ENSG00000258082", "ENSG00000265369")

lnc_rna <- lnc_rna[,c((which(colnames(lnc_rna) %in% cands$gene)), 6014,6015)] 

#For each patient add survival status and days since last seen 
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

lnc_rna = subset(lnc_rna, canc =='Kidney Adenocarcinoma, clear cell type')

#---------------------------------------------------------
#Run survival analysis on each gene-cancer combo
#from high_lncs 
#---------------------------------------------------------
results_cox <- as.data.frame(matrix(ncol=5)) ; colnames(results_cox) <- c("gene", "coef", "HR", "pval", "canc")

lnc_rna$status[lnc_rna$status=="alive"] <- 0
lnc_rna$status[lnc_rna$status=="deceased"] <- 1
lnc_rna$status <- as.numeric(lnc_rna$status)
lnc_rna$time <- as.numeric(lnc_rna$time)

for(i in 1:14){
  #1. Subset lnc_rna to those patients in cancer
  df <- lnc_rna[,c(i, 15:19)]
  
  #2. Add Median cutoff tag High or Low to each patient per each gene 
  df$median <- ""
  #median2 <- quantile(as.numeric(df[,1]), 0.75)
  median2 <- median(df[,1])
  if(median2 == 0){
    median2 = mean(df[,1])
  }
  for(y in 1:nrow(df)){
    genexp <- df[y,1]
    if(genexp >= median2){
      df$median[y] <- 1
      }
    if(genexp < median2){
      df$median[y] <- 0
      }
    } 

  gene <- colnames(df)[1]
  #cox
        #cox regression 
        colnames(df)[1] = "geneexpression"
        res.cox <- coxph(Surv(time, status) ~ median, data = df)
        #first check that model meets proportionality assumption
        #res.cox <- coxph(Surv(time, status) ~ df[,1], data = df)
        #testph <- cox.zph(res.cox)
        #p = testph$table[1,3]
        #if(p >= 0.05){
        row <- c(gene, summary(res.cox)$coefficients[1,c(1,2,5)], df$canc[1])
        names(row) <- names(results_cox)
        results_cox <- rbind(results_cox, row)  
#}
}

results_cox <- results_cox[-1,]

##Full - order plots by decreasing pvalue 
##+++++++++++++++++++++++++++++

pdf("Validating14_lncRNAs_fromTCGA_KIRC_Feb22018.pdf", pointsize=6, width=9, height=8)
require(gridExtra)

for(i in 1:14){
  #1. Subset lnc_rna to those patients in cancer
  df <- lnc_rna[,c(i, 15:19)]

  #2. Add Median cutoff tag High or Low to each patient per each gene 
  df$median <- ""
  #median2 <- quantile(as.numeric(df[,1]), 0.75)
  median2 <- median(df[,1])
  if(median2 == 0){
    median2 = mean(df[,1])
  }
  for(y in 1:nrow(df)){
    genexp <- df[y,1]
    if(genexp >= median2){
      df$median[y] <- 1
      }
    if(genexp < median2){
      df$median[y] <- 0
      }
    } 
    
  gene <- colnames(df)[1]
  gene = fantom$CAT_geneName[which(fantom$gene == gene)]
  #plot boxplot showing difference between the two groups and sex
  title <- paste(gene, df$canc[1], "Expression")
  colnames(df)[1] <- "Gene"
  g <- ggboxplot(df, x= "median", y="Gene", palette=mypal[c(4,1)], order=c("0", "1"), fill = "median",  add = "jitter")
  g <- g + stat_compare_means()
  g <- ggpar(g, font.legend = c(10, "plain", "black")) 
  g <- g + labs(title = title, y="FPKM", x="Median") + 
      theme(plot.title = element_text(hjust = 0.5))
      print(g)

  #cox
       
          #plot survival plot
          fit <- survfit(Surv(time, status) ~ median, data = df)
          s <- ggsurvplot(
          fit, 
          surv.median.line = "hv",
          font.main = c(16, "bold", "black"),
          font.x = c(14, "plain", "black"),
          font.y = c(14, "plain", "black"),
          font.tickslab = c(14, "plain", "black"),
          font.legend = 12,
          risk.table.fontsize = 5, 
          main = paste(gene, df$canc[1]),       
          legend.labs = c("Low Expression", "High Expression"),             # survfit object with calculated statistics.
          data = df,      # data used to fit survival curves. 
          risk.table = TRUE,       # show risk table.
          legend = "right", 
          pval = TRUE,             # show p-value of log-rank test.
          conf.int = FALSE,        # show confidence intervals for 
                            # point estimaes of survival curves.
          xlim = c(0,2000),        # present narrower X axis, but not affect
                            # survival estimates.
          break.time.by = 500,     # break X axis in time intervals by 500.
          #palette = colorRampPalette(mypal)(14), 
          palette = mypal[c(4,1)],
          ggtheme = theme_minimal(), # customize plot and risk table with a theme.
          risk.table.y.text.col = T, # colour risk table text annotations.
          risk.table.y.text = FALSE # show bars instead of names in text annotations
                            # in legend of risk table
          )
          print(s)  
}


dev.off()






