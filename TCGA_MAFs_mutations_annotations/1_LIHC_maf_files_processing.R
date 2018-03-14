###source_file.R

#Preamble#-------------------------------------------------
options(stringsAsFactors=F)

#Later on - incorporate FANTOM5 and CRISPRi lncRNAs 

#Libraries#------------------------------------------------
library(data.table)
library(plyr)
library(survival)
library(ggplot2)
library(ggfortify)
library(cluster)
library(dplyr)
library(rafalib)
library(RColorBrewer) 
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
library(tidyr)
library(cowplot)
library(broom)
library(tidyverse)
library(limma)
library(stringr)
library(scater)
require(maftools)

mypal = pal_npg("nrc", alpha = 0.7)(10)

##########################################################
#1. Concatenate all the files into one maf file 
##########################################################

files = as.list(list.files(pattern = "maf.txt"))

add_id = function(f){
ff = fread(f, data.table=F)
id = str_sub(f, start=1, end=-9)
ff$patient_id_main = id
return(ff)
}

files_wid = llply(files, add_id, .progress="text")#373 unique patients in concatenated maf file 
#convert list of dataframes into one dataframe 
df <- ldply(files_wid, data.frame)
unique(df$patient_id_main)

##########################################################
#2. Get overview of mutations with maftools
##########################################################

lihc = read.maf(maf = df)
write.mafSummary(maf = lihc, basename = 'lihc')

#Plotting MAF summary. 
pdf("LIHC_maf_summary.pdf")
plotmafSummary(maf = lihc, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

#We will draw oncoplots for top ten mutated genes.
pdf("Oncoplot_LIHC.pdf", width=12, height=11)
oncoplot(maf = lihc, top = 20, fontSize = 12)
dev.off()

#Lets plot lollipop plot for KRAS, which is one of the most frequent mutated gene in PAAD
pdf("TTN_lollopop_LIHC.pdf", width=15)
TTN.lpop = lollipopPlot(maf = lihc, gene = 'TTN', AACol = 'Protein_Change', showMutationRate = TRUE, domainLabelSize = 3, defaultYaxis = FALSE)
dev.off()

pdf("TP53_lollopop_LIHC.pdf", width=15)
TP53.lpop = lollipopPlot(maf = lihc, gene = 'TP53', AACol = 'Protein_Change', showMutationRate = TRUE, domainLabelSize = 3, defaultYaxis = FALSE)
dev.off()

pdf("CTNNB1_lollopop_LIHC.pdf", width=15)
CTNNB1.lpop = lollipopPlot(maf = lihc, gene = 'CTNNB1', AACol = 'Protein_Change', showMutationRate = TRUE, domainLabelSize = 3, defaultYaxis = FALSE)
dev.off()

lihc.sig = oncodrive(maf = lihc, AACol = 'Protein_Change', minMut = 5, pvalMethod = 'zscore')
plotOncodrive(res = lihc.sig, fdrCutOff = 0.1, useFraction = TRUE)





















