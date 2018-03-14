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

files_wid = llply(files, add_id, .progress="text")#469 unique samples in concatenated maf file 
#convert list of dataframes into one dataframe 
df <- ldply(files_wid, data.frame)

##########################################################
#2. Get overview of mutations with maftools
##########################################################
ov = read.maf(maf = df)
write.mafSummary(maf = ov, basename = 'ov')

#Plotting MAF summary. 
pdf("OV_maf_summary.pdf")
plotmafSummary(maf = ov, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

#We will draw oncoplots for top ten mutated genes.
pdf("Oncoplot_OV.pdf", width=12, height=11)
oncoplot(maf = ov, top = 10, fontSize = 12)
dev.off()




