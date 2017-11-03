
#top5_cancers_extraction_script0.R

#Karina Isaev
#July 19th, 2017

#Purpose: Using PCAWG, extract the top 5 cancer types with the 
#most patients and good distribution of samples across
#multiple histological subtypes 

#For each cancer type, identify list of 100-200 candidate 
#lncRNAs that wiill be used for further survival and co-expression
#analysis 

#Script0 - extracts cancer types and plot distributions 

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

mypal = pal_npg("nrc", alpha = 0.7)(10)

aid = c(
"1494bb16-f1f0-42a4-b10e-c383574cbc8b", 
"461df2ae-fcf1-4b93-be0a-c14954fe7c42", 
"4e7cdeda-6dc1-4f17-b853-72a68e5aa7e1", 
"890e840c-1d1d-4874-a8eb-f9d9a2b50a1c", 
"915cbb43-9e00-433d-818f-531011bea57e", 
"995a1ad2-faca-4a37-a59d-e62455985afb"
)

conversion <- fread("pcawgConversion.tsv", data.table=F)

#so aliqot has all the same IDs from submitter and so can just use aliqot ids 

###"TUMOUR SAMPLES"
z <- which(conversion$broad_tar_variant_calling_file_name_prefix %in% aid)
tum_pats <- conversion$icgc_donor_id[z]

