##################################################################################
# Active Driver Project; Script writing began on 2017.03.16
# Sublime Text 2: Use Ctrl+Shift+P and then type "Set Syntax: R" to set palette
# To mount my cluster folder in Ubuntu, double-click on "Connect to Server"
# Then type in the following: smb://storage.isilon.stg.oicr.on.ca/reimandlab
# Type in login credetials, and make sure group is set to "oicr"
# par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0) #Default values
##################################################################################

library(ggplot2)
library(reshape2)
library(survival)
library(data.table)
library(stringr)

#################################################################################
# SETUP
#################################################################################

#setwd("/.mounts/labs/reimandlab/private/users/kisaev/SUBID/gbm_ic/data/MASTER_FILES")
setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ")
dataset <- read.csv("final_ic_rnaseq.csv", header=TRUE, stringsAsFactors=FALSE, row.names=1)
dataset <- as.data.frame(t(dataset[,1:ncol(dataset)]))
name.prefix <- "(TCGA_RNASEQ)"

#setwd("/.mounts/labs/reimandlab/private/users/kisaev/SUBID/gbm_ic/data/MASTER_FILES")
#dataset <- read.csv("IC_ONLY_GSE4271_GPL96.csv", header=TRUE, stringsAsFactors=FALSE, row.names=1)
#dataset <- as.data.frame(t(dataset[,1:ncol(dataset)]))
#name.prefix <- "(GSE4271_GPL96)"

#setwd("/.mounts/labs/reimandlab/private/users/kisaev/SUBID/gbm_ic/data/MASTER_FILES")
#dataset <- read.csv("IC_ONLY_GSE4271_GPL97.csv", header=TRUE, stringsAsFactors=FALSE, row.names=1)
#dataset <- as.data.frame(t(dataset[,1:ncol(dataset)]))
#name.prefix <- "(GSE4271_GPL97)"

#setwd("/.mounts/labs/reimandlab/private/users/kisaev/SUBID/gbm_ic/data/MASTER_FILES")
#dataset <- read.csv("IC_ONLY_GSE7696.csv", header=TRUE, stringsAsFactors=FALSE, row.names=1)
#dataset <- as.data.frame(t(dataset[,1:ncol(dataset)]))
#name.prefix <- "(GSE7696)"

#setwd("/.mounts/labs/reimandlab/private/users/kisaev/SUBID/gbm_ic/data/MASTER_FILES")
#dataset <- read.csv("IC_ONLY_GSE13041_GPL96.csv", header=TRUE, stringsAsFactors=FALSE, row.names=1)
#dataset <- as.data.frame(t(dataset[,1:ncol(dataset)]))
#name.prefix <- "(GSE13041_GPL96)"

#setwd("/.mounts/labs/reimandlab/private/users/kisaev/SUBID/gbm_ic/data/MASTER_FILES")
#dataset <- read.csv("IC_ONLY_TCGA_Affymetrix_RMA.csv", header=TRUE, stringsAsFactors=FALSE, row.names=1)
#dataset <- as.data.frame(t(dataset[,1:ncol(dataset)]))
#name.prefix <- "(TCGA_AGILENT)"

#setwd("/.mounts/labs/reimandlab/private/users/kisaev/SUBID/gbm_ic/data/MASTER_FILES")
#dataset <- read.csv("IC_ONLY_TCGA_Overlap.csv", header=TRUE, stringsAsFactors=FALSE, row.names=1)
#dataset <- as.data.frame(t(dataset[,1:ncol(dataset)]))
#name.prefix <- "(TCGA_OVERLAP)"

#setwd("/.mounts/labs/reimandlab/private/users/kisaev/SUBID/gbm_ic/data/MASTER_FILES")
#dataset <- read.csv("IC_ONLY_TCGA_Non-Overlap.csv", header=TRUE, stringsAsFactors=FALSE, row.names=1)
#dataset <- as.data.frame(t(dataset[,1:ncol(dataset)]))
#name.prefix <- "(TCGA_NONOVERLAP)"

#remove duplicate patients 
dataset$pats = ""
dataset$pats = unlist(sapply(rownames(dataset), function(x){unlist(strsplit(x, "\\."))[4]}))
dataset$id = unlist(sapply(rownames(dataset), function(x){paste(unlist(strsplit(x, "\\."))[1:3], collapse="-")}))

primary = which(str_detect(dataset$pats, "01"))
#1 keep only primary tumours
dataset = dataset[primary,]

dups = dataset$id[which(duplicated(dataset$id))]
dataset = as.data.table(dataset)
dups = (dplyr::filter(dataset, id %in% dups))
dups = unique(dups$id)
#have no way to choose which one is the right one
#two smaples from same person completley different 
z = which(dataset$id %in% dups)
dataset = dataset[-z,]
dataset = as.data.frame(dataset)
rownames(dataset) = dataset$id
dataset$id = NULL
dataset$pats = NULL

source("subid_wnew_tcga_dat.R")

#Fix survival times to update them
#using the script subid_wnew_tcga_dat
for(i in 1:nrow(dataset)){
  print(dataset$os_time[i])  
  pat = dataset$patient[i]
  z = which(gbm$patient == pat)
  if(!(length(z)==0)){
  dataset$os_status[i] = gbm$OS[z]
  dataset$os_time[i] = gbm$OS.time[z]
  print(dataset$os_time[i])  
}
}

z = which(duplicated(all$patient))
pats = all$patient[z]
z = which(all$patient %in% pats)
all = all[-z,]

cancers = unique(all$type)

t = as.data.table(table(all$type))
t = filter(t, N >= 50)

all = subset(all, type %in% t$V1)
cancers = unique(all$type)

#cands -- ion channels 
#UCSC gene info
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
#z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
#ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,8]))
ucsc <- ucsc[-z,]

cands = read.csv("ION_CHANNELS_targets_and_families.csv")
colnames(ucsc)[8] = "HGNC.symbol"
cands = merge(cands, ucsc, by = "HGNC.symbol")

num_ics = which(colnames(all) %in% cands$hg19.ensGene.name2)

#################################################################################
# PART 1: SUBID
#################################################################################

get_setup = function(canc){

  #Setup: create full SubID output table
  output <- data.frame(matrix(, nrow = length(num_ics), ncol = 99))
  for (i in 1:99){
    colnames(output)[i] <- paste(i,"%", sep="")
    }
  rownames(output) <- colnames(all)[num_ics]
  output$canc = canc
  return(output)
}

setup = llply(cancers, get_setup)
names(setup) = cancers

# input: a file named "dataset" with the following structure:
# rows: each row is a patient
# columns: an os_time column, and os_status column, and a column for every gene's expression value

#Analysis: Run the CoxPH p-value optimization
subid = function(canc){

cancer = canc$canc[1]
output = canc
print(cancer)

for (j in 1:nrow(output)) {
  #for (j in 1:2) {
  print(paste(j, "out of", nrow(output)))
  canc_dat = subset(all, type == cancer)
  gene = rownames(output)[j]
  canc_surv_dat = canc_dat[,c(which(colnames(canc_dat) %in% c("OS", "OS.time", gene, "patient")))]
  canc_surv_dat$OS = as.numeric(canc_surv_dat$OS)
  canc_surv_dat$OS.time = as.numeric(canc_surv_dat$OS.time)
  print(colnames((canc_surv_dat)[4]))
  colnames(canc_surv_dat)[4] = "IC"

  #Part 1: Extract survival data for query gene, and sort based on exprssion
  survival.data <- data.frame(as.numeric(as.character(canc_surv_dat[,"OS.time"])), 
                              as.numeric(as.character(canc_surv_dat[,"OS"])), 
                              as.numeric(as.character(canc_surv_dat[,"IC"])))
  colnames(survival.data) <- c("os_time","os_status","grouping")
  survival.data <- survival.data[order(survival.data$grouping) , ]
  
  #Part 2: Run the CoxPH p-value calculator
  for (i in 1:99){
  
  print(paste(i, "of 99"))
  curent.coxph <- NULL
  
  x=round((((i)/100)*nrow(survival.data)), digits = 0)
  survival.data$grouping[1:x] <- "Low"            
  y=x+1
  survival.data$grouping[y:nrow(survival.data)] <- "High"
    
  current.coxph <- coxph(Surv(os_time, os_status) ~ grouping, data=survival.data)
  coxph.pvalue <- summary(current.coxph)$coefficients[5]
  output[j,i] <- coxph.pvalue
  }  
  
  }#end analysis of each gene individually 
  return(output)
}

#subid_res = llply(setup, subid)

#setwd("/.mounts/labs/reimandlab/private/users/idzneladze/gbm_ic/data/MASTER_OUTPUT")
#write.csv(output, file=paste("karin_SUBID", "_full_subid.csv", sep=""), quote=FALSE)
#saveRDS(subid_res, file="all_cancers_subid_results_march19.rds")

#################################################################################

subid_res = readRDS("all_cancers_subid_results_march19.rds")

#UCSC gene info
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
#z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
#ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,8]))
ucsc <- ucsc[-z,]

get_name_pcg = function(pcg){
  z = which(ucsc$hg19.ensGene.name2 == pcg)
  if(length(z)>1){
    z = z[1]
  }
  return(ucsc$hg19.ensemblToGeneName.value[z])
}

get_ensg_pcg = function(pcg){
  z = which(ucsc$hg19.ensemblToGeneName.value == pcg)
  if(length(z)>1){
    z = z[1]
  }
  return(ucsc$hg19.ensGene.name2[z])
}


get_output2 = function(canc){

output2 <- data.frame(matrix(, nrow = length(num_ics), ncol = 3))
colnames(output2) <- c("median_p", "subid_p", "subid_%")
rownames(output2) <- colnames(all)[num_ics]
print(head(output2))
print(canc$canc[1])

for (i in 1:(nrow(output2))){
  output2[i,"median_p"] <- canc[i,50]
  output2[i,"subid_p"] <- min(canc[i,20:80])
  output2[i,"subid_%"] <- min(which(canc[i,20:80]==min(canc[i,20:80]))+19)
}

output2$gene = rownames(output2)
output2$gene_name = unlist(sapply(output2$gene, get_name_pcg))
output2 = as.data.table(output2)
output2 = output2[order(median_p, subid_p)]
output2 <- output2[order(output2$subid_p), ]
output2$canc = canc$canc[1]
output2$subid_fdr = p.adjust(output2$subid_p, method="fdr")
output2$median_fdr = p.adjust(output2$median_p, method="fdr")
output2 = output2[order(subid_p)]
output2 = output2[order(median_p)]
output2 = output2[order(subid_p)]

return(output2)

}

res_list_subids = llply(subid_res, get_output2, .progress="text")
res_list_subids2 = as.data.table(ldply(res_list_subids))
res_list_subids2[,1] = NULL

#setwd("/.mounts/labs/reimandlab/private/users/idzneladze/gbm_ic/data/MASTER_OUTPUT")
#write.csv(output2, file=paste(name.prefix, "_subid.csv", sep=""), quote=FALSE)

saveRDS(res_list_subids2, file="subid_results_KI_IC_march20.rds")


