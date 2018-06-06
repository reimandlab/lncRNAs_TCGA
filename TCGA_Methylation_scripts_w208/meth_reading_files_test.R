
options(stringsAsFactors=F)
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

start_time <- Sys.time()

f = readRDS("ACC_all_methylation_.rds")

end_time <- Sys.time()

end_time - start_time

print(end_time - start_time)
print("done")

#1 - remove probes that are NA for all patients 

f$probe_there = ""
z = which(is.na(f$beta))
f$probe_there[z] = "not_there"
f$probe_there[-z] = "there"

part_1 = as.data.table(table(f$probe, f$probe_there))

#num of patients
pats = length(unique(f$patient))

#which probes NA in all 80 patients?
probes_rm = as.data.table(filter(part_1, V2 == "not_there",  N == pats))$V1
z = which(f$probe %in% probes_rm)
f = f[-z,]


#2 - remove probes that are all 0s for all patients 
f$probe_there = ""
z = which(as.numeric(f$beta) < 0.5) #check if all patients have unmethylated probe 
f$probe_there[z] = "not_there"
f$probe_there[-z] = "there"

part_2 = as.data.table(table(f$probe, f$probe_there))

probes_rm = as.data.table(filter(part_2, V2 == "not_there",  N == pats))$V1
z = which(f$probe %in% probes_rm)
f = f[-z,]

library(reshape2)

#pcgs enriched in high risk lncRNAs 
start_time <- Sys.time()

methylation_data = acast(f, patient ~ probe, function(x) {sort(as.character(x))[1]},
      value.var = 'beta', fill = 'na')

end_time <- Sys.time()

print(end_time - start_time) #time to make a matrix from beta values 

#3 - change patient IDs to shorter version 

	pats = as.data.frame(rownames(methylation_data))
	colnames(pats)[1] = "patient"
	get_source = function(id){
		source = unlist(strsplit(as.character(id), '-'))[4]
		source = substr(source, 1, 2)
		return(source)
		}
	pats$source = llply(pats$patient, get_source, .progress="text")
	#only keep primary solid tumour samples 
	pats = dplyr::filter(pats, source =="01")
	clean_tcga_id = function(id){
		s1 = unlist(strsplit(id, '-'))[1]
		s2 = unlist(strsplit(id, '-'))[2]
		s3 = unlist(strsplit(id, '-'))[3]
		return(paste(s1, s2, s3, sep="-"))
		}	
	pats$new_patient = llply(pats$patient, clean_tcga_id, .progress="text")

change_rows = function(pat){
	z = which(pats$patient == pat)
	new = pats$new_patient[z]
	return(new)
}

rownames(methylation_data) = unlist(llply(rownames(methylation_data), change_rows, .progress="text"))
methylation_data = as.data.frame(methylation_data)
methylation_data$canc = f$cancer[1]

canc = f$cancer[1]
saveRDS(methylation_data, file = paste(canc, "methylation_matrix.rds", sep="_"))



