library(data.table)
library(ggplot2)
library(ggpubr)
library(stringr)
library(dplyr)

#Data---------------------------------------------------

#[1] FMREs, 41 in total
load("dfr_encode_results_signf.rsav")
head(dfr_encode_results_signf)
fmres = dfr_encode_results_signf

#[2] coding drivers - need to filter for cancer_type==“PANCANCER” and 
#element_type==“gc19_pc.cds”, 47 in total
load("all_results_signf.rsav")
head(all_results_signf)
all_results_signf = as.data.table(all_results_signf)
all_results_signf = filter(all_results_signf, cancer_type == "PANCANCER")
all_results_signf = filter(all_results_signf, element_type == "gc19_pc.cds")
coding_drivers = all_results_signf

#[3] mutations in all CRMs, subset of these are FMREs
load("encode_merge_oct2016_mutations_sanger_neutral__PANCANCER_in_elements.rsav")
head(mutations_in_elements)
mutations_in_crms = mutations_in_elements

#[4] mutations in all CDS, subset of these are CDS drivers
load("gc19_pc.cds_oct2016_mutations__PANCANCER_in_elements.rsav")
cds_mutations = mutations_in_elements

#[5] all patients in cohort
load("patient2cancertype.rsav")
head(patient2cancer_type)

#[6] Clinical file 
clin <- fread("pcawg_specimen_histology_August2016_v6.tsv", data.table=F)
conversion <- fread("pcawgConversion.tsv", data.table=F)

#Analysis---------------------------------------------------

#1. Which patients have FMRE mutations  
z = which(mutations_in_crms$reg_id %in% fmres$id)
mutations_in_crms = mutations_in_crms[z,]
unique(mutations_in_crms$mut_patient)
z1 = which(conversion$tumor_wgs_submitter_sample_id %in% mutations_in_crms$mut_patient)
z2 = which(conversion$tumor_wgs_aliquot_id %in% mutations_in_crms$mut_patient)
z = unique(c(z1,z2))

#2. Which patients have coding mutations 
z = which(cds_mutations$reg_id %in% coding_drivers$id)
mutations_in_cds = cds_mutations[z,]


#Essentially we need a series of fisher’s exact tests

fisher.test(
	patients_with_CDS_driver_X_muts %in% all_patients, 
	patients_with_FMRE_driver_Y_muts %in% all_patients,
	alt=“greater”)

#for all pairs (CDS_driver_X, FMRE_driver_Y). 

tumor_wgs_submitter_specimen_id
 tumor_wgs_submitter_sample_id **
 tumor_wgs_aliquot_id