library(data.table)
library(ggplot2)
library(ggpubr)
library(stringr)
library(dplyr)
library(gridExtra)
library(grid)
library(Hmisc)
library(GenomicRanges)

#---TCGA directory 


#Data---------------------------------------------------

#cds mutations new june 12 
coding_drivers = fread("july2_cds_drivers.txt")

#[3] mutations in all CRMs, subset of these are FMREs

#old file:

#load("encode_merge_oct2016_mutations__PANCANCER_in_elements.rsav")
#head(mutations_in_elements)
#mutations_in_crms = mutations_in_elements

#new file June 12: 

load("july12_encode_merge__patient_element_snv_list.rsav")
mutations_in_crms = (patient_element_snv_list)

#crm mutations new june 12 
#fmre mutations 
fmres = fread("july12_fmre_drivers.txt")

#[4] mutations in all CDS, subset of these are CDS drivers
load("july12_gc19_pc.cds__patient_element_snv_list.rsav")
cds_mutations = patient_element_snv_list

#[5] all patients in cohort
load("patient2cancertype.rsav")
head(patient2cancer_type) #1844 all together 

patient_table = fread("patient_table.txt")

#[] cds_lnc_drivers

#[] july12_gc19_pc.cds__coords.rsav
load("july12_gc19_pc.cds__coords.rsav")
coords = element_coords
coords = as.data.table(coords)

#[] july12_lncrna.ncrna__coords.rsav
load("july12_lncrna.ncrna__coords.rsav")
lnc_coords = element_coords
lnc_coords = as.data.table(lnc_coords)

#[] july12_encode_merge__coords.rsav
load("july12_encode_merge__coords.rsav")
element_coords = element_coords
element_coords = as.data.table(element_coords)

#[] july12_lncrna.ncrna__patient_element_snv_list.rsav
load("july12_lncrna.ncrna__patient_element_snv_list.rsav")
lncrna_mutations = patient_element_snv_list

#[6] Clinical file 
clin <- fread("pcawg_specimen_histology_August2016_v6.tsv", data.table=F)
conversion <- fread("pcawgConversion.tsv", data.table=F)

#Analysis---------------------------------------------------

#1. Which patients have FMRE mutations  
z = which(names(mutations_in_crms) %in% fmres$id)
mutations_in_crms = mutations_in_crms[z]

#2. Which patients have coding mutations 
z = which(names(cds_mutations) %in% coding_drivers$id)
mutations_in_cds = cds_mutations[z]

#RESULTS_from_001_------------------------------------------------

results_pairs = fread("819_fmre_cds_pairs_fishers_analysis_with_lncRNAs_July12nd_KI.txt")

library(plyr)
library(dplyr)

#prognostic enahcners from TCGA paper 
prognostic = fread("enhancers_prognostic_TCGA_paper.txt")

prognostic$chr = unlist(lapply(prognostic$Prognostic_enhancer, function(x){unlist(strsplit(x, ":"))[1]}))
prognostic$coords = unlist(lapply(prognostic$Prognostic_enhancer, function(x){unlist(strsplit(x, ":"))[2]}))
prognostic$start = unlist(lapply(prognostic$coords, function(x){unlist(strsplit(x, "-"))[1]}))
prognostic$end = unlist(lapply(prognostic$coords, function(x){unlist(strsplit(x, "-"))[2]}))

#keep reference to go back to to check which cancer type it was
#prognostic in 

prognostic_reference = prognostic

prognostic = prognostic[,c("chr", "start", "end", "Disease")]
prognostic = as.data.frame(prognostic)

#convert to granges 
prognostic = makeGRangesFromDataFrame(prognostic)

#fmre coordiantes 
fmre_coord = fmres[,1:2]
fmre_coord$chr = unlist(lapply(fmre_coord$id, function(x){unlist(strsplit(x, "::"))[2]}))
fmre_coord$chr = unlist(lapply(fmre_coord$chr, function(x){unlist(strsplit(x, ":"))[1]}))

fmre_coord$coords = unlist(lapply(fmre_coord$id, function(x){unlist(strsplit(x, "::"))[2]}))
fmre_coord$coords = unlist(lapply(fmre_coord$coords, function(x){unlist(strsplit(x, ":"))[2]}))

fmre_coord$start = unlist(lapply(fmre_coord$coords, function(x){unlist(strsplit(x, "-"))[1]}))
fmre_coord$end = unlist(lapply(fmre_coord$coords, function(x){unlist(strsplit(x, "-"))[2]}))

fmre_coord = fmre_coord[,c("chr", "start", "end")]

#zsckan
z = which((fmre_coord$chr == "chr6") & (fmre_coord$start == 27870028))
zsckan =fmre_coord[z,]
zsckan = as.data.frame(zsckan)
zsckan = makeGRangesFromDataFrame(zsckan)
intersect(zsckan, prognostic) #none 

#what about all fmres? 
fmre_cords = makeGRangesFromDataFrame(fmre_coord)
intersect(fmre_cords, prognostic) #none 


#what about just how many FMREs overlap enhancers?  
enhancers = fread("enhancers_ALL_TCGA_paper.txt")
colnames(enhancers)[1] = "Prognostic_enhancer"
enhancers$chr = unlist(lapply(enhancers$Prognostic_enhancer, function(x){unlist(strsplit(x, ":"))[1]}))
enhancers$coords = unlist(lapply(enhancers$Prognostic_enhancer, function(x){unlist(strsplit(x, ":"))[2]}))
enhancers$start = unlist(lapply(enhancers$coords, function(x){unlist(strsplit(x, "-"))[1]}))
enhancers$end = unlist(lapply(enhancers$coords, function(x){unlist(strsplit(x, "-"))[2]}))
enhancers_reference = enhancers

enhancers = enhancers[,c("chr", "start", "end")]
enhancers = as.data.frame(enhancers)

#convert to granges 
enhancers = makeGRangesFromDataFrame(enhancers) #none? 

countOverlaps(fmre_cords, enhancers)
findOverlaps(fmre_cords, enhancers)
subsetByOverlaps(fmre_cords, enhancers)
















