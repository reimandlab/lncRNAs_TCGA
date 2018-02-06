library(data.table)
library(plyr)
library(dplyr)
library(stringr)

#1. See if any candidates have CNAs
lncswcnas = fread("fantom_lncrnas_wTCGA_CNAs_4cancers.bed")
lncswcnas = as.data.frame(lncswcnas)
lncswcnas$V12[lncswcnas$V12 == "Ovarian serous cystadenocarcinoma"] = "ovary"
lncswcnas$V12[lncswcnas$V12 == "Kidney renal clear cell carcinoma"] = "kidney"
lncswcnas$V12[lncswcnas$V12 == "Pancreatic adenocarcinoma"] = "pancreas"
lncswcnas$V12[lncswcnas$V12 == "Liver hepatocellular carcinoma"] = "liver"

#2. cands 
cands = readRDS("36_unique_cands_4cancers_TCGA_Feb6.rds")
colnames(lncswcnas)[4] = "gene"

#keep gene-cancer combinations, don't really care right now if gene has CNA
#for a different cancer where it's not a candidate 

genes = as.list(unique(cands$gene[which(cands$gene %in% lncswcnas$gene)]))
colnames(lncswcnas) = c("lnc_chr", "lnc_start", "lnc_end", "gene", "name", "Chromosome" , "Start" , 
	"End", "Num_Probes" , "Segment_Mean", "patient", "canc")

lncswcnas = dplyr::filter(lncswcnas, Num_Probes >=10, abs(Segment_Mean) >=0.2)

get_data = function(lnc){
	cancer = cands$canc[which(cands$gene == lnc)][1]
	dat = dplyr::filter(lncswcnas, canc == cancer, gene == lnc)
	#keep smallest start and largest end for each patient 
	for(i in 1:length(unique(dat$patient))){
		pat_dat = dplyr::filter(dat, patient == unique(dat$patient)[i])
	}

	return(dat)
}

lnc_cna_cancer_data = llply(genes, get_data, .progress="text")

#how to make sense of this? 
#plot segment? and where within it lies lncRNA? 

