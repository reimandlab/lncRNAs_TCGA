#check if FPKM similarity corresponds to rank similarity 

#source code
source("check_lnc_exp_cancers.R")
library(corrplot)

#COSMIC cancer gene census
census = read.csv("Census_allFri_Jul_13_16_55_59_2018.csv")
#get ensg
get_census_ensg = function(genes){
  glist = unlist(strsplit(genes, ","))
  z = which(str_detect(glist, "ENSG"))
  ensg = glist[z]
  return(ensg)
}
census$ensg = sapply(census$Synonyms, get_census_ensg)

#HiC data 
load("hic_data.rsav")
#remove rownames
rownames(hic_data) = c(1:nrow(hic_data))

#------FEATURES-----------------------------------------------------

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = subset(allCands, data == "TCGA") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 
allCands$combo = unique(paste(allCands$gene, allCands$cancer, sep="_"))

val_cands = read.csv("175_lncRNA_cancers_combos_23_cancer_types_july5.csv")
val_cands = as.data.table(val_cands)
val_cands = subset(val_cands, data == "PCAWG") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 
val_cands$combo = unique(paste(val_cands$gene, val_cands$cancer, sep="_"))
val_cands = subset(val_cands, top_pcawg_val == "YES") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 

#Combined into one dataframe because need to get ranks 
all <- merge(rna, pcg, by = c("patient", "Cancer"))
all = all[,1:25170]

#------------------------------------------------------------------
#Within each tissue type, rank lncRNAs by which percentile of 
#expression they fall into to then compare with PCAWG lncRNAS exp
#------------------------------------------------------------------

#1. log1p
z = which(str_detect(colnames(all), "ENSG"))	
#all[,z] <- log1p(all[,z])

#2. Get lncRNA - median within each tissue type
allCands$combo = paste(allCands$gene, allCands$cancer, sep="_")
combos = unique(allCands$combo)

#3. Want ranking seperatley for high lncRNA expression group versus low lncRNA expression group
#---------------------------------------------------------
#Function 1
#for each lnc-cancer, label patient as lncRNA-risk or non-risk 
#---------------------------------------------------------

get_lnc_canc = function(comb){

	#plot ranks of lncRNA across samples 
	lnc = unlist(strsplit(comb, "_"))[1]
	canc = unlist(strsplit(comb, "_"))[2]
	lnc_dat = subset(all, Cancer == canc)
	pats = as.list(unique(lnc_dat$patient))
	
	get_ranks = function(pat){
		#order all genes by FPKM-UQ for each patient
		pat_dat = t(subset(lnc_dat, patient == pat))
		pat_dat = as.data.frame(pat_dat[3:nrow(pat_dat),])
		colnames(pat_dat)[1] = "FPKMUQ"
		pat_dat$gene = rownames(pat_dat)
		pat_dat = as.data.table(pat_dat)
		z = which(str_detect(pat_dat$gene, "ENSG"))
		pat_dat = pat_dat[z,]
		pat_dat$FPKMUQ = as.numeric(pat_dat$FPKMUQ)
		pat_dat = pat_dat[order(FPKMUQ)]
		pat_dat$rank = 1:nrow(pat_dat)
		pat_dat$score = pat_dat$rank/nrow(pat_dat)
		#return lncRNAs FPKM-UQ, rank and score
		lnc_pat = filter(pat_dat, gene == lnc)
		lnc_pat$canc = canc
		return(lnc)
	}
	lnc_ranks = llply(pats, get_ranks, .progress="text")
	lnc_ranks = ldply(lnc_ranks)
	return(lnc_ranks)
}

all_pats_lncs_ranks = llply(combos, get_lnc_canc, .progress="text")
saveRDS(all_pats_lncs_ranks, file="all_lncRNA_cancers_ranks.rds")





