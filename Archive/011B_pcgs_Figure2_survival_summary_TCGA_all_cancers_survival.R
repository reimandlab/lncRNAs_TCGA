###---------------------------------------------------------------
###Load libraries and data - April 16th 
###---------------------------------------------------------------

source("source_code_Cox_MonteCarlo_CV_April12.R")
require(caTools)

library(survAUC)
#source("source_code_Cox_MonteCarlo_CV_Mar13.R")
require(caTools)
#check if this person is in my analysis: TCGA-61-2095
library(glmnet)
library(survcomp)
library(caret)
library(stringr)
library(EnvStats)
library(patchwork)

source("universal_LASSO_survival_script.R")

library(ggpubr)
library(ggrepel)
library(viridis)
library(patchwork)
library(caret)  
library(Rtsne)
library(data.table)

#------DATA---------------------------------------------------------
#UCSC gene info
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
#z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
#ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,8]))
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

#save RNA and PCG files locally
#saveRDS(rna, file="rna_lncRNAs_expression_data_june29.rds")
#saveRDS(pcg, file="rna_pcg_expression_data_june29.rds")

rna = readRDS("rna_lncRNAs_expression_data_june29.rds")
pcg = readRDS("rna_pcg_expression_data_june29.rds")

#------FEATURES-----------------------------------------------------

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = subset(allCands, data == "TCGA") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 
allCands$combo = unique(paste(allCands$gene, allCands$cancer, sep="_"))
cands_dups = unique(allCands$gene[which(duplicated(allCands$gene))])

val_cands = read.csv("175_lncRNA_cancers_combos_23_cancer_types_july5.csv")
val_cands = as.data.table(val_cands)
val_cands = subset(val_cands, data == "PCAWG") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 
val_cands$combo = unique(paste(val_cands$gene, val_cands$cancer, sep="_"))
val_cands = subset(val_cands, top_pcawg_val == "YES") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 


#start with only lncRNA_intergenic
#lincs = subset(fantom, CAT_geneClass == "lncRNA_intergenic")
#z = which(colnames(rna) %in% lincs$gene)
#rna = as.data.frame(rna)
#rna = rna[,c(z, (ncol(rna)-5):ncol(rna))]

###[2.] Data splitting 

###---------------------------------------------------------------
###---------------------------------------------------------------

#check if this person is in my analysis: TCGA-61-2095
library(glmnet)
library(survcomp)
library(caret)
library(stringr)
library(factoextra)
library(patchwork)
library(tidyr)

rna = as.data.frame(rna)
pcg = as.data.frame(pcg)

dim(rna)
dim(pcg)
dim(norm)
dim(met)

###---------------------------------------------------------------

#function that tests each lncRNA's survival 

#1. remove discrepancy 
z = which(pcg$vital_status == "[Discrepancy]")
pcg = pcg[-z,]

#2. list of cancers to apply function to 
cancers = as.list(unique(pcg$Cancer))

#remove cancer types with less than 50 patients 
pats_num = as.data.table(table(pcg$Cancer))
pats_num = filter(pats_num, N <50)
canc_rm = pats_num$V1

#remove those ones
cancers = cancers[which(!(cancers %in% canc_rm))]

#3. function that splits data into cancers 
get_canc = function(canc){
	canc_data = pcg[which(pcg$Cancer == canc),]
	return(canc_data)
}

canc_datas = llply(cancers, get_canc)

#4. function that calculates survival for each gene 
#det_lncs = readRDS("all_TCGA_cancers_lncRNAs_detectable_May18.rds")
#det_lncs =filter(det_lncs, status =="detectable")

canc_survival_genes = function(dato){
	#look at all lncRNAs that are expressed in at least some patients 
  z = which(str_detect(colnames(dato), "ENSG"))
  sums = apply(dato[,z], 2, sum)
  rm = names(sums[which(sums == 0)])
  z = which(colnames(dato) %in% rm)
  dato = dato[,-z]

  print(dato$type[1])

  z = which(str_detect(colnames(dato), "ENSG"))
  genes = unique(colnames(dato)[z])	

  #TEST------------------------------------------------------------------------------------------
  #genes = genes[1:100]
	canc_data_genes_analyze = dato 
	
	get_survival = function(gene){
	  print(gene)
  	results_cox <- as.data.frame(matrix(ncol=8)) ; colnames(results_cox) <- c("gene", "coef", "HR", "pval", "low95", "upper95", "risk_size", "num_patients")
  	z = which(colnames(canc_data_genes_analyze) == gene)
  	dat = canc_data_genes_analyze[,c(1,z,(ncol(canc_data_genes_analyze)-33):ncol(canc_data_genes_analyze))]
  	dat$OS.time = as.numeric(dat$OS.time)
  	dat$OS = as.numeric(dat$OS)
  	#remove NAs
  	z = which(is.na(dat$OS.time))
  	if(!(length(z) ==0)){
  	dat = dat[-z,]}
	  med_gene = median(as.numeric(dat[,2]))  	
	  dat$med = ""
	  if(med_gene ==0){
        #if median = 0 then anyone greater than zero is 1 
        l1 = which(dat[,2] > 0)
        l2 = which(dat[,2] ==0)
        dat$med[l1] = 1
        dat$med[l2] = 0
        }

    if(!(med_gene ==0)){
        l1 = which(dat[,2] >= med_gene)
        l2 = which(dat[,2] < med_gene)
        dat$med[l1] = 1
        dat$med[l2] = 0
    }

    check1 = table(dat$med)[1] >= 10
    check2 = table(dat$med)[2] >= 10
    
    if(!(dim(table(dat$med))==1)){

    if(check1 & check2){
      if(dim(table(dat$med)) ==2){
  	  res.cox <- coxph(Surv(OS.time, OS) ~ med, data = dat)
    	hr = summary(res.cox)$coefficients[1,c(2)]
      num_pat = nrow(dat)
      if(hr > 1){
        risk = length(which(dat$med ==1))
      }
      if(hr <1){
        risk = length(which(dat$med ==0))
      }

      row <- c(gene, summary(res.cox)$coefficients[1,c(1,2,5)],  summary(res.cox)$conf.int[1,c(3,4)], risk, num_pat)
     	names(row) <- names(results_cox)
    	return(row)
  	}}}} #end get_survival function

	genes_survival = llply(genes, get_survival, .progress="text")
	genes_survival_res = ldply(genes_survival, rbind)
	#fdr
	colnames(genes_survival_res) = c("gene", "coef", "HR", "pval", "low95", "upper95", "risk_size", "num_patients")
	genes_survival_res$fdr = p.adjust(as.numeric(genes_survival_res$pval), method="fdr")
	genes_survival_res$canc = dato$Cancer[1]
	genes_survival_res = as.data.table(genes_survival_res)
	genes_survival_res = genes_survival_res[order(fdr)]
	return(genes_survival_res)
}

all_cancers_genes_surv = llply(canc_datas, canc_survival_genes, .progress="text")
all_cancers_genes_surv_comb = ldply(all_cancers_genes_surv, data.frame)

saveRDS(all_cancers_genes_surv_comb, file="mRNAs_Survival_Results_prognostic_pcgs_July19.rds")

















