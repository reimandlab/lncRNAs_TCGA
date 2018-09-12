###now have methylation files for each cancer 
###with probes that overlap any of the lncRNA candidates (might not be the cancer specific candidate)

library(data.table)
library(plyr)
library(dplyr)
library(stringr)
library(ggpubr)
library(ggExtra)
library(cowplot)
library(ggthemes)
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
library(parallel)
library(limma)
library(EnvStats)

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
colnames(fantom)[1] = "gene"

#############################################################
#----EDITED SEPT 11TH KI-------------------------------------
#############################################################


mypal = pal_npg("nrc", alpha = 0.7)(10)

tss_codes = read.csv("TCGA_TissueSourceSite_Codes2017.csv")
source_codes = source = read.csv("TCGA_sample_codes.csv")

#1. cands 
cands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
colnames(cands)[7] = "canc"

#2. lncRNAs that intersected with probes
probes = fread("fantom_lncrnas_mapped_to450_probes.bed")
colnames(probes) = c("probe_chr", "cpg_start", "cpg_end", "cgid", "rand", "cpgstrand", "lncchr", "lncstart", "lncend", "rand2", "ensg", "lncstrand", "lncname")

##add cancer type to which that gene belong to 
add_canc = function(probeid){
	lnc = probes$ensg[which(probes$cgid == probeid)]
	canc = cands$canc[which(cands$gene == lnc)][1]
	return(canc)
}
probes$canc = llply(probes$cgid, add_canc, .progress="text")
z = which(is.na(probes$canc))
if(!(length(z)==0)){
probes = probes[-z,]}
probes$combo = paste(probes$ensg, probes$canc, sep="_")
z = which(probes$combo %in% cands$combo)
length(unique(probes$combo[z])) #88 combos have methylation data 

#3. methylation files
#KIRC
KIRC = readRDS("KIRC_methylation_data_lncs_cands.rds")
#LIHC
LIHC = readRDS("LIHC_methylation_data_lncs_cands.rds")
#OV
OV = readRDS("OV_methylation_data_lncs_cands.rds")
#PAAD
PAAD = readRDS("PAAD_methylation_data_lncs_cands.rds")
#LUAD
LUAD = readRDS("LUADmethylation_data_lncs_cands.rds")
#BRCA
BRCA = readRDS("BRCAmethylation_data_lncs_cands.rds")
#
ACC = readRDS(file="ACCmethylation_data_lncs_cands.rds")

ESCA = readRDS(file="ESCAmethylation_data_lncs_cands.rds")

GBM = readRDS(file="GBMmethylation_data_lncs_cands.rds")

HNSC = readRDS(file="HNSCmethylation_data_lncs_cands.rds")

LGG = readRDS(file="LGGmethylation_data_lncs_cands.rds")

MESO = readRDS(file="MESOmethylation_data_lncs_cands.rds")

READ = readRDS(file="READmethylation_data_lncs_cands.rds")

THCA = readRDS(file="THCAmethylation_data_lncs_cands.rds")

UVM = readRDS(file="UVMmethylation_data_lncs_cands.rds")

CESC = readRDS(file="CESCmethylation_data_lncs_cands.rds")

KIRP = readRDS(file="KIRPmethylation_data_lncs_cands.rds")

STAD = readRDS(file="STADmethylation_data_lncs_cands.rds")

LUSC = readRDS(file="LUSCmethylation_data_lncs_cands.rds")

BLCA =  readRDS(file="BLCAmethylation_data_lncs_cands.rds")

COAD = readRDS(file="COADmethylation_data_lncs_cands.rds")

UCEC = readRDS(file="UCECmethylation_data_lncs_cands.rds")

SARC = readRDS(file="SARCmethylation_data_lncs_cands.rds")

#combine all 
methylation_data = list(KIRC,
LIHC,
OV,
PAAD,
LUAD,
BRCA,
ACC,
ESCA,
GBM,
HNSC,
LGG,
MESO,
READ,
THCA,
UVM,
CESC,
KIRP,
STAD,
LUSC,
BLCA,
COAD,
UCEC,
SARC)

cancers_order = c("KIRC",
"LIHC",
"OV",
"PAAD",
"LUAD",
"BRCA",
"ACC",
"ESCA",
"GBM",
"HNSC",
"LGG",
"MESO",
"READ",
"THCA",
"UVM",
"CESC",
"KIRP",
"STAD",
"LUSC",
"BLCA",
"COAD",
"UCEC",
"SARC")

cands$canc[cands$canc=="Ovarian serous cystadenocarcinoma"] = "OV"
cands$canc[cands$canc=="Liver hepatocellular carcinoma"] = "LIHC"
cands$canc[cands$canc=="Kidney renal clear cell carcinoma"] = "KIRC"
cands$canc[cands$canc=="Lung squamous cell carcinoma"] = "LUSC"
cands$canc[cands$canc=="Breast invasive carcinoma"] = "BRCA"
cands$canc[cands$canc=="Lung adenocarcinoma"] = "LUAD"
cands$canc[cands$canc=="Pancreatic adenocarcinoma"] = "PAAD"
cands$canc[cands$canc=="Adrenocortical carcinoma"] = "ACC"
cands$canc[cands$canc=="Bladder Urothelial Carcinoma"] = "BLCA"
cands$canc[cands$canc=="Stomach adenocarcinoma"] = "STAD"
cands$canc[cands$canc=="Head and Neck squamous cell carcinoma"] = "HNSC"
cands$canc[cands$canc=="Brain Lower Grade Glioma"] = "LGG"
cands$canc[cands$canc=="Sarcoma"] = "SARC"
cands$canc[cands$canc=="Kidney renal papillary cell carcinoma"] = "KIRP"
cands$canc[cands$canc=="Mesothelioma"] = "MESO"
cands$canc[cands$canc=="Uterine Corpus Endometrial Carcinoma"] = "UCEC"
cands$canc[cands$canc=="Uveal Melanoma"] = "UVM"
cands$canc[cands$canc=="Cervical squamous cell carcinoma and endocervical adenocarcinoma"] = "CESC"
cands$canc[cands$canc=="Colon adenocarcinoma"] = "COAD"
cands$canc[cands$canc=="Rectum adenocarcinoma"] = "READ"
cands$canc[cands$canc=="Thyroid carcinoma"] = "THCA"
cands$canc[cands$canc=="Glioblastoma multiforme"] = "GBM"
cands$canc[cands$canc=="Esophageal carcinoma"] = "ESCA"


probes$canc[probes$canc=="Ovarian serous cystadenocarcinoma"] = "OV"
probes$canc[probes$canc=="Liver hepatocellular carcinoma"] = "LIHC"
probes$canc[probes$canc=="Kidney renal clear cell carcinoma"] = "KIRC"
probes$canc[probes$canc=="Lung squamous cell carcinoma"] = "LUSC"
probes$canc[probes$canc=="Breast invasive carcinoma"] = "BRCA"
probes$canc[probes$canc=="Lung adenocarcinoma"] = "LUAD"
probes$canc[probes$canc=="Pancreatic adenocarcinoma"] = "PAAD"
probes$canc[probes$canc=="Adrenocortical carcinoma"] = "ACC"
probes$canc[probes$canc=="Bladder Urothelial Carcinoma"] = "BLCA"
probes$canc[probes$canc=="Stomach adenocarcinoma"] = "STAD"
probes$canc[probes$canc=="Head and Neck squamous cell carcinoma"] = "HNSC"
probes$canc[probes$canc=="Brain Lower Grade Glioma"] = "LGG"
probes$canc[probes$canc=="Sarcoma"] = "SARC"
probes$canc[probes$canc=="Kidney renal papillary cell carcinoma"] = "KIRP"
probes$canc[probes$canc=="Mesothelioma"] = "MESO"
probes$canc[probes$canc=="Uterine Corpus Endometrial Carcinoma"] = "UCEC"
probes$canc[probes$canc=="Uveal Melanoma"] = "UVM"
probes$canc[probes$canc=="Cervical squamous cell carcinoma and endocervical adenocarcinoma"] = "CESC"
probes$canc[probes$canc=="Colon adenocarcinoma"] = "COAD"
probes$canc[probes$canc=="Rectum adenocarcinoma"] = "READ"
probes$canc[probes$canc=="Thyroid carcinoma"] = "THCA"
probes$canc[probes$canc=="Glioblastoma multiforme"] = "GBM"
probes$canc[probes$canc=="Esophageal carcinoma"] = "ESCA"


#4. Expression data 
rna = readRDS("5919_lncRNAs_tcga_all_cancers_March13_wclinical_data.rds")
unique(rna$type)

get_data = function(lnc){
	comb = lnc

  print(lnc)
  cancer = cands$canc[which(cands$combo == lnc)][1]
	lnc = unlist(strsplit(lnc, "_"))[1]
  dat = dplyr::filter(probes, canc == cancer, ensg == lnc)
	
	if(!(dim(dat)[1]==0)){
	z = which(cancers_order == cancer)
	meth_dat = methylation_data[[z]]
	z = which(meth_dat$probe %in% dat$cgid)
	meth_dat = meth_dat[z,]

	#break up patients into individual rows
	#also only want tumours 
	pats = as.list(unique(colnames(meth_dat)[2:ncol(meth_dat)]))
	rearrange = function(pat){
		z = which(colnames(meth_dat)==pat)
		newdat = meth_dat[,c(1,z), with=FALSE]
		colnames(newdat) = c("probe", "beta", "gene", "post", "coord")
		newdat$patient = pat
		return(newdat)
	}
	newdats = llply(pats, rearrange, .progress="text")
	dat <- ldply(newdats, data.table)
	dat$gene = lnc
	z = which(is.na(dat$beta))
	if(!(length(z)==0)){
		dat = dat[-z,]
	}

	if(!(dim(dat)[1] == 0)){

	dat$source = ""
	get_source = function(id){
		source = unlist(strsplit(as.character(id), '-'))[4]
		source = substr(source, 1, 2)
		return(source)
		}
	dat$source = llply(dat$patient, get_source, .progress="text")
	#only keep primary solid tumour samples 
	dat = dplyr::filter(dat, source =="01")
	clean_tcga_id = function(id){
		s1 = unlist(strsplit(id, '-'))[1]
		s2 = unlist(strsplit(id, '-'))[2]
		s3 = unlist(strsplit(id, '-'))[3]
		return(paste(s1, s2, s3, sep="-"))
		}	
	dat$patient = llply(dat$patient, clean_tcga_id, .progress="text")

	#exp_data = rna[which(rna$patient %in% dat$patient), ]
	exp_data = rna[which(rna$type %in% cancer),]
	#assign high or low to each patient in expression file
	z <- which(colnames(exp_data) %in% lnc)
  	if(!(length(z)==0)){
  	df = as.data.frame(exp_data)
  	df <- df[,c(1, z,(ncol(exp_data)-30):ncol(exp_data))]  

	df$median <- ""
 	median2 <- quantile(as.numeric(df[,2]), 0.5)
  	#if(median2 ==0){
    #median2 = mean(as.numeric(df[,2]))
  	#}

  	 if(median2 ==0){
		    #if median = 0 then anyone greater than zero is 1 
		    l1 = which(df[,2] > 0)
		    l2 = which(df[,2] ==0)
		    df$median[l1] = 1
		    df$median[l2] = 0
		    }

	  if(!(median2 ==0)){
		    l1 = which(df[,2] >= median2)
		    l2 = which(df[,2] < median2)
		    df$median[l1] = 1
		    df$median[l2] = 0
		}


  	gene <- colnames(df)[2]
  	df$OS <- as.numeric(df$OS)
  	df$OS.time <- as.numeric(df$OS.time)
  	df$median[df$median ==0] = "Low"
  	df$median[df$median==1] = "High"
  	df$median = factor(df$median, levels=c("Low", "High"))

  	#get summary of SCNA in lncRNA for each patient 
  	#take mean segment mean for all cnas in patient 
  	#covering that lncRNA
  	colnames(df)[2] = "geneExp"
  	df = merge(df, dat, by=c("patient"))
  	
    #is high expression or low expression associated with worse prognosis? 
    df$median <- factor(df$median, levels = c("Low", "High"))
    df$median  # notice the changed order of factor levels
  	
	dim1 = table(df$median)[1]
	dim2 = table(df$median)[2]
	both = (dim1>=1) & (dim2 >=1)

	if(both & (length(unique(df$patient)) >=10)){

	HR = summary(coxph(Surv(OS.time, OS) ~ median, data = df))$coefficients[2]	
    
    if(HR >1){
      risk = "high_expression"
    }
    if(HR <1){
      risk = "low_expression"
    }

    df$risk = risk

    #is copy number aberation associated with expression? 
  	df$geneExp = log1p(df$geneExp)

  	library("ggExtra")
  	name = probes$lncname[which(probes$ensg == lnc)][1]
	   z = which(is.na(df$beta))
	   if(length(z) >=1){
		    df = df[-z,]
	   }
	   z =length(table(df$probe))	

  print(paste((z), "number of probes"))

	#multiple probes present 	
	if(z > 1){
		results_all_probes = as.data.frame(matrix(ncol = 8)) ; colnames(results_all_probes) = c("patient", "median", "probe", "mut_status", "gene", "cancer", "risk", "name")
	
	for(k in 1:z){
	new = subset(df, df$probe == unique(df$probe)[k])
	new$geneExp = as.numeric(new$geneExp)
	new$beta = as.numeric(new$beta)

	#everything else the same as if it was just one probe
	#get wilcoxon p-value stored between low and high exp patients - get avg beta value for each group 
		wilcoxon_pval = wilcox.test(beta ~ median, data =new)$p.value  
    fc = median(new$beta[new$median == "Low"]) - median(new$beta[new$median == "High"])
    #methylation in low expression group relative to high expression group

    #label probe as methylated or not methylated for each patient 
    new$mut_status = ""
        get_mut_stat = function(beta){
          meth = beta > 0.5
          unmeth = beta < (0.5)
          if(meth){
            return("Methylated")
             }        
           if(unmeth){
            return("Unmethylated")
           } 
           if((!(meth)) & (!(unmeth))){
            return("oneMeth")
          }
      }
      new$mut_status = as.character(llply(new$beta, get_mut_stat))

      #get how many people have low/high expression and methylated/unmethylated
      new$median = as.character(new$median)
      if(HR >1){
        med_risk = "High"
      }
      if(HR <1){
        med_risk = "Low"
      }

       r = rcorr(new$beta[new$median == med_risk], new$geneExp[new$median == med_risk], type="spearman")$r[2]
       rr = r #correlation in high risk group

       rp = rcorr(new$beta[!(new$median == med_risk)], new$geneExp[!(new$median == med_risk)], type="spearman")$r[2]
       rp = rp #correlation in low risk group
        if(is.na(rp)){
            rp = 0
          }
       #overall correlation
       ro = rcorr(new$beta, new$geneExp, type="spearman")$r[2]


    #get_chisq_pval
    t = table(new$median, new$mut_status)
    chi = chisq.test(t)
    chi_pval = tidy(chi)[2]

    #figure out if they are balanced or non-balanced 
    bal = cands[cands$data=="TCGA",]
    bal = as.numeric(bal$perc_risk[bal$combo == comb])

    length_risk_pats = length(which(new$median == med_risk))

    probe = new$probe[1]
    new$risk = med_risk
    gene_name = fantom$CAT_geneName[which(fantom$gene == lnc)]   
    new$name = gene_name
    new$cancer = cancer
    pat_dat = new[,c("patient", "median", "probe", "mut_status", "gene", "cancer", "risk", "name")]

    print(lnc)
    results_all_probes = rbind(results_all_probes, pat_dat)

	}

	results_all_probes = results_all_probes[-1,]
	results = results_all_probes
	}

	#only one probe present 
	if(z==1){

		#get wilcoxon p-value stored between low and high exp patients - get avg beta value for each group 
		wilcoxon_pval = wilcox.test(beta ~ median, data =df)$p.value  
		fc = median(df$beta[df$median == "Low"]) - median(df$beta[df$median == "High"])
    #methylation in low expression group relative to high expression group

    #label probe as methylated or not methylated for each patient 
		df$mut_status = ""
    		get_mut_stat = function(beta){
       		meth = beta > 0.5
       		unmeth = beta < (0.5)
       		if(meth){
      		  return("Methylated")
      			 }        
      		 if(unmeth){
      		  return("Unmethylated")
      		 } 
      		 if((!(meth)) & (!(unmeth))){
        		return("oneMeth")
       		}
    	}
   		df$mut_status = as.character(llply(df$beta, get_mut_stat))

   		#get how many people have low/high expression and methylated/unmethylated
   		df$median = as.character(df$median)
     	if(HR >1){
        med_risk = "High"
      }
      if(HR <1){
        med_risk = "Low"
      }

   		 r = rcorr(df$beta[df$median == med_risk], df$geneExp[df$median == med_risk], type="spearman")$r[2]
   		 rr = r #correlation in high risk group

   		 rp = rcorr(df$beta[!(df$median == med_risk)], df$geneExp[!(df$median == med_risk)], type="spearman")$r[2]
   		 rp = rp #correlation in low risk group
    		if(is.na(rp)){
    	  		rp = 0
    			}
       #overall correlation
       ro = rcorr(df$beta, df$geneExp, type="spearman")$r[2]


    #get_chisq_pval
    t = table(df$median, df$mut_status)
    chi = chisq.test(t)
    chi_pval = tidy(chi)[2]

    #figure out if they are balanced or non-balanced 
    bal = cands[cands$data=="TCGA",]
    bal = as.numeric(bal$perc_risk[bal$combo == comb])

    length_risk_pats = length(which(df$median == med_risk))

    probe = df$probe[1]
	  df$risk = med_risk
    gene_name = fantom$CAT_geneName[which(fantom$gene == lnc)]   
    df$name = gene_name
    df$cancer = cancer
    pat_dat = df[,c("patient", "median", "probe", "mut_status", "gene", "cancer", "risk", "name")]
  	results = pat_dat
    print(lnc)

  }
	}

	if(!((both & (length(unique(df$patient)) >=10)))){
		results = as.data.frame(matrix(ncol = 8)) ; colnames(results) = c("patient", "median", "probe", "mut_status", "gene", "cancer", "risk", "name")
  }  
  return(results)
}
}
}
}

genes = as.list(unique(as.character(cands$combo[which(cands$combo %in% probes$combo)]))) #88/166 have methylation probes overlapping them 
lnc_meth_cancer_data = llply(genes, get_data, .progress="text")

lnc_meth_cancer_data2 = Filter(Negate(is.null), lnc_meth_cancer_data)
lnc_meth_cancer_data2 = ldply(lnc_meth_cancer_data2)
lnc_meth_cancer_data2 = as.data.table(lnc_meth_cancer_data2)
z= which(is.na(lnc_meth_cancer_data2$cancer))
lnc_meth_cancer_data2 = lnc_meth_cancer_data2[-z,] #362 lncRNA-probe pairs evaluated 


saveRDS(lnc_meth_cancer_data2, file="all_dfs_methylaion_status_exp_patients_sept12.rds")










