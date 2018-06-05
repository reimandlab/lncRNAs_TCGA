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


library(DMRcate)

#all probes coordinates
probe_coordinates = fread("hm450.hg19.manifest")
#strand is 6th column 
probe_coordinates = probe_coordinates[,c(1:3, 5,9,4)]


mypal = pal_npg("nrc", alpha = 0.7)(10)

tss_codes = read.csv("TCGA_TissueSourceSite_Codes2017.csv")
source_codes = source = read.csv("TCGA_sample_codes.csv")

#1. cands 
cands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_May4.rds")
cands = filter(cands, AnalysisType == "noFDR")
colnames(cands)[7] = "canc"

#3. methylation files
#KIRC
KIRC = readRDS("KIRC_methylation_data_ALLprobes.rds")
#LIHC
LIHC = readRDS("LIHC_methylation_data_ALLprobes.rds")
#OV
OV = readRDS("OV_methylation_data_ALLprobes.rds")
#PAAD
PAAD = readRDS("PAAD_methylation_data_ALLprobes.rds")
#LUAD
LUAD = readRDS("LUADmethylation_data_ALLprobes.rds")
#BRCA
BRCA = readRDS("BRCAmethylation_data_ALLprobes.rds")
#
ACC = readRDS(file="ACCmethylation_data_ALLprobes.rds")

ESCA = readRDS(file="ESCAmethylation_data_ALLprobes.rds")

GBM = readRDS(file="GBMmethylation_data_ALLprobes.rds")

HNSC = readRDS(file="HNSCmethylation_data_ALLprobes.rds")

LGG = readRDS(file="LGGmethylation_data_ALLprobes.rds")

MESO = readRDS(file="MESOmethylation_data_ALLprobes.rds")

READ = readRDS(file="READmethylation_data_ALLprobes.rds")

THCA = readRDS(file="THCAmethylation_data_ALLprobes.rds")

UVM = readRDS(file="UVMmethylation_data_ALLprobes.rds")

CESC = readRDS(file="CESCmethylation_data_ALLprobes.rds")

KIRP = readRDS(file="KIRPmethylation_data_ALLprobes.rds")

STAD = readRDS(file="STADmethylation_data_ALLprobes.rds")

LUSC = readRDS(file="LUSCmethylation_data_ALLprobes.rds")

BLCA =  readRDS(file="BLCAmethylation_data_ALLprobes.rds")

COAD = readRDS(file="COADmethylation_data_lncs_cands_ALLprobes.rds")

UCEC = readRDS(file="UCECmethylation_data_ALLprobes.rds")

SARC = readRDS(file="SARCmethylation_data_ALLprobes.rds")

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


#4. Expression data 
rna = readRDS("5919_lncRNAs_tcga_all_cancers_March13_wclinical_data.rds")
unique(rna$type)


get_data = function(lnc){
	cancer = cands$canc[which(cands$gene == lnc)][1]
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

	#multiple probes present 	
	if(z > 1){
		results_all_probes = as.data.frame(matrix(ncol = 12)) ; colnames(results_all_probes) = c("cancer", "gene", "num_patients", "numHighexpMethmatch", "numLowexpMethmatch", 
			"wilcoxon_pval", "risk_type", "num_risk_pats", "num_risk_pats_wmatchingMeth", 
    	"risk_group_correlation", "nonrisk_group_correlation", "probe")
	
	for(k in 1:z){
	new = subset(df, df$probe == unique(df$probe)[k])
	new$geneExp = as.numeric(new$geneExp)
	new$beta = as.numeric(new$beta)

	#everything else the same as if it was just one probe
	#get wilcoxon p-value stored between low and high exp patients - get avg beta value for each group 
		wilcoxon_pval = wilcox.test(beta ~ median, data =new)$p.value  
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
     	new$exp_meth_match = ""
    	get_meth_exp = function(row){
    	  exp = row[[34]]
    	  meth = row[[42]]
    	  if(((exp == "High") & (meth == "Unmethylated"))){
    	    return("HighExp_NoMeth")
    	  }
    	    if(((exp == "Low") & (meth == "Methylated"))){
    	    return("LowExp_Meth")
    	  }
    	    else{
       	   return("DontMatch")
       	 }
    		}
    	new$exp_meth_match = as.character(apply(new, 1, get_meth_exp))

    	#get if risk matches exp 
    	if(new$risk[1] == "high_expression"){
    	  med_risk = "High"
    	  meth_risk = "HighExp_NoMeth"
   		 }

   		 if(new$risk[1] == "low_expression"){
    	  med_risk = "Low"
    	  meth_risk = "LowExp_Meth"
   		 }

   		 r = rcorr(new$beta[which(new$median == med_risk)], new$geneExp[which(new$median == med_risk)], type="spearman")$r[2]
   		 rr = r #correlation in high risk group

   		 rp = rcorr(new$beta[which(new$median != med_risk)], new$geneExp[which(new$median != med_risk)], type="spearman")$r[2]
   		 rp = rp #correlation in low risk group
    		if(is.na(rp)){
    	  		rp = 0
    			}


    length_risk_pats = length(which(new$median == med_risk))
    risk_pats = new[which(new$median == med_risk),]

    #what is the number of pateints wtih a cna that matches risk? (within risk grou)
    length_risk_pats_wmeth = length(which(risk_pats$exp_meth_match == meth_risk))
    probe = new$probe[k]
	results = c(cancer, unique(new$gene), length(unique(new$patient)), length(which(new$exp_meth_match == "HighExp_NoMeth")), 
      length(which(new$exp_meth_match == "LowExp_Meth")), wilcoxon_pval, risk, length_risk_pats, length_risk_pats_wmeth, rr, rp, probe)
    names(results) = c("cancer", "gene", "num_patients", "numHighexpMethmatch", "numLowexpMethmatch", "wilcoxon_pval", "risk_type", "num_risk_pats", "num_risk_pats_wmatchingMeth", 
    	"risk_group_correlation", "nonrisk_group_correlation", "probe")
    new$median <- factor(new$median, levels = c("Low", "High"))
    new$median  # notice the changed order of factor levels


    	#scatter plot
    	sp = ggscatter(new, 
	   	x = "beta", y = "geneExp", conf.int = TRUE, # Add confidence interval
   				cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   				cor.method = "spearman", 
   				cor.coeff.args = list(method = "spearman", label.y = 1, label.sep = "\n"), 
               color = "median", palette = mypal[c(4,1)], facet.by="median", 
               size = 3, alpha = 0.6, add = "reg.line",                         # Add regression line
          	   ggtheme = theme_light(), ylab = "log1p(FPKM) Expression", font.x = c(15, "plain", "black"), font.y = c(15, "plain", "black"),
          font.tickslab = c(15, "plain", "black"), xlab="Beta value", title = paste(cancer, probe, gene, "Expression vs Methylation", "\npats wHigh Exp=", length(which(new$median=="High")),
				"\npats wLow Exp=" , length(which(new$median=="Low")), "\nwilcoxon p-val =", round(wilcoxon_pval, digits=5), "\nexpression risk = ", med_risk))
    	print(sp)


    	#density plot
   		p =  ggdensity(new, x = "beta", color = "median", fill="median", palette= mypal[c(4,1)], alpha=0.25, xlab="DNA Methylation",
			title = paste(cancer, probe, gene, "Expression vs Methylation", "\npats wHigh Exp=", length(which(new$median=="High")),
				"\npats wLow Exp=" , length(which(new$median=="Low")), "\nwilcoxon p-val =", round(wilcoxon_pval, digits=5), "\nexpression risk = ", med_risk))
   		p = p + geom_vline(xintercept = c(0, 0.5, 1), linetype="dotted", 
                color = "blue")
		print(p)


    results_all_probes = rbind(results_all_probes, results)


    print(lnc)
	}
	results_all_probes = results_all_probes[-1,]
	results = results_all_probes
	}

	#only one probe present 
	if(z==1){

		#get wilcoxon p-value stored between low and high exp patients - get avg beta value for each group 
		wilcoxon_pval = wilcox.test(beta ~ median, data =df)$p.value  
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
     	df$exp_meth_match = ""
    	get_meth_exp = function(row){
    	  exp = row[[34]]
    	  meth = row[[42]]
    	  if(((exp == "High") & (meth == "Unmethylated"))){
    	    return("HighExp_NoMeth")
    	  }
    	    if(((exp == "Low") & (meth == "Methylated"))){
    	    return("LowExp_Meth")
    	  }
    	    else{
       	   return("DontMatch")
       	 }
    		}
    	df$exp_meth_match = as.character(apply(df, 1, get_meth_exp))

    	#get if risk matches exp 
    	if(df$risk[1] == "high_expression"){
    	  med_risk = "High"
    	  meth_risk = "HighExp_NoMeth"
   		 }

   		 if(df$risk[1] == "low_expression"){
    	  med_risk = "Low"
    	  meth_risk = "LowExp_Meth"
   		 }

   		 r = rcorr(df$beta[df$median == med_risk], df$geneExp[df$median == med_risk], type="spearman")$r[2]
   		 rr = r #correlation in high risk group

   		 rp = rcorr(df$beta[!(df$median == meth_risk)], df$geneExp[!(df$median == meth_risk)], type="spearman")$r[2]
   		 rp = rp #correlation in low risk group
    		if(is.na(rp)){
    	  		rp = 0
    			}


    length_risk_pats = length(which(df$median == med_risk))
    risk_pats = df[which(df$median == med_risk),]

    #what is the number of pateints wtih a cna that matches risk? (within risk grou)
    length_risk_pats_wmeth = length(which(risk_pats$exp_meth_match == meth_risk))
    probe = df$probe[1]
	results = c(cancer, unique(df$gene), length(unique(df$patient)), length(which(df$exp_meth_match == "HighExp_NoMeth")), 
      length(which(df$exp_meth_match == "LowExp_Meth")), wilcoxon_pval, risk, length_risk_pats, length_risk_pats_wmeth, rr, rp, probe)
    names(results) = c("cancer", "gene", "num_patients", "numHighexpMethmatch", "numLowexpMethmatch", "wilcoxon_pval", "risk_type", "num_risk_pats", "num_risk_pats_wmatchingMeth", 
    	"risk_group_correlation", "nonrisk_group_correlation", "probe")
    df$median <- factor(df$median, levels = c("Low", "High"))
    df$median  # notice the changed order of factor levels


       	#scatter plot
    	 	sp = ggscatter(df, 
	   	x = "beta", y = "geneExp", conf.int = TRUE, # Add confidence interval
   				cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   				cor.method = "spearman", 
   				cor.coeff.args = list(method = "spearman", label.y = 1, label.sep = "\n"), 
               color = "median", palette = mypal[c(4,1)], facet.by="median", 
               size = 3, alpha = 0.6, add = "reg.line",                         # Add regression line
          	   ggtheme = theme_light(), ylab = "log1p(FPKM) Expression", font.x = c(15, "plain", "black"), font.y = c(15, "plain", "black"),
          font.tickslab = c(15, "plain", "black"), xlab="Beta value", title = paste(cancer, probe, gene, "Expression vs Methylation", "\npats wHigh Exp=", length(which(df$median=="High")),
				"\npats wLow Exp=" , length(which(df$median=="Low")), "\nwilcoxon p-val =", round(wilcoxon_pval, digits=5), "\nexpression risk = ", med_risk))
    	print(sp)


    	#density plot
   		p =  ggdensity(df, x = "beta", palette = mypal[c(4,1)], color = "median", fill="median", alpha=0.25, xlab="DNA Methylation",
			title = paste(cancer, probe, gene, "Expression vs Methylation", "\npats wHigh Exp=", length(which(df$median=="High")),
				"\npats wLow Exp=" , length(which(df$median=="Low")), "\nwilcoxon p-val =", round(wilcoxon_pval, digits=5), "\nexpression risk = ", med_risk))

   		p = p + geom_vline(xintercept = c(0, 0.5, 1), linetype="dotted", 
                color = "blue")
		print(p)
  		print(lnc)
	}
	}

	if(!((both & (length(unique(df$patient)) >=10)))){
		results = as.data.frame(matrix(ncol = 12)) ; colnames(results) = c("cancer", "gene", "num_patients", "numHighexpMethmatch", "numLowexpMethmatch", 
			"wilcoxon_pval", "risk_type", "num_risk_pats", "num_risk_pats_wmatchingMeth", 
    	"risk_group_correlation", "nonrisk_group_correlation", "probe")
	}

    return(results)
}
}
}
}



lnc_meth_cancer_data = llply(genes, get_data, .progress="text")





















































































#########################################################################################################################
#SUMMARIZE#
#########################################################################################################################

lnc_meth_cancer_data2 = as.data.frame(do.call("rbind", lnc_meth_cancer_data))
z = which(is.na(lnc_meth_cancer_data2[,1]))
lnc_meth_cancer_data2 = lnc_meth_cancer_data2[-z,]

lnc_meth_cancer_data2 = as.data.table(lnc_meth_cancer_data2)

lnc_meth_cancer_data2$wilcoxon_pval = as.numeric(as.character(lnc_meth_cancer_data2$wilcoxon_pval))
lnc_meth_cancer_data2 = lnc_meth_cancer_data2[order(wilcoxon_pval)]
lnc_meth_cancer_data2$fdr = p.adjust(lnc_meth_cancer_data2$wilcoxon_pval, method="fdr")

sig_diff = filter(lnc_meth_cancer_data2, wilcoxon_pval <=0.05)
sig_diff = as.data.frame(sig_diff)

#is probe in promoter or exon? 
colnames(probes)[11] = "gene"
probes = probes[,c(1:4, 6:9, 11:14)]
colnames(sig_diff)[12] = "cgid"
sig_diff = merge(sig_diff, probes, by=c("gene", "cgid"))

check_promoter = function(row){
	strand = row$lncstrand
	lnc_prom_start = row$lncstart-2000
	lnc_prom_end = row$lncstart+2000
	cpg = row$cpg_start+1

	#promoter
	if((cpg >lnc_prom_start) & (cpg<lnc_prom_end)){
		region="promoter"
	}

	if(!((cpg >lnc_prom_start) & (cpg<lnc_prom_end))){
		region="body"
	}
	return(region)
}
sig_diff$region = apply(sig_diff, 1, check_promoter)
sig_diff = as.data.frame(sig_diff)

#subset to promoters
sig_diff = subset(sig_diff, region == "promoter")

#these 15 genes have a significant difference in median beta values dsitributions between high 
#lncRNA expression and low lncRNA expression groups 
saveRDS(sig_diff, file="methylation_analysis_of_candidate_lncRNAs_promoter_only_wilcoxon_resutls_May31.rds")

#----------------------------------------------------
#SUMMARIZE and visualize same as copy number analysis 
#----------------------------------------------------

#if multiple probes get most correlated probe 
dups = as.data.table(table(sig_diff$gene), sig_diff$cgid)
dups = dups$V1[which(dups$N >1)]

#seperate sig_diffs

probes_keep = c()

for(i in 1:length(dups)){
	g = subset(sig_diff, gene == dups[i])
	#get risk group correlation
	cors = as.numeric(min((as.numeric(g$risk_group_correlation))))
	z = which(g$risk_group_correlation == cors)
	probe_keep = g$cgid[z]
	probes_keep = c(probes_keep, probe_keep)
}

dups_sig_diff = subset(sig_diff, gene %in% dups)
dups_sig_diff = subset(sig_diff, cgid %in% probes_keep)

sig_diff_nodups = subset(sig_diff, !(gene %in% dups))
sig_diff = rbind(dups_sig_diff, sig_diff_nodups)


#---> 15 lncRNAs in total with significant wilcoxon p-value 
#now let's look at the number of risk pateints that have matcinng methylation trends
#and how strong the correlation between beta value and expression is 


sig_diff$no_match = ""
get_nomatch = function(row){
  nomatch = as.numeric(as.character(row[[4]])) #all patients #
  high = as.numeric(as.character(row[[5]])) #of pats high match
  low = as.numeric(as.character(row[[6]])) 
  nomatch = nomatch-high-low
}
sig_diff$no_match = apply(sig_diff, 1, get_nomatch)

match_h = sig_diff[,c(1:5,7:24)]
colnames(match_h)[5] = "num_meth_exp_match"
match_h$type = "HighExpNoMeth"
match_h$num_meth_exp_match = as.numeric(as.character(match_h$num_meth_exp_match))

match_l = sig_diff[,c(1:4, 6, 7:24)]
colnames(match_l)[5] = "num_meth_exp_match"
match_l$type = "LowExpMeth"
match_l$num_meth_exp_match = as.numeric(as.character(match_l$num_meth_exp_match))

match_no = sig_diff[,c(1:4,25, 7:24)]
colnames(match_no)[5] = "num_meth_exp_match"
match_no$type = "NoExpMethMatch"

matched_sig = rbind(match_h, match_l, match_no)
matched_sig = as.data.table(matched_sig)
matched_sig = matched_sig[order(num_meth_exp_match, -type)]
matched_sig = as.data.frame(matched_sig)
order = as.character(unique(matched_sig$gene))
matched_sig$gene <- factor(matched_sig$gene, levels = order)
matched_sig$gene  # notice the changed order of factor levels

#for each lncRNA turn fractions into percentages 
matched_sig$num_meth_exp_match_patients = matched_sig$num_meth_exp_match
matched_sig$num_meth_exp_match = (as.numeric(matched_sig$num_meth_exp_match))/(as.numeric(as.character(matched_sig$num_patients)))
matched_sig$num_meth_exp_match = round(matched_sig$num_meth_exp_match, digits=3)

matched_sig = as.data.table(matched_sig)
matched_sig = matched_sig[order(risk_type, num_meth_exp_match, type)]

order = unique(matched_sig$gene)
matched_sig$gene <- factor(matched_sig$gene, levels = order)

matched_sig$num_meth_exp_match = as.numeric(matched_sig$num_meth_exp_match)
# Stacked bar plots, add labels inside bars
pdf("num_matching_Methylation_exp_tags_10cancers_15lncs.pdf", width= 8, height=5)
p1 = ggbarplot(matched_sig, x = "gene", y = "num_meth_exp_match",
  color = "type", fill = "risk_type", legend="left", ylab = "% of patients", xlab = "lncRNA",
  palette = mypal) + theme_light() +coord_flip() + scale_colour_manual(values = c("salmon", "steelblue", "snow3")) +
  scale_fill_manual(values = c("mistyrose", "aliceblue")) + ggtitle("Overall summary of CNAs overlapping lncRNA") +
   geom_hline(yintercept=0.5, linetype="dashed", color = "red")
  #label = TRUE, lab.col = "white", lab.pos = "in", lab.size=1.8)
p1 = ggpar(p1,
 font.tickslab = c(8,"plain", "black"),
 xtickslab.rt = 45, legend = "left") 
p1
dev.off()

#summarize where eahc gene cancer it's coming frmo and whether risk is high or low

#high risk
z = which((matched_sig$risk_type == "high_expression") & (matched_sig$type == "HighExpNoMeth"))
high = matched_sig[z,]

#low risk
z = which((matched_sig$risk_type == "low_expression") & (matched_sig$type == "LowExpMeth"))
low = matched_sig[z,]

pats_summary = rbind(high, low)
order = as.character(unique(matched_sig$gene))

pats_summary$per_patients_wrisk_meth = (as.numeric(as.character(pats_summary$num_risk_pats_wmatchingMeth)))/(as.numeric(as.character(pats_summary$num_risk_pats)))
pats_summary = as.data.table(pats_summary)

pats_summary = pats_summary[order(match(gene, order))] 

labels = as.character(pats_summary$cancer)

p2 = ggbarplot(pats_summary, x = "gene", y = "per_patients_wrisk_meth", order=order, ylab="% of risk \n patients wMethylation", xlab = "lncRNA", 
  color = "type", fill = "risk_type", legend="left", label=labels, lab.size = 2.5, lab.pos = c("out"), lab.vjust = 0.8, lab.hjust = 0.67, 
  palette = mypal) + theme_light() +coord_flip() + scale_colour_manual(values = c("salmon", "steelblue", "snow3")) +
  scale_fill_manual(values = c("mistyrose", "aliceblue")) + geom_hline(yintercept=0.2, linetype="dashed", color = "red") + 
  ggtitle("lncRNA Methylation in Risk Groups")
  #label = TRUE, lab.col = "white", lab.pos = "in", lab.size=1.8)
p2= ggpar(p2,
 font.tickslab = c(8,"plain", "black"),
 xtickslab.rt = 45, legend = "none") 

library(patchwork)
pdf("Methylation_lncRNAs_summary_15cands_may31.pdf", width=12, height=9)
p1 + p2 + plot_layout(ncol = 2, widths = c(2, 1))
dev.off()

#divide up the risk correlation and no risk correlation groups
cors_risk = pats_summary[,c(1:10, 12:26)]
cors_nonrisk = pats_summary[,c(1:9, 11, 12:26)]

colnames(cors_risk)[10] = "cor_spearman"
cors_risk$cor_group = "risk"
colnames(cors_nonrisk)[10] = "cor_spearman"
cors_nonrisk$cor_group = "nonrisk"

pats_summary = rbind(cors_risk, cors_nonrisk)

#plot the overall correlation between lncRNA copy number and expression in cohort 
pats_summary$cor_spearman = as.numeric(as.character(pats_summary$cor_spearman))

p3 = ggbarplot(pats_summary, x = "gene", y = "cor_spearman", order=order, ylab="Spearman correlation \n methylation wExpression", xlab = "lncRNA", 
  add = "segments", palette = mypal, color = "cor_group", position = position_dodge(0.6)) + theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +  theme_light() +
  coord_flip() + 
  geom_hline(yintercept=c(0), linetype="dashed", color = "black") + scale_colour_manual(values = c("steelblue", "salmon")) +
  ggtitle("lncRNA Methylation in Risk Groups") + geom_hline(yintercept=c(-0.15), linetype="dashed", color = "red") 
  #label = TRUE, lab.col = "white", lab.pos = "in", lab.size=1.8)
p3= ggpar(p3,
 font.tickslab = c(8,"plain", "black"),
 xtickslab.rt = 45, legend = "right") 

pdf("methylation_lncRNAs_summary_15cands_may31_wcor.pdf", width=14, height=8)
p1 + p2 + p3 + plot_layout(ncol = 3, widths = c(2, 1, 1))
dev.off()




























