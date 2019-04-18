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

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = subset(allCands, data == "TCGA") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 
allCands$combo = unique(paste(allCands$gene, allCands$cancer, sep="_"))

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


rna = readRDS("5919_lncRNAs_tcga_all_cancers_March13_wclinical_dataalldat.rds")
table(rna$type)

z = which(rna$vital_status == "[Discrepancy]")
rna = rna[-z,]
z = which(is.na(as.numeric(rna$age_at_initial_pathologic_diagnosis)))
rna = rna[-z,]
z = which(is.na(as.numeric(rna$OS.time)))
rna = rna[-z,]
z = which(as.numeric(rna$OS.time) == 0)
rna = rna[-z,]

table(rna$type)

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
     probecount = z
  print(paste((z), "number of probes"))

	#multiple probes present 	
	if(probecount > 1){
		results_all_probes = as.data.frame(matrix(ncol = 25)) ; colnames(results_all_probes) = c("cancer", "gene", "num_patients", "wilcoxon_pval", "chi_pval", "fc", "risk_type", "num_risk_pats",  
      "risk_group_correlation", "nonrisk_group_correlation", "overall_correlation", "rop", "probe", "risk_pat_bal", "stat_exp_cor", "stat_exp_pval", "ks_test", 
      "risk_met", "nonrisk_met", "other", "nonrisk_unmet", "risk_not_met", "risk_other", "nonrisk_other", "cox_p")
	
	for(k in 1:probecount){
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
          meth = beta >= 0.75
          unmeth = beta < (0.25)
          other = ((beta >= 0.25) & (beta < 0.75))
          if(meth){
            return("Methylated")
             }        
           if(unmeth){
            return("Unmethylated")
           } 
           if(other){
            return("other")
          }
      }
      new$mut_status = as.character(llply(new$beta, get_mut_stat))

      if(!(length(unique(new$mut_status))==1)){

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
       rop = rcorr(new$beta, new$geneExp, type="spearman")$P[2]

    #get_chisq_pval
    t = table(new$median, new$mut_status)
    print(t)
    chi = chisq.test(t)
    chi_pval = tidy(chi)[2]

    #figure out if they are balanced or non-balanced 
    bal = cands[cands$data=="TCGA",]
    bal = as.numeric(bal$perc_risk[bal$combo == comb])

    length_risk_pats = length(which(new$median == med_risk))

    probe = new$probe[1]
    
    new$median <- factor(new$median, levels = c("High", "Low"))
    new$median  # notice the changed order of factor levels

    print("pass2")

    library("ggExtra")
    new$risk = med_risk

    gene_name = fantom$CAT_geneName[which(fantom$gene == lnc)]   

    sp1 = ggplot(new, aes(x=beta, y=geneExp, color=median)) + ggtitle(paste(new$gene[1], gene_name, probe, cancer, "\nrisk=", new$risk[1])) + 
    geom_point() + 
    geom_smooth(method=lm, se=FALSE, fullrange=TRUE) + stat_cor() +
    geom_vline(xintercept=0.5, linetype="dashed", color = "grey") + 
    geom_vline(xintercept=0, linetype="dashed", color = "grey") + xlab("Beta Values") +
    ylab("log1p FPKM") + theme(text = element_text(size=15), axis.text = element_text(size=15))


    sp2 = ggplot(new, aes(x=beta, y=geneExp)) + ggtitle(paste(new$gene[1], gene_name, probe, cancer, "\nrisk=", new$risk[1])) + 
    geom_point() + 
    geom_smooth(method=lm, se=FALSE, fullrange=TRUE) + stat_cor() +
    geom_vline(xintercept=0.5, linetype="dashed", color = "grey") + 
    geom_vline(xintercept=0, linetype="dashed", color = "grey") + xlab("Beta Values") +
    ylab("log1p FPKM") + theme(text = element_text(size=15), axis.text = element_text(size=15))


    sp3 = ggplot(new, aes(x=median, y=beta, color=median)) + ggtitle(paste(gene_name, probe, cancer)) + 
    geom_violin() + geom_jitter(shape=16, position=position_jitter(0.2)) +
    geom_hline(yintercept=0.5, linetype="dashed", color = "grey") + stat_n_text(size=6) + 
    geom_hline(yintercept=0, linetype="dashed", color = "grey") + xlab("lncRNA expression group") +
    ylab("Beta Values") + stat_compare_means(label = "p.signif") +
    geom_boxplot(width=.1) + theme(text = element_text(size=15), axis.text = element_text(size=15)) 

    library(plyr)
    mu <- ddply(new, "median", summarise, grp.med=median(beta))
    head(mu)
    
    sp4 = ggplot(new, aes(x=beta, fill=median)) + ggtitle(paste(new$gene[1], gene_name, probe, cancer, "\nrisk=", new$risk[1])) + 
    geom_freqpoly(alpha=0.4, aes(color=median)) + geom_vline(data=mu, aes(xintercept=grp.med, color=median),
             linetype="dashed", size=2) +
    geom_vline(xintercept=0.5, linetype="dashed", color = "grey",size=0.5)+
    geom_vline(xintercept=0, linetype="dashed", color = "grey", size=0.5)+
    theme(text = element_text(size=15), axis.text = element_text(size=15))

    test = length(unique(new$mut_status))
    if(test >1){  

    #x-axis --> CNA status
    new$geneExp = as.numeric(new$geneExp)
    new$mut_status = factor(new$mut_status, levels = c("Unmethylated", "other", "Methylated"))
    #get correlation between cna_Status and gene expression
    stat_exp_cor = rcorr(new$mut_status, new$geneExp, type="spearman")$r[2]
    stat_exp_pval = rcorr(new$mut_status, new$geneExp, type="spearman")$P[2]

    text_add = paste("rho=", round(stat_exp_cor, digits=2), "pval=", round(stat_exp_pval, digits=8))
    ycord = as.numeric(summary(new$geneExp)[6]-0.75)

    #Create a custom color scale
    library(RColorBrewer)
    myColors <- brewer.pal(n = 3, name = "Set1")[c(2,1,3)]
    names(myColors) <- levels(new$mut_status)
    colScale <- scale_colour_manual(name = "mut_status",values = myColors)
    ks_test = kruskal.test(geneExp ~ mut_status, data = new)$p.value 

    sp6 = ggplot(new, aes(x=mut_status, y=geneExp, color=mut_status)) + ggtitle(paste(gene_name, cancer, probe)) + 
    geom_violin() + geom_jitter(shape=16, position=position_jitter(0.2)) +
    stat_n_text(size = 6) + 
    xlab("lncRNA Methylation status") +
    ylab("log1p(FPKM-UQ)") + stat_compare_means()+
    geom_boxplot(width=.1) + theme(text = element_text(size=15), axis.text = element_text(size=15))+
    annotate("text", x = 1.3, y =ycord , label = text_add)+
    scale_colour_manual(values=c("royalblue1", "gainsboro", "brown3"))

    pat_dat = new
    pat_dat$median = as.character(pat_dat$median)
    z = which(pat_dat$median == pat_dat$risk[1])
    pat_dat$median[z] = "RISK"
    pat_dat$median[-z] = "nonRISK"

    t = as.data.table(table(pat_dat$median, pat_dat$mut_status))
    t = as.data.table(filter(t, N >0))
    t = t[order(N)]
    
    new$group = ""

    if(HR < 1){
      new$group[new$mut_status == "Methylated"] = "RISK"
      new$group[new$mut_status %in% c("other", "Unmethylated")] = "nonRISK"
    }

    if(HR > 1){
      new$group[new$mut_status == "Unmethylated"] = "RISK"
      new$group[new$mut_status %in% c("other", "Methylated")] = "nonRISK"
    }

    print("pass3")

    if(!(length(unique(as.character(new$group)))) == 1) {

    #check if prognostic 
    cox_p = glance(coxph(Surv(OS.time, OS) ~ group, data = new))[8]
    newdat = new
    newdat$OS.time = newdat$OS.time/365
    newdat$group = factor(newdat$group, levels=c("RISK","nonRISK"))
    fit <- survfit(Surv(OS.time, OS) ~ group, data = newdat)
          s <- ggsurvplot(
          title = paste(cancer, name),
          fit, 
          xlab = "Time (Years)", 
          #surv.median.line = "hv",
          font.main = c(14, "bold", "black"),
          font.x = c(12, "plain", "black"),
          font.y = c(12, "plain", "black"),
          font.tickslab = c(11, "plain", "black"),
          font.legend = 10,
          risk.table.fontsize = 5, 
          legend.labs = c("Risk Methylation group", "Non-risk \nMethylation groups"),             # survfit object with calculated statistics.
          data = newdat,      # data used to fit survival curves. 
          risk.table = TRUE,       # show risk table.
          legend = "right", 
          pval = TRUE,             # show p-value of log-rank test.
          conf.int = FALSE,        # show confidence intervals for 
                            # point estimaes of survival curves.
          xlim = c(0,5),        # present narrower X axis, but not affect
                            # survival estimates.
          break.time.by = 1,     # break X axis in time intervals by 500.
          #palette = colorRampPalette(mypal)(14), 
          #palette = mypal[c(4,1)],
          palette = "npg", 
          #ggtheme = theme_bw(), # customize plot and risk table with a theme.
          risk.table.y.text.col = T, # colour risk table text annotations.
          risk.table.y.text = FALSE # show bars instead of names in text annotations
                            # in legend of risk table
          )
          print(s)

    print("pass 4")
    risk = pat_dat$risk[1]

    if(risk == "Low"){
      risk_met = filter(t, V1=="RISK", V2 == "Methylated")$N
      nonrisk_unmet = filter(t, V1=="nonRISK", V2 == "Unmethylated")$N
      risk_not_met = filter(t, V1=="RISK", V2 == "Unmethylated")$N
      nonrisk_met = filter(t, V1=="nonRISK", V2 == "Methylated")$N
      risk_other = filter(t, V1=="RISK", V2 == "other")$N
      nonrisk_other = filter(t, V1=="nonRISK", V2 == "other")$N
    }

    if(risk == "High"){
      risk_met = filter(t, V1=="RISK", V2 == "Methylated")$N
      nonrisk_unmet = filter(t, V1=="nonRISK", V2 == "Unmethylated")$N
      risk_not_met = filter(t, V1=="RISK", V2 == "Unmethylated")$N
      nonrisk_met = filter(t, V1=="nonRISK", V2 == "Methylated")$N
      risk_other = filter(t, V1=="RISK", V2 == "other")$N
      nonrisk_other = filter(t, V1=="nonRISK", V2 == "other")$N
    }

    if(length(risk_met) == 0){
      risk_met = 0
    }

    if(length(nonrisk_met) == 0){
      nonrisk_met = 0
    }

     if(length(risk_not_met) == 0){
      risk_not_met = 0
    }

     if(length(nonrisk_unmet) == 0){
      nonrisk_unmet = 0
    }

    if(length(risk_other) == 0){
      risk_other = 0
    }

    if(length(nonrisk_other) == 0){
      nonrisk_other = 0
    }

    other = length(unique(pat_dat$patient)) - risk_met - nonrisk_unmet - risk_not_met - nonrisk_met - risk_other - nonrisk_other

  #  results = c(cancer, unique(df$gene), length(unique(df$patient)), 
  #    wilcoxon_pval, chi_pval, fc, med_risk, length_risk_pats, rr, rp, ro, probe, bal, stat_exp_cor, stat_exp_pval, ks_test, risk_met, nonrisk_met, other,nonrisk_unmet, risk_not_met)
    results = c(cancer, unique(new$gene), length(unique(new$patient)), 
      wilcoxon_pval, chi_pval, fc, med_risk, length_risk_pats, rr, rp, ro, rop, probe, bal, stat_exp_cor, stat_exp_pval, ks_test, 
      risk_met, nonrisk_met, other, nonrisk_unmet, risk_not_met, risk_other, nonrisk_other, cox_p)
    
    names(results) = c("cancer", "gene", "num_patients", "wilcoxon_pval", "chi_pval", "fc", "risk_type", "num_risk_pats",  
      "risk_group_correlation", "nonrisk_group_correlation", "overall_correlation", "rop", "probe", "risk_pat_bal", "stat_exp_cor", "stat_exp_pval", "ks_test",
      "risk_met", "nonrisk_met", "other", "nonrisk_unmet", "risk_not_met", "risk_other", "nonrisk_other", "cox_p")

      print(sp6)
      print(lnc)
      print(probe)
    #print(lnc)
    }
    results_all_probes = rbind(results_all_probes, results)
  }
  }
  }

	results_all_probes = results_all_probes[-1,]
	results = results_all_probes
	}

	#only one probe present 
	if(probecount==1){

		#get wilcoxon p-value stored between low and high exp patients - get avg beta value for each group 
		wilcoxon_pval = wilcox.test(beta ~ median, data =df)$p.value  
		fc = median(df$beta[df$median == "Low"]) - median(df$beta[df$median == "High"])
    #methylation in low expression group relative to high expression group

    #label probe as methylated or not methylated for each patient 
		df$mut_status = ""
    		get_mut_stat = function(beta){
          meth = beta >= 0.75
          unmeth = beta < (0.25)
          other = ((beta >= 0.25) & (beta < 0.75))
          if(meth){
            return("Methylated")
             }        
           if(unmeth){
            return("Unmethylated")
           } 
           if(other){
            return("other")
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
    print(t)
    chi = chisq.test(t)
    chi_pval = tidy(chi)[2]

    #figure out if they are balanced or non-balanced 
    bal = cands[cands$data=="TCGA",]
    bal = as.numeric(bal$perc_risk[bal$combo == comb])

    length_risk_pats = length(which(df$median == med_risk))

    probe = df$probe[1]

    df$median <- factor(df$median, levels = c("High", "Low"))
    df$median  # notice the changed order of factor levels

    library("ggExtra")
    df$risk = med_risk

    library("ggExtra")
    
    gene_name = fantom$CAT_geneName[which(fantom$gene == lnc)]  

    sp1 = ggplot(df, aes(x=beta, y=geneExp, color=median)) + ggtitle(paste(df$gene[1], gene_name, probe, cancer, "\nrisk=", df$risk[1])) + 
    geom_point() + 
    geom_smooth(method=lm, se=FALSE, fullrange=TRUE) + stat_cor() +
    geom_vline(xintercept=0.5, linetype="dashed", color = "grey") + 
    geom_vline(xintercept=0, linetype="dashed", color = "grey") + xlab("Beta Values") +
    ylab("log1p FPKM") + theme(text = element_text(size=15), axis.text = element_text(size=15))


    sp2 = ggplot(df, aes(x=beta, y=geneExp)) + ggtitle(paste(df$gene[1], gene_name, probe, cancer, "\nrisk=", df$risk[1])) + 
    geom_point() + 
    geom_smooth(method=lm, se=FALSE, fullrange=TRUE) + stat_cor() +
    geom_vline(xintercept=0.5, linetype="dashed", color = "grey") + 
    geom_vline(xintercept=0, linetype="dashed", color = "grey") + xlab("Beta Values") +
    ylab("log1p FPKM") + theme(text = element_text(size=15), axis.text = element_text(size=15))


    sp3 = ggplot(df, aes(x=median, y=beta, color=median)) + ggtitle(paste(gene_name, probe, cancer)) + 
    geom_violin() + geom_jitter(shape=16, position=position_jitter(0.2)) +
    geom_hline(yintercept=0.5, linetype="dashed", color = "grey") + stat_n_text(size=6) + 
    geom_hline(yintercept=0, linetype="dashed", color = "grey") + xlab("lncRNA expression group") +
    ylab("Beta Values") + stat_compare_means(label = "p.signif") +
    geom_boxplot(width=.1) + theme(text = element_text(size=15), axis.text = element_text(size=15)) 

    library(plyr)
    mu <- ddply(df, "median", summarise, grp.med=median(beta))
    head(mu)
    
    sp4 = ggplot(df, aes(x=beta, fill=median)) + ggtitle(paste(df$gene[1], gene_name, probe, cancer, "\nrisk=", df$risk[1])) + 
    geom_freqpoly(alpha=0.4, aes(color=median)) + geom_vline(data=mu, aes(xintercept=grp.med, color=median),
             linetype="dashed", size=2) +
    geom_vline(xintercept=0.5, linetype="dashed", color = "grey",size=0.5)+
    geom_vline(xintercept=0, linetype="dashed", color = "grey", size=0.5) + 
    theme(text = element_text(size=15), axis.text = element_text(size=15))

    test = length(unique(df$mut_status))
    if(test >1){  

    #x-axis --> CNA status
    df$geneExp = as.numeric(df$geneExp)
    df$mut_status = factor(df$mut_status, levels = c("Unmethylated", "other", "Methylated"))
    #get correlation between cna_Status and gene expression
    stat_exp_cor = rcorr(df$mut_status, df$geneExp, type="spearman")$r[2]
    stat_exp_pval = rcorr(df$mut_status, df$geneExp, type="spearman")$P[2]

    text_add = paste("rho=", round(stat_exp_cor, digits=2), "pval=", round(stat_exp_pval, digits=8))
    ycord = as.numeric(summary(df$geneExp)[6]-0.75)

    #Create a custom color scale
    library(RColorBrewer)
    myColors <- brewer.pal(n = 3, name = "Set1")[c(2,1,3)]
    names(myColors) <- levels(df$mut_status)
    colScale <- scale_colour_manual(name = "mut_status",values = myColors)
    ks_test = kruskal.test(geneExp ~ mut_status, data = df)$p.value 

    sp6 = ggplot(df, aes(x=mut_status, y=geneExp, color=mut_status)) + ggtitle(paste(gene_name, cancer, probe)) + 
    geom_violin() + geom_jitter(shape=16, position=position_jitter(0.2)) +
    stat_n_text(size = 6) + 
    xlab("lncRNA Methylation status") +
    ylab("log1p(FPKM-UQ)") + stat_compare_means()+
    geom_boxplot(width=.1) + theme(text = element_text(size=15), axis.text = element_text(size=15))+
    annotate("text", x = 1.3, y =ycord , label = text_add)+
    scale_colour_manual(values=c("royalblue1", "gainsboro", "brown3"))

    #+
    #colScale +

    pat_dat = df
    pat_dat$median = as.character(pat_dat$median)
    z = which(pat_dat$median == pat_dat$risk[1])
    pat_dat$median[z] = "RISK"
    pat_dat$median[-z] = "nonRISK"

    t = as.data.table(table(pat_dat$median, pat_dat$mut_status))
    t = as.data.table(filter(t, N >0))
    t = t[order(N)]

    new = df
    new$group = ""

    if(HR < 1){
      new$group[new$mut_status == "Methylated"] = "RISK"
      new$group[new$mut_status %in% c("other", "Unmethylated")] = "nonRISK"
    }

    if(HR > 1){
      new$group[new$mut_status == "Unmethylated"] = "RISK"
      new$group[new$mut_status %in% c("other", "Methylated")] = "nonRISK"
    }


    if(!(length(unique(as.character(new$group)))) == 1){

    #check if prognostic 
    cox_p = glance(coxph(Surv(OS.time, OS) ~ group, data = new))[8]
    newdat = new
    newdat$OS.time = newdat$OS.time/365
    newdat$group = factor(newdat$group, levels=c("RISK","nonRISK"))
    fit <- survfit(Surv(OS.time, OS) ~ group, data = newdat)
          s <- ggsurvplot(
          title = paste(cancer, name),
          fit, 
          xlab = "Time (Years)", 
          #surv.median.line = "hv",
          font.main = c(14, "bold", "black"),
          font.x = c(12, "plain", "black"),
          font.y = c(12, "plain", "black"),
          font.tickslab = c(11, "plain", "black"),
          font.legend = 10,
          risk.table.fontsize = 5, 
          legend.labs = c("Risk Methylation group", "Non-risk \nMethylation groups"),             # survfit object with calculated statistics.
          data = newdat,      # data used to fit survival curves. 
          risk.table = TRUE,       # show risk table.
          legend = "right", 
          pval = TRUE,             # show p-value of log-rank test.
          conf.int = FALSE,        # show confidence intervals for 
                            # point estimaes of survival curves.
          xlim = c(0,5),        # present narrower X axis, but not affect
                            # survival estimates.
          break.time.by = 1,     # break X axis in time intervals by 500.
          #palette = colorRampPalette(mypal)(14), 
          #palette = mypal[c(4,1)],
          palette = "npg", 
          #ggtheme = theme_bw(), # customize plot and risk table with a theme.
          risk.table.y.text.col = T, # colour risk table text annotations.
          risk.table.y.text = FALSE # show bars instead of names in text annotations
                            # in legend of risk table
          )
          print(s)

    risk = pat_dat$risk[1]

    if(risk == "Low"){
      risk_met = filter(t, V1=="RISK", V2 == "Methylated")$N
      nonrisk_unmet = filter(t, V1=="nonRISK", V2 == "Unmethylated")$N
      risk_not_met = filter(t, V1=="RISK", V2 == "Unmethylated")$N
      nonrisk_met = filter(t, V1=="nonRISK", V2 == "Methylated")$N
      risk_other = filter(t, V1=="RISK", V2 == "other")$N
      nonrisk_other = filter(t, V1=="nonRISK", V2 == "other")$N
    }

    if(risk == "High"){
      risk_met = filter(t, V1=="RISK", V2 == "Methylated")$N
      nonrisk_unmet = filter(t, V1=="nonRISK", V2 == "Unmethylated")$N
      risk_not_met = filter(t, V1=="RISK", V2 == "Unmethylated")$N
      nonrisk_met = filter(t, V1=="nonRISK", V2 == "Methylated")$N
      risk_other = filter(t, V1=="RISK", V2 == "other")$N
      nonrisk_other = filter(t, V1=="nonRISK", V2 == "other")$N
    }

    if(length(risk_met) == 0){
      risk_met = 0
    }

    if(length(nonrisk_met) == 0){
      nonrisk_met = 0
    }

     if(length(risk_not_met) == 0){
      risk_not_met = 0
    }

     if(length(nonrisk_unmet) == 0){
      nonrisk_unmet = 0
    }

    if(length(risk_other) == 0){
      risk_other = 0
    }

    if(length(nonrisk_other) == 0){
      nonrisk_other = 0
    }

    other = length(unique(pat_dat$patient)) - risk_met - nonrisk_unmet - risk_not_met - nonrisk_met - risk_other - nonrisk_other

    results = c(cancer, unique(df$gene), length(unique(df$patient)), 
      wilcoxon_pval, chi_pval, fc, med_risk, length_risk_pats, rr, rp, ro, rop, probe, bal, stat_exp_cor, stat_exp_pval, ks_test, risk_met, nonrisk_met, 
      other,nonrisk_unmet, risk_not_met, risk_other, nonrisk_other, cox_p)
    
    names(results) = c("cancer", "gene", "num_patients", "wilcoxon_pval", "chi_pval", "fc", "risk_type", "num_risk_pats",  
      "risk_group_correlation", "nonrisk_group_correlation", "overall_correlation", "rop", "probe", "risk_pat_bal", "stat_exp_cor", "stat_exp_pval", "ks_test", 
      "risk_met", "nonrisk_met", "other", "nonrisk_unmet", "risk_not_met", "risk_other", "nonrisk_other", "cox_p")
    
    results = as.data.frame(results)

      #print(sp1)
      #print(sp2)
      #print(sp3)
      print(sp6)
      print(lnc)
      print("pass")
  	print(lnc)
  }
  }
	}
  }

	if(!((both & (length(unique(df$patient)) >=10)))){
		results = as.data.frame(matrix(ncol = 18)) ; colnames(results) = c("cancer", "gene", "num_patients", "wilcoxon_pval", "chi_pval", "fc", "risk_type", "num_risk_pats",  
      "risk_group_correlation", "nonrisk_group_correlation", "overall_correlation", "rop", "probe", "risk_pat_bal", "stat_exp_cor", "stat_exp_pval", "ks_test", "cox_p")
  }  
  return(results)
}
}
}
}

#pdf("candidate_lncRNAs_methylation_versus_Expression_only_NOFDR_candidates_Nov1.pdf")
#genes = as.list(unique(as.character(cands$combo[which(cands$combo %in% probes$combo)]))) #88/166 have methylation probes overlapping them 
#lnc_meth_cancer_data = llply(genes, get_data, .progress="text")
#dev.off()

#lnc_meth_cancer_data2 = Filter(Negate(is.null), lnc_meth_cancer_data)
#lnc_meth_cancer_data2 = ldply(lnc_meth_cancer_data2)
#lnc_meth_cancer_data2 = as.data.table(lnc_meth_cancer_data2)
#z= which(is.na(lnc_meth_cancer_data2$cancer))
#lnc_meth_cancer_data2 = lnc_meth_cancer_data2[-z,] #255 lncRNA-probe pairs evaluated 
#lnc_meth_cancer_data2$combo = paste(lnc_meth_cancer_data2$gene, lnc_meth_cancer_data2$cancer)
#saveRDS(lnc_meth_cancer_data2, file="new_results_methylation_Nov1.rds")
#73 unique lncRNA-cancer pairs evaluated 

lnc_meth_cancer_data2 = readRDS("new_results_methylation_Nov1.rds")
cands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
colnames(cands)[7] = "canc"

#---------PROCESS RESULTS-----------------------------------------------------------------------------------------------------

lnc_meth_cancer_data2$combo = paste(lnc_meth_cancer_data2$gene, lnc_meth_cancer_data2$cancer, sep="_") #55 combos evaluated 
lnc_meth_cancer_data2$wilcoxon_pval = as.numeric(as.character(lnc_meth_cancer_data2$wilcoxon_pval))
lnc_meth_cancer_data2 = lnc_meth_cancer_data2[order(wilcoxon_pval)]
lnc_meth_cancer_data2$wilcoxon_pval = as.numeric(lnc_meth_cancer_data2$wilcoxon_pval)
lnc_meth_cancer_data2$ks_test = as.numeric(as.character(lnc_meth_cancer_data2$ks_test))
lnc_meth_cancer_data2$stat_exp_cor = as.numeric(as.character(lnc_meth_cancer_data2$stat_exp_cor))
lnc_meth_cancer_data2$stat_exp_pval = as.numeric(as.character(lnc_meth_cancer_data2$stat_exp_pval))

#get fdr by cancer type - if only one test then there weren't multiple tests done
dats = split(lnc_meth_cancer_data2, by="cancer")
add_fdr = function(dat){
  if(dim(dat)[1]>5){
  dat$fdr = p.adjust(dat$ks_test, method="fdr")
  }
  if(dim(dat)[1] <=5){
    dat$fdr = dat$ks_test
  }
  return(dat)
}
lnc_meth_cancer_data2 = llply(dats, add_fdr)
lnc_meth_cancer_data2 = ldply(lnc_meth_cancer_data2)
sig_diff = filter(lnc_meth_cancer_data2, fdr <=0.05) #31 with sig wilcoxon pval 
sig_diff = as.data.table(sig_diff)
sig_diff = as.data.table(filter(sig_diff, abs(overall_correlation) >= 0.2))

#keep only ones with sig correlation
sig_diff = as.data.table(filter(sig_diff, rop < 0.05)) #17 sig positive correlation and fdr sig difference in dist

#plot just the sig ones 
sig_diff$combo = paste(sig_diff$gene, sig_diff$cancer, sep="_") #31/73 unique combos ~ 42% 

#z = which(sig_diff$other < 0)
#sig_diff = sig_diff[-z,] #23/73 unique combos ~ 32% 

#---------FIGURE SUMMARY FOR PAPER--------------------------------------------------------------------------------------------
canc_conv = readRDS("canc_conv.rds")
colnames(canc_conv)[2] = "canc"
cands = merge(cands, canc_conv, by="canc")
cands$combo = paste(cands$gene, cands$type, sep="_")
sig_diff = merge(sig_diff, cands, by=colnames(cands)[which(colnames(cands) %in% colnames(sig_diff))])
sig_diff = subset(sig_diff, data == "TCGA")
sig_diff$stat_exp_cor = unlist(sig_diff$stat_exp_cor)
sig_diff$stat_exp_cor = round(sig_diff$stat_exp_cor, digits=2)
sig_diff$ks_test = -log10(sig_diff$ks_test)

sig_diff$HR = as.numeric(sig_diff$HR)
sig_diff$stat[sig_diff$HR > 1] = "Unfavourable"
sig_diff$stat[sig_diff$HR < 1] = "Favourable"

#order 
sig_diff = as.data.table(sig_diff)
sig_diff = sig_diff[order(stat, abs(overall_correlation))]
sig_diff$CAT_geneName = factor(sig_diff$CAT_geneName, levels=unique(sig_diff$CAT_geneName))
sig_diff$canc = factor(sig_diff$canc, levels=unique(sig_diff$canc))
sig_diff$stat = factor(sig_diff$stat, levels=c("Unfavourable", "Favourable"))

#x-axis = cancer
#y-axis = lncRNA 

library(patchwork)

sig_diff$cor[sig_diff$stat_exp_cor < 0] = "Negative" 
sig_diff$cor[sig_diff$stat_exp_cor > 0] = "Positive" 
#z = which(is.na(sig_diff$cor))
#sig_diff = sig_diff[-z,]

#which have multiple probes per lncRNA? --> keep strongest absolute correlated one
t = as.data.table(table(sig_diff$combo))
t = as.data.table(filter(t, N > 1))
dups = unique(t$V1)

keep_dat = as.data.frame(matrix(ncol=ncol(sig_diff))) ; colnames(keep_dat) = colnames(sig_diff)
for(i in 1:length(dups)){
  dat = as.data.table(filter(sig_diff, combo == dups[i]))
  dat = dat[order(-(abs(stat_exp_cor)), -(abs(risk_group_correlation)))]
  keep = dat[1,]
  keep_dat = rbind(keep_dat, keep)
}
keep_dat = keep_dat[-1,]
z = which(!(sig_diff$combo %in% dups))
keep_dat = rbind(keep_dat, sig_diff[z,])
sig_diff = keep_dat
sig_diff = sig_diff[order(stat, -(abs(stat_exp_cor)))]

sig_diff$CAT_geneName = factor(sig_diff$CAT_geneName, levels=unique(sig_diff$CAT_geneName))
sig_diff$canc = factor(sig_diff$canc, levels=unique(sig_diff$canc))
sig_diff$stat = factor(sig_diff$stat, levels=c("Unfavourable", "Favourable"))
sig_diff = as.data.table(filter(sig_diff, stat_exp_pval < 0.05))
write.csv(sig_diff, file="29_lncRNAs_wmethylation_relationship.csv", quote=F, row.names=F)

sig_diff = as.data.table(filter(sig_diff, abs(stat_exp_cor) >= 0.2))

table(sig_diff$cor)
summary(abs(as.numeric(sig_diff$dist_toTSS))[sig_diff$cor == "Negative"])
summary(abs(as.numeric(sig_diff$dist_toTSS))[sig_diff$cor == "Positive"])

pdf("Methylation_figure_partB_sep27.pdf", width=10, height=8)

g = ggplot(sig_diff, aes(canc, CAT_geneName)) +
  geom_tile(aes(fill = ks_test)) + geom_point(aes(size=abs(stat_exp_cor), color=cor))+
    scale_fill_gradient(low = "lightgrey", high = "black", na.value = 'transparent') +
    scale_colour_manual(values = c("blue", "red")) + 
    xlab("Cancer") + ylab("lncRNA") + theme_bw() +
    theme(legend.position="bottom", legend.box = "horizontal", 
      legend.text=element_text(size=9), legend.title=element_text(size=6)) +
    guides(colour = guide_legend(order = 1),
         fill = guide_colourbar(order = 2),
         size = guide_legend(order = 3))

g = ggpar(g,
 font.tickslab = c(8,"plain", "black"),
 xtickslab.rt = 45)

#covariate favourable vs unfavoubrale
cov = ggplot(sig_diff, aes(CAT_geneName, 1)) + geom_tile(aes(fill = stat)) +
theme_void() + coord_flip() + theme(legend.position="none")

cov + g + plot_layout(ncol = 2, widths = c(1,15))

ggplot(sig_diff, aes(CAT_geneName, 1)) + geom_tile(aes(fill = stat)) +
theme_void() + coord_flip() 
 
dev.off()


##visualize probes 
probes_list = unique(sig_diff$probe)

get_plot_dist_probe_tss = function(pb){

probe_dat = filter(probes, cgid == pb)
probe_cord = probe_dat[,c("probe_chr", "cpg_start", "cpg_end", "cgid", "cpgstrand")]
lnc_cord = probe_dat[,c("lncchr", "lncstart", "lncend", "ensg", "lncstrand")]
colnames(lnc_cord) = c("probe_chr", "cpg_start", "cpg_end", "cgid", "cpgstrand")
probe_dat = rbind(probe_cord, lnc_cord)

probe_dat = melt(probe_dat)
probe_dat$yaxis = c(1.4,1.5)

g = ggplot(probe_dat, aes(x=value, y=yaxis, color=cgid)) +
  geom_line(aes(linetype=cgid))+
  scale_colour_manual(values = c("brown", "purple"))+
  geom_point(size=0) +
  geom_line(arrow = arrow(length=unit(0.20,"cm"), ends="first", type = "closed"), size = 1)+ 
  scale_y_continuous(breaks = c(0,1), labels = c(0,1))
print(g)
lprobetss = probe_dat$value[2] - probe_dat$value[1]
return(c(pb, probe_dat$cgid[2], lprobetss))

}

pdf("ucec_probe_figure.pdf", height=2)
probe_dat = llply(probes_list, get_plot_dist_probe_tss)
dev.off()

probe_dat = ldply(probe_dat)
saveRDS(probe_dat, "summary_probes_genes_oct10.rds")
colnames(probe_dat) = c("probe", "gene", "dist_toTSS")
sig_diff = merge(sig_diff, probe_dat, by= c("gene", "probe"))

#check cor vs dist
corsit = sig_diff[,c("probe", "gene", "cor", "dist_toTSS")]
corsit = corsit[order(cor)]

summary(abs(as.numeric(corsit$dist_toTSS))[corsit$cor == "Negative"])
summary(abs(as.numeric(corsit$dist_toTSS))[corsit$cor == "Positive"])
wilcox.test(abs(as.numeric(corsit$dist_toTSS))[corsit$cor == "Negative"], abs(as.numeric(corsit$dist_toTSS))[corsit$cor == "Positive"])

###############################
###stacked barplot#############
###############################

sig_diff$combo = paste(sig_diff$CAT_geneName, sig_diff$canc)
sig_diff$met_impact = (unlist(sig_diff$risk_met) + unlist(sig_diff$nonrisk_unmet)) /sig_diff$num_patients
sig_diff = sig_diff[order(-met_impact)]
sig_diff$combo = paste(sig_diff$CAT_geneName, sig_diff$type)

for_plot = sig_diff[,c("combo", "risk_met", "nonrisk_met", "other", "nonrisk_unmet", "risk_not_met", "risk_other", "nonrisk_other", "stat", "cor")]
cols = c("risk_met", "nonrisk_met", "other", "nonrisk_unmet", "risk_not_met", "risk_other", "nonrisk_other")

for(i in 1:nrow(for_plot)){
  #check what methylation risk group has 
  cor_t = for_plot$cor[i]
  risk = as.character(for_plot$stat[i])
  #check which one is large 
  check = for_plot$risk_met[i] >= for_plot$risk_not_met[i]
  
  if(check){
    #met group is larger check how many risk = met how many non risk = unmet
    #everyone else is other 
    risk_wmet = for_plot$risk_met[i]
    nonrisk_wmet = for_plot$nonrisk_unmet[i]
    other = for_plot$nonrisk_met[i] + for_plot$risk_not_met[i] + for_plot$risk_other[i] + for_plot$nonrisk_other[i]
  }

  if(!(check)){
    #unmet group is larger, how many risk = unmet, how many non risk = met
    #everyone else other
    risk_wmet = for_plot$risk_not_met[i]
    nonrisk_wmet = for_plot$nonrisk_met[i]
    other = for_plot$nonrisk_unmet[i] + for_plot$risk_met[i] +  for_plot$risk_other[i] + for_plot$nonrisk_other[i]
  }

  for_plot$riskwmet[i] = risk_wmet
  for_plot$nonrisk_wmet[i] = nonrisk_wmet
  for_plot$other_met[i] = other

}

#for_plot = melt(for_plot, measure.vars = cols)

#barplot = as.data.table(filter(for_plot, value >0))
#barplot = barplot[order(value)]
#barplot$variable = factor(barplot$variable, levels = c("risk_met", "risk_not_met", "nonrisk_met", "nonrisk_unmet", "risk_other", "nonrisk_other"))
#barplot$combo = factor(barplot$combo, levels=sig_diff$combo)

#barplot = merge(barplot, sig_diff, by=c("combo"))
barplot = for_plot
barplot$met_impact = barplot$riskwmet  + barplot$nonrisk_wmet 
barplot$total_pats = barplot$riskwmet  + barplot$nonrisk_wmet + barplot$other_met
barplot$met_impact = barplot$met_impact/barplot$total_pats

barplot = melt(barplot, measure.vars = c("riskwmet", "nonrisk_wmet", "other_met"))
barplot = as.data.table(filter(barplot, value >0))
barplot = barplot[,c("combo", "variable", "value", "met_impact", "total_pats")]
barplot$values = barplot$value/barplot$total_pats
barplot = barplot[order(-met_impact)]
barplot$combo = factor(barplot$combo, levels=unique(barplot$combo))
barplot$variable = factor(barplot$variable, levels = c("other_met", "nonrisk_wmet", "riskwmet"))

pdf("final_methylation_figure_partC.pdf", width=9, height=6)
g = ggplot(barplot, aes(combo, values))
g + geom_col(aes(fill = variable)) + scale_fill_manual(values=c("gainsboro",  "royalblue1", "brown3"), name="Methylation Type",
                       breaks=c("other_met", "nonrisk_wmet", "riskwmet"),
                       labels=c("Other", "NonRisk wMethylation", "Risk wMethylation"))+
theme_classic()+
theme(axis.text.x = element_text(size=9, angle=45, hjust=1),
          axis.text.y = element_text(size=12), legend.position="top") + xlab("lncRNA-cancer") + ylab("% of patients")
dev.off() 

saveRDS(barplot, file="methylation_data_23_candidates_final_figure_oct10.rds")

date = Sys.Date()

write.table(barplot, file=paste(date, "methylation_data_23_candidates_final_figure.csv", sep="_"), row.names=FALSE, na="", col.names = FALSE, sep=",")





