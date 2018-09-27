#Libraries#------------------------------------------------
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

mypal = pal_npg("nrc", alpha = 0.7)(10)

#1. See if any candidates have CNAs
lncswcnas = fread("fantom_lncrnas_wTCGA_CNAs_23cancers.bed")
lncswcnas = as.data.frame(lncswcnas)

#2. cands 
cands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
colnames(cands)[7] = "canc"

colnames(lncswcnas)[6] = "gene"
colnames(lncswcnas)[15] = "canc"

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

lncswcnas$canc[lncswcnas$canc=="ovary"] = "OV"
lncswcnas$canc[lncswcnas$canc=="liver"] = "LIHC"
lncswcnas$canc[lncswcnas$canc=="kidney"] = "KIRC"
lncswcnas$canc[lncswcnas$canc=="brca"] = "BRCA"
lncswcnas$canc[lncswcnas$canc=="luad"] = "LUAD"
lncswcnas$canc[lncswcnas$canc=="pancreas"] = "PAAD"

#turn cancers to capitals 
lncswcnas = lncswcnas %>% 
      mutate(canc = toupper(canc))

#keep gene-cancer combinations, don't really care right now if gene has CNA
#for a different cancer where it's not a candidate 
lncswcnas$combo = ""
lncswcnas$combo = paste(lncswcnas$gene, lncswcnas$canc, sep="_")
cands$combo = paste(cands$gene, cands$canc, sep="_")

colnames(lncswcnas) = c("lnc_chr", "lnc_start", "lnc_end", "width", "type", "gene" , "type_lnc" , "name", "name2", "Chromosome", "Start", 
	"End", "Num_Probes" , "Segment_Mean", "canc", "rm" ,"patient", "combo")

lncswcnas$rm = NULL

#HOW MANY THOUGH HAVE SEGMENTS AMP/DEL GREATER THAN 30 MEGABASE? 
lncswcnas$length = lncswcnas$End-lncswcnas$Start
lncswcnas$length = lncswcnas$length/1000000
lncswcnas = subset(lncswcnas, length <20)

genes = as.list(unique(as.character(cands$combo[which(cands$combo %in% lncswcnas$combo)]))) #152/173 have CNAs overlapping them 
#with segments that are shorter than 5 MB

rna = readRDS("5919_lncRNAs_tcga_all_cancers_March13_wclinical_data.rds")
unique(rna$type)

# Create the function.
getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}

get_data = function(lnc){
	 
  comb = lnc

  cancer = cands$canc[which(cands$combo == lnc)][1]
	
  dat = dplyr::filter(lncswcnas, canc == cancer, gene == cands$gene[which(cands$combo == lnc)][1])
	dat$canc = NULL

  if((dim(dat)[1]) > 30){

	exp_data = subset(rna, type == cancer)
	#assign high or low to each patient in expression file

  lnc = cands$gene[which(cands$combo == lnc)][1]

	z <- which(colnames(exp_data) %in% lnc)
  	if(!(length(z)==0)){
  	df = as.data.frame(exp_data)
  	df <- df[,c(1, z,(ncol(exp_data)-32):ncol(exp_data))]  

	df$median <- ""
 	median2 <- quantile(as.numeric(df[,2]), 0.5)
  	
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

  	#get summary of SCNA in lncRNA for each patient 
  	#take mean segment mean for all cnas in patient 
  	#covering that lncRNA
  	df = merge(df, dat, by=c("patient"))
  	colnames(df)[2] = "geneexp"
  	
    #is high expression or low expression associated with worse prognosis? 
    df$median <- factor(df$median, levels = c("Low", "High"))
    df$median  # notice the changed order of factor levels

    test1 = table(df$median)[1]
    test2 = table(df$median)[2]

    if((test1 > 10) & (test2 > 10)){

    HR = summary(coxph(Surv(OS.time, OS) ~ median, data = df))$coefficients[2]
    
    fit <- survfit(Surv(OS.time, OS) ~ median, data = df)

    if(HR >1){
      risk = "high_expression"
    }
    if(HR <1){
      risk = "low_expression"
    }

    df$risk = risk

    #is copy number aberation associated with expression? 
  	df$geneexp = log1p(df$geneexp)

    #df = df[,c(1:2, 3, 36:47)]
    df$V1 = NULL
    df$cna_status = ""
    get_cna_stat = function(seg_mean){
       dup = seg_mean >=0.3
       del = seg_mean < (-0.3)
       if(dup){
        return("DUP")
       }        
       if(del){
        return("DEL")
       } 
       if((!(dup)) & (!(del))){
        return("noCNA")
       }
    }
    df$cna_status = as.character(llply(df$Segment_Mean, get_cna_stat))

    df$median = as.character(df$median)
    
    #overall correlation 
    #risk 
    #what is number of patients with higher risk?
    if(df$risk[1] == "high_expression"){
      med_risk = "High"
    }

    if(df$risk[1] == "low_expression"){
      med_risk = "Low"
    }

    c_cor = (table(df$median)[1] >=5) & (table(df$median)[2] >=5) 

    if(c_cor){
    r = rcorr(df$Segment_Mean[df$median == med_risk], df$geneexp[df$median == med_risk], type="spearman")$r[2]
    rr = r #correlation in high risk group

    rp = rcorr(df$Segment_Mean[!(df$median == med_risk)], df$geneexp[!(df$median == med_risk)], type="spearman")$r[2]
    rp = rp #correlation in low risk group
    if(is.na(rp)){
      rp = 0
    }
    }

    if(!(c_cor)){
      rr = "cant"
      rp="cant"
    }

    #overall correlation
    ro = rcorr(df$Segment_Mean, df$geneexp, type="spearman")$r[2]

    #get_wilcoxon_pval 
    wilcoxon_pval = wilcox.test(Segment_Mean ~ median, data =df)$p.value  

    #get_chisq_pval
    t = table(df$median, df$cna_status)
    chi = chisq.test(t)
    chi_pval = tidy(chi)[2]

    length_risk_pats = length(which(df$med == med_risk))

    #figure out if they are balanced or non-balanced 
    bal = cands[cands$data=="TCGA",]
    bal = as.numeric(bal$perc_risk[bal$combo == comb])
    df$length = (df$End - df$Start)/1000000
    avg_length = mean(df$length)

    results = c(cancer, unique(df$gene), unique(df$name2), length(unique(df$patient)), 
      wilcoxon_pval, chi_pval, 
      risk, length_risk_pats, rr, rp, ro, bal, avg_length)
    
    names(results) = c("cancer", "gene", "name", "num_patients", "wilcoxon_pval", "chi_pval", 
      "risk_type", "num_risk_pats", "risk_group_correlation", "nonrisk_group_correlation", "overall_correlation", "balance_risk_pats", 
      "avg_length")
    
    df$median <- factor(df$median, levels = c("High", "Low"))
    df$median  # notice the changed order of factor levels
  	df$risk = risk
    
    print(lnc)
    #plot length of segments 
    print(hist(df$length))
    pat_dat = df[,c("patient", "median", "cna_status", "name", "name2", "combo", "length", "Segment_Mean", "Num_Probes")]
    pat_dat$cancer = cancer
    return(pat_dat)
}
}
}
}

#each df for each combo - saved patients with CNAs and CNA status 
lnc_cna_cancer_data = llply(genes, get_data, .progress="text")
lnc_cna_cancer_data2 = as.data.frame(do.call("rbind", lnc_cna_cancer_data))
lnc_cna_cancer_data2 = as.data.table(lnc_cna_cancer_data2)
saveRDS(lnc_cna_cancer_data2, file="all_dfs_CNAs_status_exp_patients_sept26.rds")

r = readRDS("all_dfs_CNAs_status_exp_patients_sept26.rds")
r$combo = paste(r$name2, r$cancer)

print("done sir")

get_dens_plot = function(lnc){

  dat = as.data.table(filter(r, combo == lnc))
  
  library(plyr)
  mu <- ddply(dat, "median", summarise, grp.med=median(Segment_Mean))
  head(mu)

  sp5 = ggplot(dat, aes(x=Segment_Mean, fill=median)) + ggtitle(dat$name2[1]) + 
    geom_density(alpha = 0.5) + geom_vline(xintercept = 0.3, linetype="dotted", 
                color = "grey") + geom_vline(xintercept = -0.3, linetype="dotted", color = "grey") +
    theme(text = element_text(size=6), axis.text = element_text(size=6), legend.position="none")
  return(sp5)
}

get_boxplot_plot = function(lnc){

  dat = as.data.table(filter(r, combo == lnc))
  
  sp5 = ggplot(dat, aes(x=median, y=Segment_Mean, color=median)) + 
    geom_violin() + geom_jitter(shape=16, position=position_jitter(0.2)) +
    geom_hline(yintercept=0.3, linetype="dashed", color = "grey") + stat_n_text(size = 3) + 
    geom_hline(yintercept=-0.3, linetype="dashed", color = "grey") + xlab("lncRNA expression group") +
    ylab("Segment Mean SCNA") + stat_compare_means(label = "p.signif")+
    geom_boxplot(width=.1) + theme(text = element_text(size=6), axis.text = element_text(size=6), legend.position="none")

  return(sp5)
}

sig_diff = readRDS("sig_diff_CNAs_sept26.rds")
sig_diff = as.data.table(sig_diff)
sig_diff = sig_diff[order(fdr)]
sig_diff$combo = paste(sig_diff$name, sig_diff$cancer)

lncs = unique(as.character(sig_diff$combo)) #65 lncRNAs have both gene expression and copy number data available 
#with enough pateints for analysis (>30 patients) & segment mean of max 10mb

list_cna_plots = llply(lncs, get_dens_plot)

library(gridExtra)
n <- length(list_cna_plots)
nCol <- floor(sqrt(n))

arranged = align_plots(plotlist = list_cna_plots, align="hv", axis="tblr")

pdf("sig_cnas_density_plots.pdf", height=10)
do.call("grid.arrange", c(arranged, nrow=7, ncol=2))
dev.off()
dev.off()

list_cna_boxplots = llply(lncs, get_boxplot_plot)

library(gridExtra)
n <- length(list_cna_boxplots)
nCol <- floor(sqrt(n))

arranged = align_plots(plotlist = list_cna_boxplots, align="hv", axis="tblr")

pdf("sig_cnas_boxplots_plots.pdf", height=13, width=10)
do.call("grid.arrange", c(arranged, nrow=7, ncol=2))
dev.off()
dev.off()












