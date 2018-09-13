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
lncswcnas = subset(lncswcnas, length <30)

genes = as.list(unique(as.character(cands$combo[which(cands$combo %in% lncswcnas$combo)]))) #153/173 have CNAs overlapping them 
#with segments that are shorter than 5 MB

rna = readRDS("5919_lncRNAs_tcga_all_cancers_March13_wclinical_data.rds")
unique(rna$type)

# Create the function.
getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}

pat_dat_all = as.data.frame(matrix(ncol = 7)) ; colnames(pat_dat_all) = c("patient", "median", "cna_status", "name", "name2", "combo", "cancer")

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

    library("ggExtra")

    sp1 = ggplot(df, aes(x=Segment_Mean, y=geneexp, color=median)) + ggtitle(paste(df$name2[1], cancer)) + 
    geom_point() + 
    geom_smooth(method=lm, se=FALSE, fullrange=TRUE) + stat_cor() +
    geom_vline(xintercept=0.2, linetype="dashed", color = "grey") + 
    geom_vline(xintercept=-0.2, linetype="dashed", color = "grey") + xlab("Segment Mean SCNA") +
    ylab("log1p FPKM") + theme(text = element_text(size=15), axis.text = element_text(size=15))


    sp2 = ggplot(df, aes(x=Segment_Mean, y=geneexp)) + ggtitle(paste(df$name2[1], cancer)) + 
    geom_point() + 
    geom_smooth(method=lm, se=FALSE, fullrange=TRUE) + stat_cor() +
    geom_vline(xintercept=0.2, linetype="dashed", color = "grey") + 
    geom_vline(xintercept=-0.2, linetype="dashed", color = "grey") + xlab("Segment Mean SCNA") +
    ylab("log1p FPKM") + theme(text = element_text(size=15), axis.text = element_text(size=15))


    sp3 = ggplot(df, aes(x=median, y=Segment_Mean, color=median)) + ggtitle(paste(df$name2[1], cancer)) + 
    geom_violin() + geom_jitter(shape=16, position=position_jitter(0.2)) +
    geom_hline(yintercept=0.2, linetype="dashed", color = "grey") + stat_n_text(size = 6) + 
    geom_hline(yintercept=-0.2, linetype="dashed", color = "grey") + xlab("lncRNA expression group") +
    ylab("Segment Mean SCNA") + stat_compare_means(label = "p.signif")+
    geom_boxplot(width=.1) + theme(text = element_text(size=15), axis.text = element_text(size=15))

    library(plyr)
    mu <- ddply(df, "median", summarise, grp.med=median(Segment_Mean))
    head(mu)
    
    sp4 = ggplot(df, aes(x=Segment_Mean, fill=median)) + ggtitle(paste(df$name2[1], cancer)) + 
    geom_freqpoly(alpha=0.4, aes(color=median)) + geom_vline(data=mu, aes(xintercept=grp.med, color=median),
             linetype="dashed") + theme(text = element_text(size=15), axis.text = element_text(size=15))

      print(sp1)
      print(sp2)
      print(sp3)
      print(sp4)
      print(lnc)
    #plot length of segments 
    print(hist(df$length))
    pat_dat = df[,c("patient", "median", "cna_status", "name", "name2", "combo")]
    pat_dat$cancer = cancer
    pat_dat_all = rbind(pat_dat_all, pat_dat)
    return(results)
}
}
}
pdf("candidate_lncRNAs_CNA_versus_Expression_Sept7.pdf")
lnc_cna_cancer_data = llply(genes, get_data, .progress="text")
dev.off()

lnc_cna_cancer_data2 = as.data.frame(do.call("rbind", lnc_cna_cancer_data))
lnc_cna_cancer_data2 = as.data.table(lnc_cna_cancer_data2)
saveRDS(lnc_cna_cancer_data2, file="new_results_CNAs_Sept7.rds")

#---------PROCESS RESULTS-----------------------------------------------------------------------------------------------------

lnc_cna_cancer_data2$wilcoxon_pval = as.numeric(as.character(lnc_cna_cancer_data2$wilcoxon_pval))
lnc_cna_cancer_data2 = lnc_cna_cancer_data2[order(wilcoxon_pval)]
lnc_cna_cancer_data2$wilcoxon_pval = as.numeric(lnc_cna_cancer_data2$wilcoxon_pval)
lnc_cna_cancer_data2$cancer = as.character(lnc_cna_cancer_data2$cancer)

#get fdr by cancer type - if only one test then there weren't multiple tests done
dats = split(lnc_cna_cancer_data2, by="cancer")
add_fdr = function(dat){
  if(dim(dat)[1]>5){
  dat$fdr = p.adjust(dat$wilcoxon_pval, method="fdr")
  }
  if(dim(dat)[1] <=5){
    dat$fdr = dat$wilcoxon_pval
  }
  return(dat)
}
lnc_cna_cancer_data2 = llply(dats, add_fdr)
lnc_cna_cancer_data2 = ldply(lnc_cna_cancer_data2)
sig_diff = filter(lnc_cna_cancer_data2, fdr <=0.05)

#plot just the sig ones 
sig_diff$combo = paste(sig_diff$gene, sig_diff$cancer, sep="_")
genes = as.list(unique(sig_diff$combo))
pdf("FDR_sig_candidate_lncRNAs_CNA_versus_Expression_Sept_JUST_SIG_ONES.pdf")
lnc_cna_cancer_data = llply(genes, get_data, .progress="text")
dev.off()


#---------FIGURE SUMMARY FOR PAPER--------------------------------------------------------------------------------------------

sig_diff = merge(sig_diff, cands, by=colnames(cands)[which(colnames(cands) %in% colnames(sig_diff))])
sig_diff = subset(sig_diff, data == "TCGA")
sig_diff$overall_correlation = unlist(sig_diff$overall_correlation)
sig_diff$overall_correlation = round(sig_diff$overall_correlation, digits=2)
sig_diff$wilcoxon_pval = -log10(sig_diff$wilcoxon_pval)

sig_diff$HR = as.numeric(sig_diff$HR)
sig_diff$stat[sig_diff$HR > 1] = "Unfavourable"
sig_diff$stat[sig_diff$HR < 1] = "Favourable"

#remove negative correlations 
sig_diff = subset(sig_diff, overall_correlation >0)

#order 
sig_diff = as.data.table(sig_diff)
sig_diff = sig_diff[order(stat, -(abs(overall_correlation)))]
sig_diff$CAT_geneName = factor(sig_diff$CAT_geneName, levels=unique(sig_diff$CAT_geneName))
sig_diff$canc = factor(sig_diff$canc, levels=unique(sig_diff$canc))
sig_diff$stat = factor(sig_diff$stat, levels=c("Unfavourable", "Favourable"))


sig_diff$cor[sig_diff$overall_correlation < 0] = "Negative" 
sig_diff$cor[sig_diff$overall_correlation > 0] = "Positive" 


#x-axis = cancer
#y-axis = lncRNA 

library(patchwork)

pdf("CNA_figure_partA.pdf", width=10, height=7)

g = ggplot(sig_diff, aes(canc, CAT_geneName)) +
  geom_tile(aes(fill = wilcoxon_pval)) + geom_point(aes(size=abs(overall_correlation), color=cor))+
    scale_fill_gradient(low = "white", high = "black", na.value = 'transparent') +
    scale_colour_manual(values = c("red", "red")) + 
    xlab("Cancer") + ylab("lncRNA") + theme_bw() +
    theme(legend.position="bottom", legend.box = "horizontal", 
      legend.text=element_text(size=9), legend.title=element_text(size=6))
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








