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

mypal = pal_npg("nrc", alpha = 0.7)(10)


#1. See if any candidates have CNAs
lncswcnas = fread("fantom_lncrnas_wTCGA_CNAs_23cancers.bed")
lncswcnas = as.data.frame(lncswcnas)

#2. cands 
cands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
colnames(cands)[7] = "canc"

colnames(lncswcnas)[4] = "gene"
colnames(lncswcnas)[11] = "canc"

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
genes = as.list(unique(as.character(cands$gene[which(cands$gene %in% lncswcnas$gene)]))) #146/166 have CNAs overlapping them 
colnames(lncswcnas) = c("lnc_chr", "lnc_start", "lnc_end", "gene", "name", "Chromosome" , "Start" , 
	"End", "Num_Probes" , "Segment_Mean", "canc", "rm", "patient")
lncswcnas$rm = NULL

rna = readRDS("5919_lncRNAs_tcga_all_cancers_March13_wclinical_data.rds")
unique(rna$type)

get_data = function(lnc){
	cancer = cands$canc[which(cands$gene == lnc)][1]
	dat = dplyr::filter(lncswcnas, canc == cancer, gene == lnc)
	dat$canc = NULL

	exp_data = subset(rna, type == cancer)
	#assign high or low to each patient in expression file
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

    df = df[,c(1:2, 3, 36:47)]
    df$V1 = NULL
    df$cna_status = ""
    get_cna_stat = function(seg_mean){
       dup = seg_mean >=0.2
       del = seg_mean < (-0.2)
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
    df$exp_cna_match = ""
    get_cna_exp = function(row){
      exp = row[[3]]
      cna = row[[15]]
      if(((exp == "High") & (cna == "DUP"))){
        return("HighExp_Duplicated")
      }
        if(((exp == "Low") & (cna == "DEL"))){
        return("LowExp_Deleted")
      }
        else{
          return("DontMatch")
        }
    }
    df$exp_cna_match = as.character(apply(df, 1, get_cna_exp))

    #overall correlation 
    #risk 
    #what is number of patients with higher risk?
    if(df$risk[1] == "high_expression"){
      med_risk = "High"
      cna_risk = "HighExp_Duplicated"
    }

    if(df$risk[1] == "low_expression"){
      med_risk = "Low"
      cna_risk = "LowExp_Deleted"
    }

    r = rcorr(df$Segment_Mean[df$median == med_risk], df$geneexp[df$median == med_risk], type="spearman")$r[2]
    rr = r #correlation in high risk group

    rp = rcorr(df$Segment_Mean[!(df$median == med_risk)], df$geneexp[!(df$median == med_risk)], type="spearman")$r[2]
    rp = rp #correlation in low risk group
    if(is.na(rp)){
      rp = 0
    }

    #get_wilcoxon_pval 
    wilcoxon_pval = wilcox.test(Segment_Mean ~ median, data =df)$p.value  

    length_risk_pats = length(which(df$med == med_risk))
    risk_pats = df[which(df$med == med_risk),]

    #what is the number of pateints wtih a cna that matches risk? (within risk grou)
    length_risk_pats_wcna = length(which(risk_pats$exp_cna_match == cna_risk))

    results = c(cancer, unique(df$gene), unique(df$name), length(unique(df$patient)), length(which(df$exp_cna_match == "HighExp_Duplicated")), 
      length(which(df$exp_cna_match == "LowExp_Deleted")), wilcoxon_pval, risk, length_risk_pats, length_risk_pats_wcna, rr, rp)
    names(results) = c("cancer", "gene", "name", "num_patients", "numHighCNAmatch", "numLowCNAmatch", "wilcoxon_pval", "risk_type", "num_risk_pats", "num_risk_pats_wmatchingCNA", "risk_group_correlation", "nonrisk_group_correlation")
    df$median <- factor(df$median, levels = c("Low", "High"))
    df$median  # notice the changed order of factor levels
  	df$risk = risk

    library("ggExtra")
	    sp = ggscatter(df, main = paste(df$gene[1], df$type[1]), 
	   	x = "Segment_Mean", y = "geneexp",
               color = "median", palette = mypal[c(4,1)],
               size = 3, alpha = 0.6, add = "reg.line",                         # Add regression line
          	   ggtheme = theme_light(), ylab = "log1p(FPKM) Expression", font.x = c(15, "plain", "black"), font.y = c(15, "plain", "black"),
          font.tickslab = c(15, "plain", "black"), xlab="Segment Mean SCNA") + stat_cor() 
    	print(sp)
	     xplot = ggboxplot(df, main= paste(df$name[1], cancer, "CNA vs Exp", "n=", length(unique(df$patient))),
		    x = "median", y = "Segment_Mean", legend.title = "Expression Tag", font.x = c(15, "plain", "black"), font.y = c(15, "plain", "black"),
          font.tickslab = c(15, "plain", "black"),
                  fill = "median", palette = mypal[c(4,1)], order=(c("Low", "High")), ggtheme = theme_light(), xlab="Expression", ylab="Segment Mean SCNA")+rotate()
    	xplot= xplot + stat_compare_means(label = "p.signif", label.x = 1.5)
    	print(xplot)
       yplot <- ggboxplot(df, x = "median", y = "geneexp", main = df$name[1], font.x = c(15, "plain", "black"), font.y = c(15, "plain", "black"),
          font.tickslab = c(15, "plain", "black"), fill = "median", palette = mypal[c(4,1)], order=(c("Low", "High")), ggtheme = theme_light(),
                    xlab="Expression", ylab="log1p(FPKM) Expression")
	     yplot= yplot + stat_compare_means(label = "p.signif", label.x = 1.5)
    
      print(yplot)
      print(lnc)
    return(results)
}
}
pdf("candidate_lncRNAs_CNA_versus_Expression_Aug13.pdf", height=5.5, width=8.3)
lnc_cna_cancer_data = llply(genes, get_data, .progress="text")
dev.off()

lnc_cna_cancer_data2 = as.data.frame(do.call("rbind", lnc_cna_cancer_data))
lnc_cna_cancer_data2 = as.data.table(lnc_cna_cancer_data2)


#---------PROCESS RESULTS-----------------------------------------------------------------------------------------------------

lnc_cna_cancer_data2$wilcoxon_pval = as.numeric(as.character(lnc_cna_cancer_data2$wilcoxon_pval))
lnc_cna_cancer_data2 = lnc_cna_cancer_data2[order(wilcoxon_pval)]
lnc_cna_cancer_data2$wilcoxon_pval = as.numeric(lnc_cna_cancer_data2$wilcoxon_pval)
lnc_cna_cancer_data2$fdr = p.adjust(lnc_cna_cancer_data2$wilcoxon_pval, method="fdr")

sig_diff = filter(lnc_cna_cancer_data2, wilcoxon_pval <=0.05)

lnc_cna_cancer_data2 = lnc_cna_cancer_data2[order(fdr, -numHighCNAmatch, -numLowCNAmatch)]

sig_diff = as.data.frame(sig_diff)
sig_diff$no_match = ""
get_nomatch = function(row){
  nomatch = as.numeric(as.character(row[[4]])) #all patients #
  high = as.numeric(as.character(row[[5]])) #of pats high match
  low = as.numeric(as.character(row[[6]]))
  nomatch = nomatch-high-low
}
sig_diff$no_match = apply(sig_diff, 1, get_nomatch)

match_h = sig_diff[,c(1:5,8, 9:13)]
colnames(match_h)[5] = "num_cna_exp_match"
match_h$type = "HighExpDup"
match_h$num_cna_exp_match = as.numeric(as.character(match_h$num_cna_exp_match))

match_l = sig_diff[,c(1:4, 6, 8, 9:13)]
colnames(match_l)[5] = "num_cna_exp_match"
match_l$type = "LowExpDel"
match_l$num_cna_exp_match = as.numeric(as.character(match_l$num_cna_exp_match))

match_no = sig_diff[,c(1:4,14,8, 9:13)]
colnames(match_no)[5] = "num_cna_exp_match"
match_no$type = "NoExpCNAMatch"

matched_sig = rbind(match_h, match_l, match_no)
matched_sig = as.data.table(matched_sig)
matched_sig = matched_sig[order(num_cna_exp_match, -type)]
matched_sig = as.data.frame(matched_sig)
order = as.character(unique(matched_sig$gene))
matched_sig$gene <- factor(matched_sig$gene, levels = order)
matched_sig$gene  # notice the changed order of factor levels

#for each lncRNA turn fractions into percentages 
matched_sig$num_cna_exp_match_patients = matched_sig$num_cna_exp_match
matched_sig$num_cna_exp_match = (as.numeric(matched_sig$num_cna_exp_match))/(as.numeric(as.character(matched_sig$num_patients)))
matched_sig$num_cna_exp_match = round(matched_sig$num_cna_exp_match, digits=3)

matched_sig = as.data.table(matched_sig)
matched_sig = matched_sig[order(risk_type, num_cna_exp_match, type)]

order = unique(matched_sig$name)
matched_sig$name <- factor(matched_sig$name, levels = order)

# Stacked bar plots, add labels inside bars
pdf("num_matching_CNA_exp_tags_23cancers_57lncs.pdf", width= 8, height=5)
p1 = ggbarplot(matched_sig, x = "name", y = "num_cna_exp_match",
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
z = which((matched_sig$risk_type == "high_expression") & (matched_sig$type == "HighExpDup"))
high = matched_sig[z,]

#low risk
z = which((matched_sig$risk_type == "low_expression") & (matched_sig$type == "LowExpDel"))
low = matched_sig[z,]

pats_summary = rbind(high, low)
order = as.character(unique(matched_sig$name))

saveRDS(matched_sig, file="lncRNA_CNA_summary_Aug13.rds")

pats_summary$per_patients_wrisk_cna = (as.numeric(as.character(pats_summary$num_risk_pats_wmatchingCNA)))/(as.numeric(as.character(pats_summary$num_risk_pats)))

pats_summary = pats_summary[order(match(name, order))] 

labels = as.character(pats_summary$cancer)

p2 = ggbarplot(pats_summary, x = "name", y = "per_patients_wrisk_cna", order=order, ylab="% of risk \n patients wCNA", xlab = "lncRNA", 
  color = "type", fill = "risk_type", legend="left", label=labels, lab.size = 2.5, lab.pos = c("out"), lab.vjust = 0.5, lab.hjust = 0.5, 
  palette = mypal) + theme_light() +coord_flip() + scale_colour_manual(values = c("salmon", "steelblue", "snow3")) +
  scale_fill_manual(values = c("mistyrose", "aliceblue")) + geom_hline(yintercept=0.2, linetype="dashed", color = "red") + 
  ggtitle("lncRNA CNAs in Risk Groups")
  #label = TRUE, lab.col = "white", lab.pos = "in", lab.size=1.8)
p2= ggpar(p2,
 font.tickslab = c(8,"plain", "black"),
 xtickslab.rt = 45, legend = "none") 

library(patchwork)
pdf("CNAs_lncRNAs_summary_30cands_may30.pdf", width=12, height=9)
p1 + p2 + plot_layout(ncol = 2, widths = c(2, 1))
dev.off()

#divide up the risk correlation and no risk correlation groups
cors_risk = pats_summary[,c(1:9, 11:14)]
cors_nonrisk = pats_summary[,c(1:8, 10, 11:14)]

colnames(cors_risk)[9] = "cor_spearman"
cors_risk$cor_group = "risk"
colnames(cors_nonrisk)[9] = "cor_spearman"
cors_nonrisk$cor_group = "nonrisk"

pats_summary = rbind(cors_risk, cors_nonrisk)

#plot the overall correlation between lncRNA copy number and expression in cohort 
pats_summary$cor_spearman = as.numeric(as.character(pats_summary$cor_spearman))

p3 = ggbarplot(pats_summary, x = "name", y = "cor_spearman", order=order, ylab="Spearman correlation \n CNA wExpression", xlab = "lncRNA", 
  add = "segments", palette = mypal, color = "cor_group", position = position_dodge(0.6)) + theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +  theme_light() +
  coord_flip() + 
  geom_hline(yintercept=c(0), linetype="dashed", color = "black") + scale_colour_manual(values = c("steelblue", "salmon")) +
  ggtitle("lncRNA CNAs in Risk Groups") + geom_hline(yintercept=c(0.15), linetype="dashed", color = "red") 
  #label = TRUE, lab.col = "white", lab.pos = "in", lab.size=1.8)
p3= ggpar(p3,
 font.tickslab = c(8,"plain", "black"),
 xtickslab.rt = 45, legend = "right") 

pdf("CNAs_lncRNAs_summary_30cands_may30_wcor.pdf", width=14, height=8)
p1 + p2 + p3 + plot_layout(ncol = 3, widths = c(2, 1, 1))
dev.off()



##############----------------------------------------------------------------------------------------------------################
##############----------------------------------RISK--------------------------------------------------------------################
##############----------------------------------------------------------------------------------------------------################


###save a new dataframe with each patient labelled as high risk/low risk, carry CNA or deletion or no CNA at all 
###also get average correlation values 


#1. for each patient figure out if they are low risk or high risk group 
head(sig_diff)

#keep gene-cancer combinations, don't really care right now if gene has CNA
#for a different cancer where it's not a candidate 
genes = as.list(unique(as.character(cands$gene[which(cands$gene %in% lncswcnas$gene)]))) #146/166 have CNAs overlapping them 

rna = readRDS("5919_lncRNAs_tcga_all_cancers_March13_wclinical_data.rds")
unique(rna$type)

get_data = function(lnc){
  cancer = cands$canc[which(cands$gene == lnc)][1]
  dat = dplyr::filter(lncswcnas, canc == cancer, gene == lnc)
  dat$canc = NULL

  exp_data = subset(rna, type == cancer)
  #assign high or low to each patient in expression file
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
    if(HR >1){
      risk = "high_expression"
    }
    if(HR <1){
      risk = "low_expression"
    }
    df$risk = risk
    #is copy number aberation associated with expression? 
    df$geneexp = log1p(df$geneexp)
    df = df[,c(1:2, 3, 36:47)]
    df$V1 = NULL
    df$cna_status = ""
    get_cna_stat = function(seg_mean){
       dup = seg_mean >=0.1
       del = seg_mean < (-0.1)
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
    df$exp_cna_match = ""
    get_cna_exp = function(row){
      exp = row[[3]]
      cna = row[[15]]
      if(((exp == "High") & (cna == "DUP"))){
        return("HighExp_Duplicated")
      }
        if(((exp == "Low") & (cna == "DEL"))){
        return("LowExp_Deleted")
      }
        else{
          return("DontMatch")
        }
    }
    df$exp_cna_match = as.character(apply(df, 1, get_cna_exp))
    df$cancer = cancer
    df$length_cna_segment = df$End-df$Start
    #convert to megabase 1mil
    df$length_cna_segment = df$length_cna_segment/1000000
    return(df)
}
}

lnc_cna_cancer_data_files_save = llply(genes, get_data, .progress="text")
lnc_cna_cancer_data_files_save = as.data.frame(do.call("rbind", lnc_cna_cancer_data_files_save))
lnc_cna_cancer_data_files_save = as.data.table(lnc_cna_cancer_data_files_save)

#plot summarize length of copy number aberations 
#by cancer type 

pdf("summary_all_146lncRNAs_copynumber_May30.pdf", width=9)
g = ggviolin(lnc_cna_cancer_data_files_save, x = "cancer", y  = "length_cna_segment", xlab="length CNA (Kb)", color="risk", draw_quantiles = 0.5, palette=mypal) + theme_light()
 ggpar(g,
 font.tickslab = c(8,"plain", "black"),
 xtickslab.rt = 45, legend = "right") 

dev.off()

#only subset of 44 lncRNAs with associations between cna and expression 
smaller_subset = subset(lnc_cna_cancer_data_files_save, gene %in% pats_summary$gene)

pdf("summary_all_44_lncRNAs_copynumber_May30.pdf", width=9)
g = ggviolin(smaller_subset, x = "cancer", y  = "length_cna_segment", xlab="length CNA (Kb)", color="risk", draw_quantiles = 0.5, palette=mypal) + theme_light() + ggtitle("Only lncRNAs with sig Wilcoxon")
 ggpar(g,
 font.tickslab = c(8,"plain", "black"),
 xtickslab.rt = 45, legend = "right") 
dev.off()



##############----------------------------------------------------------------------------------------------------################
##############----------------------------------FINAL FILE SAVE---------------------------------------------------################
##############----------------------------------------------------------------------------------------------------################

saveRDS(lnc_cna_cancer_data_files_save, file="146_lncRNAs_wCNA_data_patient_status_May30th.rds")














