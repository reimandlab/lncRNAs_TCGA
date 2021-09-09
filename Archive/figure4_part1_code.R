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

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = subset(allCands, data == "TCGA") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 
allCands$combo = unique(paste(allCands$gene, allCands$cancer, sep="_"))

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

cands = as.data.table(filter(cands, data == "TCGA"))

colnames(lncswcnas) = c("lnc_chr", "lnc_start", "lnc_end", "width", "type", "gene" , "type_lnc" , "name", "name2", "Chromosome", "Start", 
	"End", "Num_Probes" , "Segment_Mean", "canc", "rm" ,"patient", "combo")

lncswcnas$rm = NULL

#HOW MANY THOUGH HAVE SEGMENTS AMP/DEL GREATER THAN 30 MEGABASE? 
lncswcnas$length = lncswcnas$End-lncswcnas$Start
lncswcnas$length = lncswcnas$length/1000000
lncswcnas = subset(lncswcnas, length <20) #10 KB

genes = as.list(unique(as.character(cands$combo[which(cands$combo %in% lncswcnas$combo)]))) #148/168 have CNAs overlapping them 
#with segments that are shorter than 20 MB

rna = readRDS("5919_lncRNAs_tcga_all_cancers_March13_wclinical_data.rds")
table(rna$type)

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

  if((dim(dat)[1]) > 20){

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

    if((test1 >=10) & (test2 >= 10)){

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
        return("AMP")
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

    sp5 = ggplot(df, aes(x=Segment_Mean, fill=median)) + ggtitle(paste(df$name2[1], cancer)) + 
    geom_density(alpha = 0.5) + geom_vline(data=mu, aes(xintercept=grp.med, color=median),
             linetype="dashed") + theme(text = element_text(size=15), axis.text = element_text(size=15))

    #x-axis --> CNA status
    df$geneexp = as.numeric(df$geneexp)
    df$cna_status = factor(df$cna_status, levels = c("DEL", "noCNA", "AMP"))
    #get correlation between cna_Status and gene expression
    stat_exp_cor = rcorr(df$cna_status, df$geneexp, type="spearman")$r[2]
    stat_exp_pval = rcorr(df$cna_status, df$geneexp, type="spearman")$P[2]

    text_add = paste("rho=", round(stat_exp_cor, digits=2), "pval=", round(stat_exp_pval, digits=8))
    ycord = as.numeric(summary(df$geneexp)[6]-0.75)

    #Create a custom color scale
    library(RColorBrewer)
    myColors <- brewer.pal(n = 3, name = "Set1")[c(2,3,1)]
    names(myColors) <- levels(df$cna_status)
    colScale <- scale_colour_manual(name = "cna_status",values = myColors)
    
    test = length(unique(df$cna_status))
    if(test >1){  

    ks_test = kruskal.test(geneexp ~ cna_status, data = df)$p.value 

    sp6 = ggplot(df, aes(x=cna_status, y=geneexp, color=cna_status)) + ggtitle(paste(df$name2[1], cancer)) + 
    geom_violin() + geom_jitter(shape=16, position=position_jitter(0.2)) +
    stat_n_text(size = 6) + 
    xlab("lncRNA CNA status") +
    ylab("log1p(FPKM-UQ)") + stat_compare_means()+
    geom_boxplot(width=.1) + theme(text = element_text(size=15), axis.text = element_text(size=15))+
    scale_colour_manual(values=c("royalblue1", "gainsboro", "brown3")) +
    annotate("text", x = 1.3, y =ycord , label = text_add)

      #print(sp1)
      #print(sp2)
      #print(sp3)
      #print(sp6)
      print(lnc)
    #plot length of segments 
    print(hist(df$length))
    pat_dat = df[,c("patient", "median", "cna_status", "name", "name2", "combo", "risk")]
    pat_dat$cancer = cancer
    pat_dat$median = as.character(pat_dat$median)
    pat_dat$median[pat_dat$median == "Low"] = "low_expression"
    pat_dat$median[pat_dat$median == "High"] = "high_expression"
    z = which(pat_dat$median == pat_dat$risk[1])
    pat_dat$median[z] = "RISK"
    pat_dat$median[-z] = "nonRISK"

    t = as.data.table(table(pat_dat$median, pat_dat$cna_status))
    t = as.data.table(filter(t, N >0))
    t = t[order(N)]

    if(risk == "low_expression"){
      risk_cna = filter(t, V1=="RISK", V2 == "DEL")$N
      nonrisk_cna = filter(t, V1=="nonRISK", V2 == "AMP")$N
    }

    if(risk == "high_expression"){
      risk_cna = filter(t, V1=="RISK", V2 == "AMP")$N
      nonrisk_cna = filter(t, V1=="nonRISK", V2 == "DEL")$N
    }

    if(length(risk_cna) == 0){
      risk_cna = 0
    }

    if(length(nonrisk_cna) == 0){
      nonrisk_cna = 0
    }

    other = length(unique(df$patient)) - risk_cna - nonrisk_cna

    results = c(cancer, unique(df$gene), unique(df$name2), length(unique(df$patient)), 
      wilcoxon_pval, chi_pval, 
      risk, length_risk_pats, rr, rp, ro, bal, avg_length, stat_exp_cor, stat_exp_pval, ks_test, risk_cna, nonrisk_cna, other)
    
    names(results) = c("cancer", "gene", "name", "num_patients", "wilcoxon_pval", "chi_pval", 
      "risk_type", "num_risk_pats", "risk_group_correlation", "nonrisk_group_correlation", "overall_correlation", "balance_risk_pats", 
      "avg_length", "stat_exp_cor", "stat_exp_pval", "ks_test", "risk_cna", "nonrisk_cna", "other")

    return(results)
}
}
}
}
}
#pdf("candidate_lncRNAs_CNA_versus_Expression_Nov1_new_plots.pdf")

lnc_cna_cancer_data = llply(genes, get_data, .progress="text")
#dev.off()

#lnc_cna_cancer_data2 = as.data.frame(do.call("rbind", lnc_cna_cancer_data))
#lnc_cna_cancer_data2 = as.data.table(lnc_cna_cancer_data2)
#saveRDS(lnc_cna_cancer_data2, file="new_results_CNAs_Sept27.rds")

lnc_cna_cancer_data2 = readRDS("new_results_CNAs_Sept27.rds") #evaluated 71 lncRNA-cancer combos

#---------PROCESS RESULTS-----------------------------------------------------------------------------------------------------

lnc_cna_cancer_data2$wilcoxon_pval = as.numeric(as.character(lnc_cna_cancer_data2$wilcoxon_pval))
lnc_cna_cancer_data2 = lnc_cna_cancer_data2[order(wilcoxon_pval)]
lnc_cna_cancer_data2$wilcoxon_pval = as.numeric(lnc_cna_cancer_data2$wilcoxon_pval)
lnc_cna_cancer_data2$cancer = as.character(lnc_cna_cancer_data2$cancer)
lnc_cna_cancer_data2$stat_exp_cor = as.numeric(as.character(lnc_cna_cancer_data2$stat_exp_cor))
lnc_cna_cancer_data2$stat_exp_pval = as.numeric(as.character(lnc_cna_cancer_data2$stat_exp_pval))
lnc_cna_cancer_data2$ks_test = as.numeric(as.character(lnc_cna_cancer_data2$ks_test))

#get fdr by cancer type - if only one test then there weren't multiple tests done
dats = split(lnc_cna_cancer_data2, by="cancer")
add_fdr = function(dat){
  if(dim(dat)[1]>10){
  dat$fdr = p.adjust(dat$ks_test, method="fdr")
  }
  if(dim(dat)[1] <=10){
    dat$fdr = dat$ks_test
  }
  return(dat)
}
lnc_cna_cancer_data2 = llply(dats, add_fdr)
lnc_cna_cancer_data2 = ldply(lnc_cna_cancer_data2)
lnc_cna_cancer_data2 = as.data.table(lnc_cna_cancer_data2)

lnc_cna_cancer_data2 = as.data.table(filter(lnc_cna_cancer_data2, abs(stat_exp_cor) >= 0.2))

#check which have sig overall correlation and sig kruskal wallis 
sig_diff = as.data.table(filter(lnc_cna_cancer_data2, fdr <=0.05))
#positive correlation
sig_diff = as.data.table(filter(sig_diff, stat_exp_cor > 0))

#keep only ones with sig correlation
sig_diff = as.data.table(filter(sig_diff, stat_exp_pval < 0.05)) #17 sig positive correlation and fdr sig difference in dist

saveRDS(sig_diff, file="sig_diff_CNAs_sept27.rds")

#plot just the sig ones 
sig_diff$combo = paste(sig_diff$gene, sig_diff$cancer, sep="_")
sig_diff = sig_diff[order(fdr)]
genes = as.list(unique(sig_diff$combo))

#pdf("FDR_sig_candidate_lncRNAs_CNA_versus_Expression_Sept_JUST_SIG_ONES_sep28.pdf")
#lnc_cna_cancer_data = llply(genes, get_data, .progress="text")
#dev.off()

#---------FIGURE SUMMARY FOR PAPER--------------------------------------------------------------------------------------------
sig_diff$gene = as.character(sig_diff$gene)
sig_diff = merge(sig_diff, cands, by=colnames(cands)[which(colnames(cands) %in% colnames(sig_diff))])
sig_diff = subset(sig_diff, data == "TCGA")
sig_diff$stat_exp_cor = unlist(sig_diff$stat_exp_cor)
sig_diff$stat_exp_cor = round(sig_diff$stat_exp_cor, digits=2)
sig_diff$ks_test = -log10(sig_diff$ks_test)

sig_diff$HR = as.numeric(sig_diff$HR)
sig_diff$stat[sig_diff$HR > 1] = "Unfavourable"
sig_diff$stat[sig_diff$HR < 1] = "Favourable"

#remove negative correlations 
sig_diff = subset(sig_diff, stat_exp_cor >0)

#order 
sig_diff = as.data.table(sig_diff)
sig_diff = sig_diff[order(stat, -(abs(stat_exp_cor)))]
sig_diff$CAT_geneName = factor(sig_diff$CAT_geneName, levels=unique(sig_diff$CAT_geneName))
sig_diff$canc = factor(sig_diff$canc, levels=unique(sig_diff$canc))
sig_diff$stat = factor(sig_diff$stat, levels=c("Unfavourable", "Favourable"))

sig_diff$cor[sig_diff$stat_exp_cor < 0] = "Negative" 
sig_diff$cor[sig_diff$stat_exp_cor > 0] = "Positive" 
sig_diff = as.data.table(filter(sig_diff, stat_exp_pval < 0.05)) #17 sig positive correlation and fdr sig difference in dist

#x-axis = cancer
#y-axis = lncRNA 

library(patchwork)

pdf("CNA_figure_partA_sept27.pdf", width=10, height=7)

g = ggplot(sig_diff, aes(canc, CAT_geneName)) +
  geom_tile(aes(fill = ks_test)) + geom_point(aes(size=abs(stat_exp_cor), color=cor))+
    scale_fill_gradient(low = "lightgrey", high = "black", na.value = 'transparent') +
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


###############################
###stacked barplot#############
###############################

sig_diff$total_cnas = unlist(sig_diff$risk_cna) + unlist(sig_diff$nonrisk_cna) + unlist(sig_diff$other)
sig_diff = sig_diff[order(-total_cnas)]
sig_diff$combo = paste(sig_diff$CAT_geneName, sig_diff$canc)
sig_diff$cna_impact = (unlist(sig_diff$risk_cna) + unlist(sig_diff$nonrisk_cna))/sig_diff$total_cnas
sig_diff = sig_diff[order(-cna_impact)]


barplot_1 = sig_diff[,c("gene", "name", "cancer", "canc", "CAT_geneName", "ks_test", "stat_exp_cor", "cor", "num_risk", 
  "perc_risk", "risk_cna")]
barplot_1$type = "risk_cna"
colnames(barplot_1)[ncol(barplot_1)-1] = "cna_counts"
barplot_1$cna_counts = unlist(barplot_1$cna_counts)/sig_diff$total_cnas

barplot_2 = sig_diff[,c("gene", "name", "cancer", "canc", "CAT_geneName", "ks_test", "stat_exp_cor", "cor", "num_risk", 
  "perc_risk", "nonrisk_cna")]
barplot_2$type = "nonrisk_cna"
colnames(barplot_2)[ncol(barplot_2)-1] = "cna_counts"
barplot_2$cna_counts = unlist(barplot_2$cna_counts)/sig_diff$total_cnas

barplot_3 = sig_diff[,c("gene", "name", "cancer", "canc", "CAT_geneName", "ks_test", "stat_exp_cor", "cor", "num_risk", 
  "perc_risk", "other")]
barplot_3$type = "other"
colnames(barplot_3)[ncol(barplot_3)-1] = "cna_counts"
barplot_3$cna_counts = unlist(barplot_3$cna_counts)/sig_diff$total_cnas

barplot = rbind(barplot_1, barplot_2, barplot_3)

barplot = barplot[order(type)]
barplot$combo = paste(barplot$CAT_geneName, barplot$canc)
barplot$type = factor(barplot$type, levels = c("other", "nonrisk_cna", "risk_cna"))
barplot$combo = factor(barplot$combo, levels=sig_diff$combo)

pdf("final_cna_figure_partA.pdf", width=9, height=6)
g <- ggplot(barplot, aes(combo, cna_counts))
# Number of cars in each class:
g + geom_col(aes(fill = type)) + scale_fill_manual(values=c("gainsboro", "royalblue1", "brown3"), name="CNA Type",
                       breaks=c("other", "nonrisk_cna", "risk_cna"),
                       labels=c("Other", "Non-Risk wCNA", "Risk wCNA"))+
theme_bw()+
theme(axis.text.x = element_text(size=10, angle=45, hjust=1),
          axis.text.y = element_text(size=14), legend.position="top") + xlab("lncRNA-cancer") + ylab("% of patients")
dev.off()

saveRDS(sig_diff, file="cna_data_16_candidates_final_figure_oct10.rds")

data1<-data.frame(lapply(sig_diff, as.character), stringsAsFactors=FALSE)
date = Sys.Date()

write.table(data1, file=paste(date, "lncRNAs_sig_associated_with_CNAs.csv", sep="_"), row.names=FALSE, na="", col.names = FALSE, sep=",")




































