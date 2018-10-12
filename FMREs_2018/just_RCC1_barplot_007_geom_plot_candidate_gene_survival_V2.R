library(survAUC)
#source("source_code_Cox_MonteCarlo_CV_Mar13.R")
require(caTools)
#check if this person is in my analysis: TCGA-61-2095
library(glmnet)
library(survcomp)
library(caret)
library(stringr)

source("universal_LASSO_survival_script.R")

library(ggpubr)
library(ggrepel)
library(viridis)
library(patchwork)
library(caret)  
library(Rtsne)

#------FEATURES-----------------------------------------------------

#UCSC gene info
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,8]))
ucsc <- ucsc[-z,]

cands = fread("FMRE_difexp_target_genes.txt")

#--------This script ------------------------------------------------

#Candidate gene survival analysis, Vol2 
#These are differentially expressed genes in mutated samples in PCAWG
#here we summarize the survival results for each of these candidates 
#within each cancer type 
#Author = Karina Isaev
#Date = July 3rd, 2018

#Can you please perform median dichotomized survival analysis for 
#each of these genes across all cancer types in TCGA?

#Two expected outputs:

###----1. a geom_tile plot of gene X cancer type, HR to color the tile, asterisks to show significance of 
#HR from CoxPH (. 0.1 * 0.05 ** 0.01 *** 0.001 **** 0.0001). Genes ranked by significance from table. 
#We should somehow color or annotate gene names in plot to show those with positive and negative coefficients 
#shown in the table (pos coefficient - gene up regulated with mutations). 

###----2. all gene X cancer type KM plots, one file per gene, p-value ranked. 


#------PROCESS CANDIDATES--------------------------------------------

#1. get genes ids 
clean_gene = function(ge){
  newge = unlist(strsplit(ge, "::"))[2]
  #get ensembl id
  z = which(ucsc$hg19.ensemblToGeneName.value == newge)
  newge = ucsc$hg19.ensGene.name2[z]
  return(newge)
}

cands$gene = unlist(llply(cands$gene, clean_gene))

#--------------------------------------------------------------------
#DATA - results 
#--------------------------------------------------------------------

tcga_results1 = readRDS(file="TCGA_FMRE_diff_exp_PCGS_survival_results_Oct10.rds")
cancers_conv = rna[,which(colnames(rna) %in% c("type", "Cancer"))]
cancers_conv = cancers_conv[!duplicated(cancers_conv), ]
colnames(cancers_conv)[2] = "cancer"
tcga_results1 = merge(tcga_results1, cancers_conv, by="cancer")
tcga_results1 = tcga_results1[order(fdr_pval)]
tcga_results1$HR = as.numeric(tcga_results1$HR)

t = as.data.table(table(rna$type))
t
t = as.data.table(t, N > 50)
t
t = t[order(N)]
t
t = as.data.table(filter(t, N >50))
tcga_results1 = as.data.table(filter(tcga_results1, type %in% t$V1))

#Remove with unrealistic larger HRs
z = which(tcga_results1$upper95 == "Inf")
tcga_results1 = tcga_results1[-z,]
z = which(tcga_results1$HR >=3)
tcga_results1 = tcga_results1[order(pval)]
write.csv(tcga_results1, file="TCGA_FMRE_diff_exp_PCGs_survival_results_Oct10.csv", quote=F, row.names=F)


tcga_results1$HR = log2(tcga_results1$HR)

#replace fdr values with **
tcga_results1$fdr_text = ""
tcga_results1$fdr_text[tcga_results1$pval <= 0.05] = "*"
tcga_results1$fdr_text[tcga_results1$pval <= 0.01] = "**"
tcga_results1$fdr_text[tcga_results1$pval <= 0.001] = "***"
tcga_results1$fdr_text[tcga_results1$pval <= 0.0001] = "****"

#keep only those with ***
#z = which(tcga_results1$fdr_text == "")
#tcga_results1 = tcga_results1[-z,]

#order by increasing fdr
tcga_results1 = tcga_results1[order(pval)]
cancers = unique(tcga_results1$type)
tcga_results1$type = factor(tcga_results1$type, levels = cancers)
genes = unique(tcga_results1$name)
tcga_results1$name = factor(tcga_results1$name, levels = genes)

#--------------------------------------------------------------------
#PLOT - results 
#--------------------------------------------------------------------

#x = Cancer type
#y = Hazard Ratio for RCC1 
#color = p-value 

rcc1 = as.data.table(filter(tcga_results1, name=="RCC1"))
rcc1 = rcc1[order(pval)]
rcc1$type = factor(rcc1$type, levels=rcc1$type)

pdf("rcc1_hrs_cancers_all.pdf")
g <- ggplot(rcc1, aes(type, HR))
# Number of cars in each class:
g + geom_col(aes(fill = -log10(pval))) +
theme_bw()+
theme(axis.text.x = element_text(size=10, angle=45, hjust=1),
          axis.text.y = element_text(size=10), legend.position="top") + xlab("Cancer Type") + ylab("log2(Hazard Ratio)")+
scale_fill_gradient(low = "grey", high = "black")
dev.off()

#ordered by Hazard ratios
rcc1 = as.data.table(filter(tcga_results1, name=="RCC1"))
rcc1 = rcc1[order(HR)]
rcc1$type = factor(rcc1$type, levels=rcc1$type)

pdf("rcc1_hrs_cancers_all_ordered_by_Hrs.pdf")
g <- ggplot(rcc1, aes(type, HR))
# Number of cars in each class:
g + geom_col(aes(fill = -log10(pval))) +
theme_bw()+
theme(axis.text.x = element_text(size=10, angle=45, hjust=1),
          axis.text.y = element_text(size=10), legend.position="top") + xlab("Cancer Type") + ylab("log2(Hazard Ratio)")+
scale_fill_gradient(low = "grey", high = "black")
dev.off()

#ordered normal pvalues 
rcc1 = as.data.table(filter(tcga_results1, name=="RCC1"))
rcc1 = rcc1[order(HR)]
rcc1$type = factor(rcc1$type, levels=rcc1$type)

pdf("rcc1_hrs_cancers_all_ordered_by_hrs_normal_pvalues.pdf")
g <- ggplot(rcc1, aes(type, HR))
# Number of cars in each class:
g + geom_col(aes(fill = (pval))) +
theme_bw()+
theme(axis.text.x = element_text(size=10, angle=45, hjust=1),
          axis.text.y = element_text(size=10), legend.position="top") + xlab("Cancer Type") + ylab("log2(Hazard Ratio)")+
scale_fill_gradient(low = "black", high = "grey")
dev.off()


#just sig ones
rcc1 = as.data.table(filter(tcga_results1, name=="RCC1", pval < 0.05))
rcc1 = rcc1[order(pval)]
rcc1$type = factor(rcc1$type, levels=rcc1$type)

pdf("rcc1_hrs_cancers_just_sig_ones.pdf")
g <- ggplot(rcc1, aes(type, HR))
# Number of cars in each class:
g + geom_col(aes(fill = -log10(pval))) +
theme_bw()+
theme(axis.text.x = element_text(size=10, angle=45, hjust=1),
          axis.text.y = element_text(size=10), legend.position="top") + xlab("Cancer Type") + ylab("log2(Hazard Ratio)")+
scale_fill_gradient(low = "grey", high = "black")
dev.off()

##############
###this one###
##############

#just sig ones ordered by HRs 
rcc1 = as.data.table(filter(tcga_results1, name=="RCC1", pval < 0.05))
rcc1 = rcc1[order(-HR)]
rcc1$type = factor(rcc1$type, levels=rcc1$type)
rcc1$factor = "slateblue4"

pdf("rcc1_hrs_cancers_just_sig_ones_ordered_by_HR.pdf")
g <- ggplot(rcc1, aes(type, HR, color="factor"))
# Number of cars in each class:
g + geom_col(aes(fill = -log10(pval))) +
theme_bw()+
scale_color_manual(values = c("slateblue4")) +
theme(axis.text.x = element_text(size=10, angle=45, hjust=1),
          axis.text.y = element_text(size=10), legend.position="top") + xlab("Cancer Type") + ylab("log2(Hazard Ratio)")+
scale_fill_gradient(low = "grey", high = "black")+
guides(color = FALSE)
dev.off()




#just sig ones ordered by HRs and normal pvalues
rcc1 = as.data.table(filter(tcga_results1, name=="RCC1", pval < 0.05))
rcc1 = rcc1[order(HR)]
rcc1$type = factor(rcc1$type, levels=rcc1$type)

pdf("rcc1_hrs_cancers_just_sig_ones_ordered_by_HR_normal_pval.pdf")
g <- ggplot(rcc1, aes(type, HR))
# Number of cars in each class:
g + geom_col(aes(fill = (pval))) +
theme_bw()+
theme(axis.text.x = element_text(size=10, angle=45, hjust=1),
          axis.text.y = element_text(size=10), legend.position="top") + xlab("Cancer Type") + ylab("log2(Hazard Ratio)")+
scale_fill_gradient(low = "black", high = "grey")
dev.off()





