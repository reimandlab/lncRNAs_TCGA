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

tcga_results1 = readRDS(file="TCGA_FMRE_diff_exp_PCGS_survival_results_July3.rds")
cancers_conv = rna[,which(colnames(rna) %in% c("type", "Cancer"))]
cancers_conv = cancers_conv[!duplicated(cancers_conv), ]
colnames(cancers_conv)[2] = "cancer"
tcga_results1 = merge(tcga_results1, cancers_conv, by="cancer")
tcga_results1 = tcga_results1[order(fdr_pval)]
tcga_results1$HR = as.numeric(tcga_results1$HR)

#Remove with unrealistic larger HRs
z = which(tcga_results1$upper95 == "Inf")
tcga_results1 = tcga_results1[-z,]
z = which(tcga_results1$HR >=3)

tcga_results1$HR = log2(tcga_results1$HR)

#replace fdr values with **
tcga_results1$fdr_text = ""
tcga_results1$fdr_text[tcga_results1$fdr_pval <= 0.05] = "*"
tcga_results1$fdr_text[tcga_results1$fdr_pval <= 0.01] = "**"
tcga_results1$fdr_text[tcga_results1$fdr_pval <= 0.001] = "***"
tcga_results1$fdr_text[tcga_results1$fdr_pval <= 0.0001] = "****"

#keep only those with ***
#z = which(tcga_results1$fdr_text == "")
#tcga_results1 = tcga_results1[-z,]

#order by increasing fdr
tcga_results1 = tcga_results1[order(fdr_pval)]
cancers = unique(tcga_results1$type)
tcga_results1$type = factor(tcga_results1$type, levels = cancers)
genes = unique(tcga_results1$name)
tcga_results1$name = factor(tcga_results1$name, levels = genes)

#--------------------------------------------------------------------
#PLOT - results 
#--------------------------------------------------------------------

#x = gene
#y = cancer

#color of genes names by their coefficeint
tcga_results1 = merge(tcga_results1, cands, by="gene")
tcga_results1 = tcga_results1[order(fdr_pval)]

#remove cancer types with no *
empty = as.data.table(table(tcga_results1$type, tcga_results1$fdr_text))
empty = as.data.table(filter(empty, N >0))
empty = as.data.table(filter(empty, V2 != ""))

tcga_results1 = subset(tcga_results1, type %in% empty$V1)

a <- ifelse(unique(tcga_results1$mut_coef_lm) > 0, "red", "blue")

pdf("10_diff_exp_cands_survival_results_TCGA_july3.pdf")
g = ggplot(tcga_results1, aes(name, type)) +
  geom_tile(aes(fill = HR)) +
  geom_text(aes(label = fdr_text), size=5) +
    scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0, na.value = 'transparent') +
    xlab("Gene") + ylab("Cancer") + theme_bw()
ggpar(g,
 font.tickslab = c(7,"plain", "black"),
 xtickslab.rt = 45, legend.title="log2(HR)") + 
 theme(axis.text.x = element_text(colour = a)) 

dev.off()

#Split into two - one heatmap for genes with -ve coef, one for +ve coef

neg = subset(tcga_results1, mut_coef_lm <0)
pos = subset(tcga_results1, mut_coef_lm > 0)

a <- ifelse(unique(neg$mut_coef_lm) > 0, "red", "blue")
n = ggplot(neg, aes(name, type)) +
  geom_tile(aes(fill = HR)) +
  geom_text(aes(label = fdr_text), size=5) +
    scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0, na.value = 'transparent') +
    xlab("FMRE target") + ylab("Cancer type, TCGA") + theme_bw()
n = ggpar(n,
 font.tickslab = c(7,"plain", "black"),
 xtickslab.rt = 45, legend="none") + ggtitle("Down-regulation of \nFMRE target gene") + 
 theme(axis.text.x = element_text(colour = "black"), plot.title=element_text(size=9), axis.title=element_text(size=9),
  axis.text=element_text(size=8)) +
   theme(plot.title = element_text(hjust = 0.5))


a <- ifelse(unique(pos$mut_coef_lm) > 0, "red", "blue")
p = ggplot(pos, aes(name, type)) +
  geom_tile(aes(fill = HR)) +
  geom_text(aes(label = fdr_text), size=5) +
    scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0, na.value = 'transparent') +
    xlab("FMRE target") + ylab("Cancer") + theme_bw()
p = ggpar(p,
 font.tickslab = c(7,"plain", "black"),
 xtickslab.rt = 45, legend.title="log2(HR)") + ggtitle("Up-regulation of \nFMRE target gene")+
theme(axis.text.x = element_text(colour = "black"), plot.title=element_text(size=9), axis.title=element_text(size=9),
  axis.text=element_text(size=8)) +
  theme(plot.title = element_text(hjust = 0.5)) +
   theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())


pdf("10_diff_exp_cands_survival_results_TCGA_july3_split_two.pdf")
n + p + plot_layout(ncol = 2, widths = c(1, 3)) 
dev.off()


