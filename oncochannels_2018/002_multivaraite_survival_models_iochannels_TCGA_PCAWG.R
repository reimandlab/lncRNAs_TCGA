
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#author: Karina Isaev, karin.isaev@gmail.com 
#date updated: Sept 21, 2018
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


#-------------------------------------------------------------------
#this script uses data from TCGA (gene expression and clinical) to evaluate 
#the prognositc value of ion channels across different cancer types 
#working directory (source files, data files)
#/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ

#------Load libraries and scripts-----------------------------------

library(survAUC)
#source("source_code_Cox_MonteCarlo_CV_Mar13.R")
require(caTools)
#check if this person is in my analysis: TCGA-61-2095
library(glmnet)
library(survcomp)
library(caret)
library(stringr)

#this script below prepares the RNA and clinical files for analysis 
source("universal_LASSO_survival_script.R")

library(ggpubr)
library(ggrepel)
library(viridis)
library(patchwork)
library(caret)  
library(Rtsne)

#------DATA---------------------------------------------------------
#UCSC gene info
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
#z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
#ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,8]))
ucsc <- ucsc[-z,]
colnames(ucsc)[8] = "HGNC.symbol"

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

#remove cancer types with less than 50 patients 
pats_num = as.data.table(table(rna$Cancer))
pats_num = filter(pats_num, N <50)
canc_rm = pats_num$V1
rna = rna[-which(rna$Cancer %in% canc_rm),]

#Combined into one dataframe because need to get ranks 
com = colnames(pcg)[which(colnames(pcg) %in% colnames(rna))]
all <- merge(rna, pcg, by = com)

#------FEATURES-----------------------------------------------------

#cands -- ion channels 
cands = read.csv("ION_CHANNELS_targets_and_families.csv")
cands = merge(cands, ucsc, by = "HGNC.symbol")

#outlier gene exp 
load("expr_discr_cc.rsav")

#run fisher's test to see if any ion channels within cancer types 
#are signifcantly assoicated with each other



gbm = expr_discr_cc[[1]]
lgg = expr_discr_cc[[27]]


#add surv info 
genes_test = c("CATSPER1", "SCN9A", "AQP9", "KCNN4")


#add IDH mutation status 


#-------------------------------------------------------------------
#------SUMMARY------------------------------------------------------
#-------------------------------------------------------------------

tcga_results1 = readRDS("TCGA_ION_CHANNEL_results_Sept21.rds")
tcga_results1$data = "TCGA"
tcga_results1$combo = paste(tcga_results1$gene, tcga_results1$cancer, sep="_")

#check which models violate the PH assumption
#to those models add age * survival time interaction 
which(tcga_results1$global_test_ph <= 0.05) 
tcga_results1$num_risk = as.numeric(tcga_results1$num_risk)

#remove those with risk group less than 15 patients 
z = which(tcga_results1$num_risk <15)
tcga_results1 = tcga_results1[-z,]

#plot sig KM plots
sig = as.data.table(filter(tcga_results1, pval <= 0.05))
source("check_lnc_exp_cancers.R")

#convert back to cancer codes
canc_conv = rna[,c("type", "Cancer")]
canc_conv = unique(canc_conv)
colnames(canc_conv)[2] = "cancer"
sig = merge(sig, canc_conv, by="cancer")
sig = sig[order(fdr_pval)]


#230 sig after fdr 
#plot summary 

#1. # of oncochannels (FDR < 0.05) per cancer type 
sig$name = unlist(llply(sig$gene, get_name_pcg))
#269 unique genes across 25/28 cancer types 
#169 unique genes across 13/28 cancers FDR < 0.05
#474 Favourable, 507 Unfavourable 

fdr_sig = as.data.table(filter(sig, fdr_pval < 0.05))
#5 genes appeared in > 3 cancer types :
#P2RX5
#KCNK5
#GJC1
#CLCN5
#P2RX6 across ACC, KIRC, LGG, MESO and UVM 

#cancer types with most FDR sig genes: LGG, KIRC, ACC, KIRP and MESO 

#2. # of oncochannels (p-val < 0.05) per cancert type
sig$HR = as.numeric(sig$HR)
sig$status[sig$HR <1] = "Favourable"
sig$status[sig$HR >1] = "Unfavourable"

sum = as.data.table(table(sig$type, sig$status))
sum = sum[order(V2, N)]
sum=as.data.table(filter(sum, N >1))

canc_order = as.data.table(table(sig$type))
canc_order = canc_order[order(-N)]
canc_order=as.data.table(filter(canc_order, N >1))

sum$V1 = factor(sum$V1, levels = canc_order$V1)

pdf("onco_channels_summary_prog_pval0.05.pdf", width=9)
g = ggbarplot(sum, "V1", "N",
  fill = "V2", color = "V2", palette = "Paired",
  label = TRUE, lab.col = "white", lab.pos = "in")
ggpar(g,font.tickslab = c(10,"plain", "black"),
 xtickslab.rt = 45)
dev.off()

#3. oncochannels summary (which ones most prognostic) 

fdr_sig = as.data.table(filter(sig,fdr_pval<0.05))
z = which(is.na(fdr_sig$type))
if(!(length(z)==0)){
  fdr_sig = fdr_sig[-z,]
}

sum = as.data.table(table(fdr_sig$type, fdr_sig$status))
sum = sum[order(V2, N)]
sum=as.data.table(filter(sum, N >1))

canc_order = as.data.table(table(fdr_sig$type))
canc_order = canc_order[order(-N)]
canc_order=as.data.table(filter(canc_order, N >1))

sum$V1 = factor(sum$V1, levels = canc_order$V1)

pdf("onco_channels_summary_prog_fdr0.05.pdf", width=9)
g = ggbarplot(sum, "V1", "N",
  fill = "V2", color = "V2", palette = "Paired",
  label = TRUE, lab.col = "white", lab.pos = "in")
ggpar(g,font.tickslab = c(10,"plain", "black"),
 xtickslab.rt = 45)
dev.off()


#4. summary FDR genes
#geom_tile
#x-axis ion channel
#y-axis cancer type
#fill=HR 

fdr_sig$type = factor(fdr_sig$type, levels = canc_order$V1)
fdr_sig$HR = log2(fdr_sig$HR)
z = which(is.na(fdr_sig$type))
if(!(length(z)==0)){
  fdr_sig = fdr_sig[-z,]
}

canc_order = as.data.table(table(fdr_sig$name))
canc_order = canc_order[order(-N)]
canc_order=as.data.table(filter(canc_order, N >=1))
fdr_sig$name = factor(fdr_sig$name, levels = canc_order$V1)

pdf("onco_channels_summary_prog_fdr0.05_geom_tile_genes.pdf", width=9, height=6)

g = ggplot(fdr_sig, aes(name, type)) +
  geom_tile(aes(fill = HR)) +
    scale_fill_gradient2(low = "blue", midpoint = 0, high = "red") +
    xlab("Ion Channel") + ylab("Cancer type") + theme_bw()+
    theme(legend.position="bottom")
ggpar(g,
 font.tickslab = c(4,"plain", "black"),
 xtickslab.rt = 70, legend.title="Hazard Ratio")
dev.off()













