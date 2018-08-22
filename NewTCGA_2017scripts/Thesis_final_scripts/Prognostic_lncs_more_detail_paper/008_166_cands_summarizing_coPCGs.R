library(survAUC)
#source("source_code_Cox_MonteCarlo_CV_Mar13.R")
require(caTools)
#check if this person is in my analysis: TCGA-61-2095
library(glmnet)
library(survcomp)
library(caret)
library(stringr)
library(EnvStats)
library(patchwork)

source("universal_LASSO_survival_script.R")

library(ggpubr)
library(ggrepel)
library(viridis)
library(patchwork)
library(caret)  
library(Rtsne)
library(data.table)

#------DATA---------------------------------------------------------
#UCSC gene info
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
#z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
#ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,8]))
ucsc <- ucsc[-z,]

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

#save RNA and PCG files locally
#saveRDS(rna, file="rna_lncRNAs_expression_data_june29.rds")
#saveRDS(pcg, file="rna_pcg_expression_data_june29.rds")

rna = readRDS("rna_lncRNAs_expression_data_june29.rds")
pcg = readRDS("rna_pcg_expression_data_june29.rds")

#------FEATURES-----------------------------------------------------

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = subset(allCands, data == "TCGA") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 
allCands$combo = unique(paste(allCands$gene, allCands$cancer, sep="_"))

val_cands = read.csv("175_lncRNA_cancers_combos_23_cancer_types_july5.csv")
val_cands = as.data.table(val_cands)
val_cands = subset(val_cands, data == "PCAWG") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 
val_cands$combo = unique(paste(val_cands$gene, val_cands$cancer, sep="_"))
val_cands = subset(val_cands, top_pcawg_val == "YES") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 

#--------This script ------------------------------------------------

#summarize results from co-expression analysis of PCGs
#how many per lcnRNA
#how many pathways per lncRNA
#how many cancer genes, RBPs, TFs...

#--------------------------------------------------------------------

#write function that adds tag to whole data group 
#and does survival analysis on whole group

cands_dups = unique(allCands$gene[which(duplicated(allCands$gene))])

#--------------------------------------------------------------------
#RESULTS-------------------------------------------------------------
#--------------------------------------------------------------------

#coexp = readRDS("new_PCG_enrichment_min_5FPKMmedian_july20.rds")
#new file using diff exp to get PCGs
coexp = readRDS("diff_expressed_PCGs_lncRNA_risk_groups_Aug21.rds")
#summary how many PCGs per lncRNA
coexp$combo = paste(coexp$lnc, coexp$cancer)
sum = as.data.table(table(coexp$combo))
sum
sum = sum[order(N)]
sum$V1 = factor(sum$V1, levels = unique(sum$V1))
pdf("summary_lncRNA_DE_PCGs_148_pairs_aug21.pdf", height=12)
p<-ggplot(data=sum, aes(x=V1, y=N)) +
  geom_bar(stat="identity", color="blue", fill="white")  
# Horizontal bar plot
p + coord_flip() + ylab("Number of sig DE PCGs") + xlab("lncRNA-cancer pair")+
geom_vline(xintercept = 500, linetype="dashed", color = "red")+
theme(axis.text.y = element_text(size=4)) 

dev.off()

#z = which(coexp$lnc %in% cands_dups)
#if(!(length(z))==0){
#coexp$lnc[z] = paste(coexp$lnc[z], coexp$canc[z], sep="_")
#}

#coexp$combo=""
#z = which(str_detect(coexp$lnc, "_"))
#coexp$combo[-z] = paste(coexp$lnc[-z], coexp$canc[-z], sep="_")
#coexp$combo[z] = coexp$lnc[z]
#there was one extra combination from before that shouldn't be part of this
#z = which(coexp$combo %in% allCands$combo)
#coexp = coexp[z,]

#can run FDR on all PCGs for all Cancer types 
#175 * 20,0000 = enormous FDR test 

#1. Get FDR by cancer type 
#cancers = as.list(unique(coexp$canc))
#canc_fdr = function(cancer){
#  dat = subset(coexp, canc == cancer)
#	dat$pvalue = as.numeric(dat$pvalue)
#  dat$fdr = p.adjust(dat$pvalue, method = "fdr")
#	z = which(is.na(dat$pvalue))
#  if(!length(z)==0){
#    dat = dat[-z,]
#  }
#  dat = as.data.table(dat)
#	dat = dat[order(fdr)]
#	return(dat)
#}

library(dplyr)
library(plyr)
#canc_dats = llply(cancers, canc_fdr) 
#canc_dats = ldply(canc_dats, data.frame)

#plot scatter plot, FC versus p-value
#coexp$fdr = -log10(coexp$fdr)
#canc_dats$mean_diff = as.numeric(canc_dats$mean_diff)
#z = which(canc_dats$mean_diff == 0)
#canc_dats = canc_dats[-z,]
#canc_dats = as.data.table(canc_dats)

#2. Summarize per lncRNA/cancer, how many PCGs upregulated in risk group
#and how many upregulated in non-risk group 
#(>0 --> more expressed in risk group, <0, more expressed in low risk group)

#keep sig fdr, add risk/non-risk tag
coexp = as.data.table(coexp)
coexp$combo = paste(coexp$lnc, coexp$cancer, sep="_")

coexp$risk_type = ""
#add lncRNA - risk status whether its high exp or low exp
get_risk = function(combo){
  hr = as.numeric(allCands$HR[which(allCands$combo == combo)])
  if(hr >1){
    risk = "high_exp"
  }
  if(hr <1){
    risk = "low_exp"
  }
  return(risk)
}
coexp$risk_type = llply(coexp$combo, get_risk, .progress="text")
coexp$risk_type = unlist(coexp$risk_type)

coexp$pcg_risk[(coexp$logFC >= log(2))&(coexp$risk_type=="high_exp")] = "upregulated_in_risk"
coexp$pcg_risk[(coexp$logFC <= log(0.5))&(coexp$risk_type=="low_exp")] = "upregulated_in_risk"
coexp$pcg_risk[(coexp$logFC >= log(2))&(coexp$risk_type=="low_exp")] = "downregulated_in_risk"
coexp$pcg_risk[(coexp$logFC <= log(0.5))&(coexp$risk_type=="high_exp")] = "downregulated_in_risk"

canc_dats = coexp 

#fold change of two equals to 0.51
#fold change of 1/2 equals to -0.51
summary = as.data.table(table(canc_dats$combo, canc_dats$pcg_risk))
summary = filter(summary, N >0)
summary = as.data.table(summary)
colnames(summary) = c("combo", "Risk", "NumPCGs")

#re-order
summary = summary[order(NumPCGs)]
summary$canc = ""
for(i in 1:nrow(summary)){
	print(i)
  canc = unique(canc_dats$canc[which(canc_dats$combo %in% summary$combo[i])])
	summary$canc[i] = canc
}

lncorder = unique(summary$combo)
summary$combo = factor(summary$combo, levels = lncorder)

#Change lnc IDs to gene names 
summary$name = ""
for(i in 1:nrow(summary)){
  z2 = which(allCands$combo == unlist(summary[i,1]))
  name = allCands$CAT_geneName[c(z2)]
  summary$name[i] = name
}

lncorder = unique(summary$name)
summary$name = factor(summary$name, levels = lncorder)

#############################################################
####keep those with at least 10 sig PCGs#####################
#############################################################

summary = as.data.table(filter(summary, NumPCGs >= 10))

####Need to re-level cancer types - not in right order!!!! fixed? 

#color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
#mypal = sample(color, 23)
mypal = readRDS(file="palette_23_cancer_types.rds")

cancers_conv = rna[,which(colnames(rna) %in% c("type", "Cancer"))]
cancers_conv = cancers_conv[!duplicated(cancers_conv), ]
colnames(cancers_conv)[2] = "canc"
saveRDS(cancers_conv, file="cancers_conv_july23.rds")
summary = merge(summary, cancers_conv, by="canc")
summary = summary[order(Risk, NumPCGs)]
lncorder = unique(summary$name)
summary$name = factor(summary$name, levels = lncorder)
#canc_order = unique(summary$type)
#summary$type = factor(summary$type, levels = canc_order)
lncorder = unique(summary$name)
summary$name = factor(summary$name, levels = lncorder)
riskorder = c("upregulated_in_risk", "downregulated_in_risk")
summary$Risk = factor(summary$Risk, levels = riskorder)

##---------Main plot with barplot risk vs non risk -----------------------

#for plot only keep those lnc-canc combos with at least 50 pcgs 
summary = as.data.table(summary)

#add the candidates that validated
z = which(summary$combo %in% val_cands$combo)
summary$val = ""
summary$val[z] = "*"
summary$val[-z] = ""

ratios = ggplot(data=summary, aes(x=name, y=NumPCGs, color=Risk)) +
  geom_bar(stat="identity") + geom_text(aes(label = summary$val), color="black", size=5) + 
  theme_bw() + 
  coord_flip() +
  theme(axis.title.y=element_blank())

ratios = ggpar(ratios, legend = "right", font.ytickslab=c(4, "plain", "black"))

##---------Covariate heatmap for cancer type-------------------------------------------

#just canc data 174 rows instead of 348
z = which(duplicated(summary$name))
justcanc = summary[-z,]
justcanc$name = factor(justcanc$name, levels = lncorder)

cancers = ggplot(justcanc, aes(name, 0.2)) +
    geom_tile(aes(fill = type)) + geom_text(aes(label = type), size=1.5) +
    theme_classic() + scale_fill_manual(values=mypal) +
    coord_flip()

cancers = ggpar(cancers, legend = "none") + theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

#combine
pdf("pcg_summary_lncrna_cand_risk_nonrisk_july3.pdf")
cancers + ratios + plot_layout(ncol = 2, widths = c(1, 12)) 
dev.off()

#summarize which ones didn't have any sig pcgs
z2 = which(allCands$combo %in% summary$combo)
nosig = allCands[-c(z2),]

##--------------------------------------------------------------------------------------
#subset canc_dats to only include lncRNA-cancer-PCGs as in summary
##--------------------------------------------------------------------------------------

z = which(canc_dats$combo %in% summary$combo)
canc_dats = canc_dats[z,]

####-------SAVE PROCESSED CO-EXPRESSION RESULTS-----------------------------------------

saveRDS(canc_dats, file="coexpression_results_processed_july24.rds")

##---------What kinds of genes are enriched per lncRNA-----------------------------------

cancer_genes = fread("list_of_cancer_genes_mmc10-2.txt", header=F)
colnames(cancer_genes)[1] = "name"
colnames(ucsc)[6] = "V1"
colnames(ucsc)[8] = "name"
cancer_genes = merge(cancer_genes, ucsc, by="name")

rbps = fread("list_of_RBPs_mmc2-5.txt", header=F)
colnames(rbps)[1] = "name"
rbps = merge(rbps, ucsc, by="name")

TFS = fread("list_of_TFs_mmc2-5.txt", header=F)
colnames(TFS)[1] = "name"
TFS = merge(TFS, ucsc, by="name")

TF_targets = fread("TF_targets_mmc5-2.txt", header=F)
colnames(TF_targets) = c("TF", "target")

protein_at = fread("proteinatlas.tsv")
sub_loc = fread("protein_atlas_subcellular_location.tsv")
#keep only reliable
sub_loc = as.data.table(filter(sub_loc, Reliability == "Approved"))
patho = fread("protein_atlas_pathology.tsv")

census = read.csv("Census_allFri_Jul_13_16_55_59_2018.csv")

#how many duplicated PCGs
#pcgs_sum = as.data.table(table(canc_dats$pcg, canc_dats$risk_type, canc_dats$canc))
pcgs_sum = as.data.table(table(canc_dats$ID, canc_dats$pcg_risk))
pcgs_sum = as.data.table(filter(pcgs_sum, N >=1))
pcgs_sum = pcgs_sum[order(N)]
pcgs_sum$both_risk_groups = ""

#which pcgs both in risk and non-risk
both_risks = as.data.table(table(pcgs_sum$V1))
both_risks = as.data.table(filter(both_risks, N >1))

pcgs_sum$both_risk_groups = ""
z = which(pcgs_sum$V1 %in% both_risks$V1)
# 6758/19029 pcgs appear in both risk groups 
pcgs_sum$both_risk_groups[z] = "BOTH_risks"
z = which(pcgs_sum$both_risk_groups == "BOTH_risks")
both_risks = unique(pcgs_sum$V1[z])
#what type of PCGs are they?
#875/4229 are cancer genes
#1/4229 RBP
#172/4229 are Transcription factors 

#plot summary how many cancer types pcg
#is in high risk group
risk = as.data.table(filter(pcgs_sum, V2 == "upregulated_in_risk"))
risk$groupy = cut(risk$N, breaks =c(0,1,5, 10, 15, 20,25, 31))

nonrisk = as.data.table(filter(pcgs_sum, V2 == "downregulated_in_risk"))
nonrisk$groupy = cut(nonrisk$N, breaks =c(0,1,5, 10, 15, 20,25))

#SUMMARIZE
# Change line color and fill color
riskplot = ggplot(risk, aes(x=groupy))+
  geom_histogram(color="darkblue", fill="darksalmon", stat="count") + theme_bw()

nonriskplot= ggplot(nonrisk, aes(x=groupy))+
  geom_histogram(color="darkblue", fill="lightblue", stat="count") + theme_bw()

#risk + nonrisk

both = rbind(risk, nonrisk)

pdf("summary_pcgs_coexpressed_mutliple_cancer_types_24.pdf", width=5, height=5)
g = ggplot(both, aes(x=groupy, fill=V2)) +
  scale_fill_brewer(palette="Dark2") +
  geom_histogram(position="dodge", stat="count", alpha=0.8)+
  theme(legend.position="top") +
  theme_bw() +
  xlab("Number of lncRNA candidates") + ylab("Number of PCGs")
ggpar(g, legend.title="lncRNA group")
dev.off()

#distribution of PCGs that appeared in both groups
bothrisks = as.data.table(filter(both, both_risk_groups=="BOTH_risks"))
pdf("just_pcgs_both_risks_summary_pcgs_coexpressed_mutliple_cancer_types_24.pdf", width=5, height=5)
g = ggplot(bothrisks, aes(x=groupy, fill=V2)) +
  scale_fill_brewer(palette="Dark2") +
  geom_histogram(position="dodge", stat="count", alpha=0.8)+
  theme(legend.position="top") +
  theme_bw() +
  xlab("Number of Cancers") + ylab("Number of PCGs")
ggpar(g, legend.title="lncRNA group")
dev.off()

#4229 appear in both risk and non-risk groups (across all cancers)

#Save files for gprofiler 
saveRDS(both, file="pcgs_enriched_in_risk_groups_non_lncRNA_risk_groups_pcg_analysis_july24.rds")

#---------BY CANCER TYPE ANALYSIS-----------------------------------------------------------------

#PCGs enriched in risk and non-risk by cancer type 
pcgs_sum = as.data.table(table(canc_dats$canc, canc_dats$ID, canc_dats$pcg_risk))
pcgs_sum = as.data.table(filter(pcgs_sum, N >=1))
pcgs_sum = pcgs_sum[order(N)]
colnames(pcgs_sum) = c("cancer", "pcg", "group", "NumTimes")

#summarize how many PCGs appear in all lncRNA risk groups 
#how many cands are there per each cancer type?
lncs_canc = as.data.table(table(allCands$cancer))
colnames(lncs_canc) = c("cancer", "Num_lncs_cands")

pcgs_sum = merge(pcgs_sum, lncs_canc, by="cancer")
pcgs_sum = pcgs_sum[order(NumTimes)]

#add gene names
colnames(ucsc)[6] = "pcg"
pcgs_sum = merge(pcgs_sum, ucsc, by="pcg")
pcgs_sum = pcgs_sum[order(NumTimes)]
colnames(cancers_conv)[2] = "cancer"
pcgs_sum = merge(pcgs_sum, cancers_conv, by="cancer")
pcgs_sum = pcgs_sum[order(NumTimes)]
pcgs_sum$frac_cands = pcgs_sum$NumTimes/pcgs_sum$Num_lncs_cands

#plot just the ones that appear in all pcgs risks groups 
#and nonrisk groups
pcgs_sum = filter(pcgs_sum, NumTimes == Num_lncs_cands)
pcgs_sum = as.data.table(pcgs_sum)

#dont really know what ^ is for

##---------HRs versus # of PCGS-------------------------------------------------------------------
res_tog = merge(allCands, summary, by= "combo")
res_tog$HR = as.numeric(res_tog$HR)
res_tog$HR = log2(res_tog$HR)
mypal = readRDS(file="palette_32_cancer_types.rds")

# Basic scatter plot
pdf("summary_coexpressed_risk_non_risk_wHR_500.pdf")
ggplot(res_tog, aes(x=NumPCGs, y=HR, shape=Risk)) + geom_point(size=1.5) +
geom_hline(yintercept = 0, linetype="dashed", color = "red") + 
geom_vline(xintercept = 500, linetype="dashed", color = "red") +
theme_bw()+
geom_label_repel(data=filter(res_tog, NumPCGs >=500), aes(label=name, fill=type), color = 'darkslategrey',size=1)+
scale_fill_manual(values=mypal)
dev.off()

pdf("summary_coexpressed_risk_non_risk_wHR_750.pdf")
ggplot(res_tog, aes(x=NumPCGs, y=HR, shape=Risk)) + geom_point(size=1.5) +
geom_hline(yintercept = 0, linetype="dashed", color = "red") + 
geom_vline(xintercept = 750, linetype="dashed", color = "red") +
theme_bw()+
geom_label_repel(data=filter(res_tog, NumPCGs >=750), aes(label=name, fill=type), color = 'black',size=2)+
scale_fill_brewer(palette="Paired")
dev.off()

#summarize # of PCGs DE in favourable lncRNAs versus unfavourable
res_tog$lnc_risk[res_tog$HR > log2(1)] = "Unfavourable"
res_tog$lnc_risk[res_tog$HR < log2(1)] = "Favourable"
t = table(res_tog$lnc_risk, res_tog$Risk)
chisq.test(t)
#no sig association

#add layer about risk, are they all "binary" lncs or "balanced" lncs?

res_tog$perc_risk = as.numeric(res_tog$perc_risk)
res_tog$perc_risk_label = ""
res_tog$perc_risk_label[which(res_tog$perc_risk < 0.49)] = "SmallRisk"
res_tog$perc_risk_label[which(res_tog$perc_risk > 0.51)] = "HighRisk"
res_tog$perc_risk_label[which(res_tog$perc_risk_label == "")] = "BalRisk"

pdf("summary_coexpressed_risk_non_risk_wHR_750_wcolor.pdf")
ggplot(res_tog, aes(x=NumPCGs, y=HR, shape=Risk)) + geom_point(size=1.5) +
geom_hline(yintercept = 0, linetype="dashed", color = "red") + 
geom_vline(xintercept = 750, linetype="dashed", color = "red") +
theme_bw()+
geom_label_repel(data=filter(res_tog, NumPCGs >=750), aes(label=name, fill=type, color = perc_risk_label), size=2)+
scale_fill_brewer(palette="Paired") +
scale_color_manual(values=c("black", "blue", "darkslategrey"))

dev.off()

res_tog$HR = as.numeric(res_tog$HR)
res_tog$lnc_stat = ""
res_tog$lnc_stat[which(res_tog$HR < 0)] = "Favourable"
res_tog$lnc_stat[which(res_tog$HR > 0)] = "Unfavourable"

saveRDS(res_tog, file="summary_pcg_analysis_wHRs_jul2y24.rds")

#144/173 have significant pcg signatures 









