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

coexp = readRDS("all_results_for_each_cancer_from_coexpression_analysis_june27th_allCands.rds")

#can run FDR on all PCGs for all Cancer types 
#175 * 20,0000 = enormous FDR test 

#1. Get FDR by cancer type 
cancers = as.list(unique(coexp$canc))
canc_fdr = function(cancer){
	dat = subset(coexp, canc == cancer)
	dat$fdr = p.adjust(dat$pvalue, method = "fdr")
	dat = as.data.table(dat)
	dat = dat[order(fdr)]
	return(dat)
}

library(dplyr)
library(plyr)
canc_dats = llply(cancers, canc_fdr) 
canc_dats = ldply(canc_dats, data.frame)

#2. Summarize per lncRNA/cancer, how many PCGs upregulated in risk group
#and how many upregulated in non-risk group 
#(>0 --> more expressed in risk group, <0, more expressed in low risk group)

#keep sig fdr, add risk/non-risk tag
canc_dats = as.data.table(canc_dats)
canc_dats = as.data.table(filter(canc_dats, fdr <=0.01))
canc_dats$risk_type = ""
canc_dats$risk_type[canc_dats$mean_diff < (log1p(2)-log1p(4))] = "NonRisk"
canc_dats$risk_type[canc_dats$mean_diff >= (log1p(4)-log1p(2))] = "Risk"
canc_dats = canc_dats[which(!(canc_dats$risk_type =="")),]

#fold change of two equals to 0.51
#fold change of 1/2 equals to -0.51

summary = as.data.table(table(canc_dats$lnc, canc_dats$risk_type))
summary = filter(summary, N >0)
summary = as.data.table(summary)
colnames(summary) = c("lnc", "Risk", "NumPCGs")

#re-order
summary = summary[order(NumPCGs)]
summary$canc = ""
for(i in 1:nrow(summary)){
	canc = unique(canc_dats$canc[which(canc_dats$lnc %in% summary$lnc[i])])
	summary$canc[i] = canc
}

lncorder = unique(summary$lnc)
summary$lnc = factor(summary$lnc, levels = lncorder)

#Change lnc IDs to gene names 
summary$name = ""
for(i in 1:nrow(summary)){
  z1 = which(allCands$gene == unlist(summary[i,1]))
  z2 = which(allCands$combo == unlist(summary[i,1]))
  name = allCands$CAT_geneName[c(z1,z2)]
  summary$name[i] = name
}

lncorder = unique(summary$name)
summary$name = factor(summary$name, levels = lncorder)

#add tag to non-unique ones
cands_dups = unique(allCands$CAT_geneName[which(duplicated(allCands$CAT_geneName))])
z = which(summary$name %in% cands_dups)
summary$name = as.character(summary$name)
for(y in 1:length(z)){
  print(y)
  g = paste(summary$name[z[y]], y, sep="_")
  summary$name[z[y]] = g
}

####Need to re-level cancer types - not in right order!!!! fixed? 

#color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
#mypal = sample(color, 23)
mypal = readRDS(file="palette_23_cancer_types.rds")

cancers_conv = rna[,which(colnames(rna) %in% c("type", "Cancer"))]
cancers_conv = cancers_conv[!duplicated(cancers_conv), ]
colnames(cancers_conv)[2] = "canc"
summary = merge(summary, cancers_conv, by="canc")
summary = summary[order(Risk, NumPCGs)]
lncorder = unique(summary$lnc)
summary$lnc = factor(summary$lnc, levels = lncorder)
#canc_order = unique(summary$type)
#summary$type = factor(summary$type, levels = canc_order)
lncorder = unique(summary$name)
summary$name = factor(summary$name, levels = lncorder)
riskorder = c("Risk", "NonRisk")
summary$Risk = factor(summary$Risk, levels = riskorder)

##---------Main plot with barplot risk vs non risk -----------------------

ratios = ggplot(data=summary, aes(x=name, y=NumPCGs, color=Risk)) +
  geom_bar(stat="identity") +
  theme_bw() + 
  coord_flip() +
  theme(axis.title.y=element_blank())

ratios = ggpar(ratios, legend = "right", font.ytickslab=c(2, "plain", "black"))

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
z1 = which(allCands$gene %in% summary$lnc)
z2 = which(allCands$combo %in% summary$lnc)
nosig = allCands[-c(z1, z2),]

#this one is only one that didnt have any
#ENSG00000240889 0.933567971848771 2.54356838690527 0.0008554614
#              low95          upper95                  cancer     fdr_pval data
#1: 1.46929431922299 4.40329759274187 Glioblastoma multiforme


##---------What kinds of genes are enriched per lncRNA-----------------------------------

cancer_genes = fread("list_of_cancer_genes_mmc10-2.txt", header=F)
rbps = fread("list_of_RBPs_mmc2-5.txt", header=F)
TFS = fread("list_of_TFs_mmc2-5.txt", header=F)
TF_targets = fread("TF_targets_mmc5-2.txt", header=F)
colnames(TF_targets) = c("TF", "target")
protein_at = fread("proteinatlas.tsv")
sub_loc = fread("protein_atlas_subcellular_location.tsv")
#keep only reliable
sub_loc = as.data.table(filter(sub_loc, Reliability == "Approved"))
patho = fread("protein_atlas_pathology.tsv")

#how many duplicated PCGs
pcgs_sum = as.data.table(table(canc_dats$pcg, canc_dats$risk_type))
pcgs_sum = as.data.table(filter(pcgs_sum, N >1))
pcgs_sum = pcgs_sum[order(N)]
pcgs_sum$both_risk_groups = ""

#which pcgs both in risk and non-risk
both_risks = as.data.table(table(pcgs_sum$V1))
both_risks = as.data.table(filter(both_risks, N >1))

pcgs_sum$both_risk_groups = ""
z = which(pcgs_sum$V1 %in% both_risks$V1)
#7212/19175 in both risk and nonrisk
pcgs_sum$both_risk_groups[z] = "BOTH_risks"

#plot summary how many cancer types pcg
#is in high risk group

risk = as.data.table(filter(pcgs_sum, V2 == "Risk"))
risk$groupy = cut(risk$N, breaks =c(1,5, 10, 15, 20,25, 30))

nonrisk = as.data.table(filter(pcgs_sum, V2 == "NonRisk"))
nonrisk$groupy = cut(nonrisk$N, breaks =c(1,5, 10, 15, 20,25, 30))

#SUMMARIZE
# Change line color and fill color
riskplot = ggplot(risk, aes(x=groupy))+
  geom_histogram(color="darkblue", fill="darksalmon", stat="count") + theme_bw()

nonriskplot= ggplot(nonrisk, aes(x=groupy))+
  geom_histogram(color="darkblue", fill="lightblue", stat="count") + theme_bw()

#risk + nonrisk

both = rbind(risk, nonrisk)

pdf("summary_pcgs_coexpressed_mutliple_cancer_types.pdf")
ggplot(both, aes(x=groupy, fill=V2)) +
  scale_fill_brewer(palette="Dark2") +
  geom_histogram(position="dodge", stat="count", alpha=0.8)+
  theme(legend.position="top") +
  theme_bw() +
  xlab("Number of Cancers") + ylab("Number of PCGs")
dev.off()

#distribution of PCGs that appeared in both groups
bothrisks = as.data.table(filter(both, both_risk_groups=="BOTH_risks"))
pdf("just_pcgs_both_risks_summary_pcgs_coexpressed_mutliple_cancer_types.pdf")
ggplot(bothrisks, aes(x=groupy, fill=V2)) +
  scale_fill_brewer(palette="Dark2") +
  geom_histogram(position="dodge", stat="count", alpha=0.8)+
  theme(legend.position="top") +
  theme_bw() +
  xlab("Number of Cancers") + ylab("Number of PCGs")
dev.off()

#4229 appear in both risk and non-risk groups (across all cancers)

##---------pathways enriched by PCGs that appear in at least 2 cancers and in risk group
#remove those that are in both risk and non-risk ...? 

risk = as.data.table(filter(risk, both_risk_groups == "", N >=10))
genes = risk$V1
combined_paths <- gprofiler(genes, organism = "hsapiens", exclude_iea=TRUE, ordered_query= TRUE, min_set_size=5, max_set_size = 200, min_isect_size=2, correction_method="fdr")
print(dim(combined_paths)[1])

if(!(dim(combined_paths)[1]==0)){
#only keep GO or REACTOME
reac <- grep("REAC", combined_paths$term.id)
go <- grep("GO", combined_paths$term.id)
combined_paths <- combined_paths[c(reac, go), ]
combined_paths <- combined_paths[,c(9,12, 3, 3, 1, 14)]
colnames(combined_paths) <- c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Genes")
combined_paths$Phenotype[combined_paths$Phenotype==1] = "1"
combined_paths$Phenotype[combined_paths$Phenotype==2] = "-1"
write.table(combined_paths, sep= "\t", file=paste(colnames(d)[6], d$canc[1], "PathwaysUsingtALL_DEgenesFeb16_linearLASSO.txt", sep="_"), quote=F, row.names=F)
}

##---------pathways enriched by PCGs that appear in at least 2 cancers and in non-risk group
#remove those that are in both risk and non-risk ...? 

nonrisk = as.data.table(filter(nonrisk, both_risk_groups == "", N >=10))
genes = nonrisk$V1










