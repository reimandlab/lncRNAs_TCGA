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
canc_dats = as.data.table(filter(canc_dats, fdr <=0.05))
canc_dats$risk_type = ""
canc_dats$risk_type[canc_dats$mean_diff <0] = "NonRisk"
canc_dats$risk_type[canc_dats$mean_diff >0] = "Risk"

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

lncorder = summary$lnc

##---------Main plot with barplot risk vs non risk -----------------------

g = ggbarplot(summary, x="lnc", y="NumPCGs", color="Risk") + 
	coord_flip() +
    theme(axis.title.y=element_blank())

g = ggpar(g, font.ytickslab = c(3, "plain", "black"), legend = "right") 
main_plot = g

##---------Covariate heatmap for cancer type-------------------------------------------

#color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
#mypal = sample(color, 23)
mypal = readRDS(file="palette_23_cancer_types.rds")

cancers_conv = rna[,which(colnames(rna) %in% c("type", "Cancer"))]
cancers_conv = cancers_conv[!duplicated(cancers_conv), ]
colnames(cancers_conv)[2] = "canc"
summary = merge(summary, cancers_conv, by="canc")

cancers = ggplot(summary, aes(lnc, 0.2)) +
    geom_tile(aes(fill = canc)) + geom_text(aes(label = type), size=1.5) +
    theme_classic() + scale_fill_manual(values=mypal) +
    coord_flip()

cancers = ggpar(cancers, legend = "none") + theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

#combine
cancers + main_plot + plot_layout(ncol = 2, widths = c(1, 12)) 


























