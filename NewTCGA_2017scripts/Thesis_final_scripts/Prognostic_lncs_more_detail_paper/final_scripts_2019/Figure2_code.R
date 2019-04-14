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
require("powerSurvEpi")
library(SIBER)
library(EnvStats)

#------FEATURES-----------------------------------------------------

#cands -- should be this file
cands = readRDS("genes_keep_100CV_No_FDR_May2nd2018.rds")

#--------This script ------------------------------------------------

#just make KM plots for TCGA 
#whatever data is available for PCAWG
#make them KM plots as well 
#just get list of genes that are significant in both data sets
#also check Cox PH assumptions within each data-set
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

get_name = function(ensg){
    z = which(fantom$CAT_geneID == ensg)
    return(fantom$CAT_geneName[z][1])
}

#--------------------------------------------------------------------

res = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
res = as.data.table(filter(res, data=="TCGA"))
canc_conv = readRDS("canc_conv.rds")
colnames(canc_conv)[2] = "cancer"
res = merge(res, canc_conv, by="cancer")
t = as.data.table(table(res$type))
t = t[order(N)]
res$type = factor(res$type, levels =t$V1)
res = res[order(type, HR)]

#need to order genes wtihin cancer type 
res$combo = paste(res$CAT_geneName, res$type)
res = res[order(type, HR)]
order = unique(res$CAT_geneName)

res$CAT_geneName = factor(res$CAT_geneName, levels=order)

#x = lnc
#stratify by type
#y-axis = HR
#colour = lncRNA type CAT_geneClass

res = res[1:20,]

res$HR = log2(res$HR)

g = ggplot(data=res, aes(x=CAT_geneName, y=HR, fill=CAT_geneClass, order = -HR)) + 
  geom_bar(stat="identity") + facet_grid(~ type, scale="free", space = "free")+
  geom_hline(yintercept=0, linetype="dashed", color = "red") + theme_light()
ggpar(g, xtickslab.rt=45, font.tickslab=c(6, "plain", "black"),
	legend = "bottom", legend.title = "lncRNA type",
 font.legend = c(5, "plain	", "black")) + scale_fill_brewer(palette="Dark2")








