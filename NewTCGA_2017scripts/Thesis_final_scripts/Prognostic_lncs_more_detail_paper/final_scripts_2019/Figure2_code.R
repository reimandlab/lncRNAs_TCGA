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
t = t[order(-N)]
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

#res = res[1:20,]

res$HR = log2(res$HR)
res$fdr = -log10(res$fdr)
res$fdr = as.numeric(res$fdr)

pdf("lncRNA_candidates_final_figure2B.pdf", width=15, height=6)

g = ggplot(data=res, aes(x=CAT_geneName, y=HR, order = -HR)) + 
  geom_bar(stat="identity", aes(fill=fdr)) + facet_grid(~ type, scale="free", space = "free")+
  geom_hline(yintercept=0, linetype="dashed", color = "red") + theme_minimal()
ggpar(g, xtickslab.rt=90, font.tickslab=c(7, "plain", "black"),
	legend = "bottom", legend.title = "Wald test, adjusted \n-log10(p-value)",
 font.legend = c(10, "plain	", "black")) + scale_fill_gradient(low = "black", high = "white")+
theme(strip.text.x = element_text(size = 8, colour = "Black", angle=90)) + xlab("lncRNA") + ylab("log2(Hazard Ratio)")

dev.off()






