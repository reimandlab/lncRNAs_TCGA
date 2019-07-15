options(stringsAsFactors=F)

#load libraries 
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table", 
              "plyr", "ggpubr",
              "ggrepel", "stringr", "maftools", "magrittr", "ggExtra", "cowplot")
lapply(packages, require, character.only = TRUE)

date = Sys.Date()


#figure 2 
res = readRDS("Documents/final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
res = as.data.table(filter(res, data=="TCGA"))
canc_conv = readRDS("Documents/canc_conv.rds")
colnames(canc_conv)[2] = "cancer"
z=which(res$gene_name == "HOXA-AS4")
res$gene_name[z] = "HOXA10-AS"

res = merge(res, canc_conv, by="cancer")
t = as.data.table(table(res$type))
t = t[order(-N)]
res$type = factor(res$type, levels =t$V1)
res = res[order(type, HR)]

#need to order genes wtihin cancer type 
res$combo = paste(res$gene_name, res$type)
res = res[order(type, HR)]
order = unique(res$gene_name)

res$gene_name = factor(res$gene_name, levels=order)

#x = lnc
#stratify by type
#y-axis = HR
#colour = lncRNA type CAT_geneClass

#res = res[1:20,]
res$HR = as.numeric(res$HR)
res$HR = log2(res$HR)
res$fdr = -log10(res$fdr)

pdf("lncRNA_candidates_final_figure2B.pdf", width=14, height=6)

g = ggplot(data=res, aes(x=gene_name, y=HR, order = -HR)) + 
  geom_bar(stat="identity", aes(fill=fdr),colour="black") + facet_grid(~ type, scale="free", space = "free")+
  geom_hline(yintercept=0, linetype="dashed", color = "red") + theme_minimal()
ggpar(g, xtickslab.rt=90, font.tickslab=c(6, "plain", "black"),
      legend = "bottom", legend.title = "Wald test, -log10(FDR)",
      font.legend = c(10, "plain", "black")) + scale_fill_gradient(low = "white", high = "black")+
  theme(strip.text.x = element_text(size = 9, colour = "Black", angle=90), axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1)) + 
  xlab("lncRNA") + ylab("log2(Hazard Ratio)")

dev.off()



