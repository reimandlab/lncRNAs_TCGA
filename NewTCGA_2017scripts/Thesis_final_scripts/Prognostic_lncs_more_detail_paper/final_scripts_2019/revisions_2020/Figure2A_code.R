options(stringsAsFactors=F)

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table",
              "plyr", "ggpubr",
              "ggrepel", "stringr", "maftools", "magrittr", "ggExtra", "cowplot")
lapply(packages, require, character.only = TRUE)

date = Sys.Date()

setwd("Documents/lncRNAs/Jan2021")

#figure 2 - univaraite
res = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
res = as.data.table(filter(res, data=="TCGA"))
write.csv(res, file="main_170_survival_annotations_file.csv", quote=F, row.names=F)

t = as.data.table(table(res$canc_type))
t = t[order(-N)]
res$type = factor(res$canc_type, levels =t$V1)
res$HR=as.numeric(res$HR)
res = res[order(type, HR)]

#need to order genes wtihin cancer type
res$combo = paste(res$gene_symbol, res$canc_type)
res = res[order(canc_type, HR)]

order = unique(res$gene_symbol)
res$gene_symbol = factor(res$gene_symbol, levels=order)

order = unique(res$combo)
res$combo = factor(res$combo, levels=order)

#x = lnc
#stratify by type
#y-axis = HR
#colour = lncRNA type CAT_geneClass

#res = res[1:20,]
res$HR = as.numeric(res$HR)
res$HR = log2(res$HR)
res$fdr_pval = -log10(res$fdr_pval)
res$low95 = as.numeric(res$low95)
res$low95 = log2(res$low95)
res$upper95 = as.numeric(res$upper95)
res$upper95 = log2(res$upper95)

pdf("Figure2_related/Figure2A_unadjusted_HR.pdf", width=14, height=6)

g = ggplot(data=res, aes(x=combo, y=HR, order = -HR)) +
  geom_bar(stat="identity", aes(fill=fdr_pval),colour="black") +
  geom_errorbar(aes(ymin=low95, ymax=upper95),
                  width=.25,  size=.2,                   # Width of the error bars
                  position=position_dodge(.9), color="darkred") +
  facet_grid(~ type, scale="free", space = "free")+
  geom_hline(yintercept=0, color = "black") + theme_bw()
ggpar(g, xtickslab.rt=90, font.tickslab=c(6, "plain", "black"),
      legend = "bottom", legend.title = "Wald test, -log10(FDR)",
      font.legend = c(10, "plain", "black")) + scale_fill_gradient(low = "white", high = "black")+
  theme(strip.text.x = element_text(size = 9, colour = "Black", angle=90), axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1)) +
  xlab("lncRNA") + ylab("log2(Hazard Ratio)")

dev.off()

#make histogram for risk group per lncRNA
pdf("Figure2_related/all_lncRNAs_cands_risk_group_summary_histogram.pdf", width=5, height=5)
gghistogram(res, x="perc_risk", color="black", fill="#00AFBB") + xlim(0,1)+theme_bw()+
xlab("Percent of patients in risk group") + ylab("Number of lncRNAs")
dev.off()

#make histogram for percent of patients with zero expression
pdf("Figure2_related/all_lncRNAs_cands_percent_zero_expression.pdf", width=5, height=5)
res$perc_zeroes = as.numeric(res$perc_zeroes)
gghistogram(res, x="perc_zeroes", color="black", fill="#00AFBB") + xlim(0,1)+theme_bw()+
xlab("Percent of patients with zero expression") + ylab("Number of lncRNAs")
dev.off()

#barplot
res$perc_zeroes = as.numeric(res$perc_zeroes)
res$perc_on = 1 - res$perc_zeroes
res_off = res[,c("gene_symbol", "type", "perc_zeroes", "perc_on")]
res_off = melt(res_off)

pdf("fraction_lncRNA_turned_on_off.pdf", width=8,height=3)
g = ggbarplot(res_off, "gene_symbol", "value",
  fill = "variable", color = "variable",
  palette = c("#00AFBB", "#E7B800"))
g = ggpar(g,
 font.tickslab = c(5,"plain", "black"),
 xtickslab.rt = 90)+ylab("% of patients with no detected expression")
g + facet_grid(~type, scales = "free", space = "free") + theme(text = element_text(size=5),
axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=3))
dev.off()

#figure 2 - multivaraite
res = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
res = as.data.table(filter(res, data=="TCGA"))
t = as.data.table(table(res$canc_type))
t = t[order(-N)]
res$type = factor(res$canc_type, levels =t$V1)
res$HR = as.numeric(res$hr_adjusted)
res = res[order(type, HR)]

#need to order genes wtihin cancer type
res$combo = paste(res$gene_symbol, res$canc_type)
res = res[order(canc_type, HR)]

order = unique(res$gene_symbol)
res$gene_symbol = factor(res$gene_symbol, levels=order)

order = unique(res$combo)
res$combo = factor(res$combo, levels=order)

#x = lnc
#stratify by type
#y-axis = HR
#colour = lncRNA type CAT_geneClass

#res = res[1:20,]
res$HR = as.numeric(res$HR)
res$HR = log2(res$HR)
res$fdr_pval = -log10(res$fdr_pval_adjusted)
res$low95 = as.numeric(res$HR_adj_low95)
res$low95 = log2(res$low95)
res$upper95 = as.numeric(res$HR_adj_high95)
res$upper95 = log2(res$upper95)

pdf("Figure2_related/lncRNA_candidates_multivariate_figure2A.pdf", width=14, height=6)

g = ggplot(data=res, aes(x=combo, y=HR, order = -HR)) +
  geom_bar(stat="identity", aes(fill=fdr_pval),colour="black") +
  geom_errorbar(aes(ymin=low95, ymax=upper95),
                  width=.25,  size=.2,                   # Width of the error bars
                  position=position_dodge(.9), color="darkred") +
  facet_grid(~ type, scale="free", space = "free")+
  geom_hline(yintercept=0, color = "black") + theme_bw()
ggpar(g, xtickslab.rt=90, font.tickslab=c(6, "plain", "black"),
      legend = "bottom", legend.title = "Wald test, -log10(FDR)",
      font.legend = c(10, "plain", "black")) + scale_fill_gradient(low = "white", high = "black")+
  theme(strip.text.x = element_text(size = 9, colour = "Black", angle=90), axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1)) +
  xlab("lncRNA") + ylab("log2(Hazard Ratio)")

dev.off()
