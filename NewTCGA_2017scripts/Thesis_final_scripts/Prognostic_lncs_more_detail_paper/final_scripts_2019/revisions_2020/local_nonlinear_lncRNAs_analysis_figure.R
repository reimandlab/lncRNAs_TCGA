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

#make histogram for percent of patients with zero expression
res$perc_zeroes = as.numeric(res$perc_zeroes)
res$median_exp = as.numeric(res$median_exp)
res$median_nonzero = as.numeric(res$median_nonzero)

nonlinear = filter(res, median_exp == 0)
print(dim(nonlinear))
#[1] 97 31

#% non-silent tumors,
#hazard ratio,
#median RNA abundance in expressed-group used in X-axis, Y-axis, dot size, dot color gradient
nonlinear$hazard = ""
nonlinear$hazard[nonlinear$HR >0] = "HighExpHazard"
nonlinear$hazard[nonlinear$HR <0] = "LowExpHazard"
nonlinear$median_nonzero = log1p(nonlinear$median_nonzero)
nonlinear$name = ""
nonlinear$name[nonlinear$gene_symbol == "HOXA10-AS"] = "HOXA10-AS"

#colours
colours_palette=readRDS("23_cancers_color_palette.rds")

myColors=colours_palette$color
names(myColors)=colours_palette$cancer
myColors
colScale <- scale_colour_manual(name = "type",values = myColors)
colScale <- scale_fill_manual(name = "type",values = myColors)

pdf("Nonlinear_lncRNAs_hazard_ratio_vs_non_zero_median_expression.pdf", height=5, width=6)
ggscatter(nonlinear, x = "HR", y = "median_nonzero", label = "name",
   shape=21, fill = "type", size="perc_zeroes") + ylab("log1p(FPKM-UQ)")+
   xlab("log2(HR)")+colScale+guides(fill=FALSE) + ylim(c(0, max(nonlinear$median_nonzero)))
dev.off()
