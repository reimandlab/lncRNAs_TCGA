setwd("/Users/kisaev/remote10")

library(data.table)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(plyr)
library(EnvStats)
library(ggsci)

#summmary results 

res = readRDS("results_analysis_Feb26.rds")
res = as.data.table(ldply(res))
res$combo = paste(res$canc, res$tis)
res$combo2 = paste(res$gene, res$canc, sep="_")

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
#save only the ones that came from the noFDR appraoch 
allCands = filter(allCands, data=="TCGA", fdr_pval <=0.05)

#save pairs for sup table 3
ss3 = unique(res[,c("canc", "tis")])
write.csv(ss3, file="sup_table3_gtex_tis_cancers_studied.csv", quote=F, row.names=F)

#keep only sig ones
res$median_difference = as.numeric(res$median_difference)
res = as.data.table(filter(res, fdr < 0.05, abs(median_difference) >= 0.25))

#figure 1b 
res$sign = ""
res$sign[res$median_difference > 0] = "upreg"
res$sign[res$median_difference < 0] = "downreg"
fig = unique(res[,c("gene", "canc", "combo2", "sign")])

fig1 = as.data.table(table(fig$canc, fig$sign))
fig1 = fig1[order(N)]
summ = as.data.table(table(fig$canc))
summ = summ[order(N)]
fig1$V1 = factor(fig1$V1, levels= summ$V1)
fig1$V2 = factor(fig1$V2, levels = c("upreg", "downreg"))

#fig 1c----------------------------------------------
pdf("Manuscript_FIGURE1C.pdf", width=9)
g = ggbarplot(fig1, "V1", "N",
          fill = "V2", color = "V2", 
          label = TRUE, palette = "npg", lab.size = 2, 
          position = position_dodge(0.9))
ggpar(g,
      font.xtickslab = c(7,"plain", "black"),
      xtickslab.rt = 45) + labs(x="Cancer type", y="Number of differentially ranked genes")
dev.off()

#cancer based sum
summar = unique(res[,c("gene", "canc")])
canc_sum = as.data.table(table(summar$canc))
canc_sum = canc_sum[order(N)]

gene_sum = as.data.table(table(summar$gene))
gene_sum = gene_sum[order(N)]
fig1b = as.data.table(table(gene_sum$N))

#plot 1b 
#fig 1c----------------------------------------------
pdf("Manuscript_FIGURE1B.pdf", width=9)
g = ggbarplot(fig1b, "V1", "N", fill="dodgerblue4", 
              label = TRUE, lab.size = 2)
ggpar(g,
      font.xtickslab = c(7,"plain", "black")) + labs(x="# of cancer type", y="# of differentially ranked genes")
dev.off()

#summarize how many genes up/down reg in all cancers
sum_genes = as.data.table(table(res$gene, res$canc, res$sign))
sum_genes = as.data.table(filter(sum_genes, N > 0))
sum_genes2 = as.data.table(table(sum_genes$V1, sum_genes$V3))
sum_genes2 = as.data.table(filter(sum_genes2, N >=1))
sum_genes2 = sum_genes2[order(N)]
sum_genes2$combo = paste(sum_genes2$V1, sum_genes2$V2)
t = as.data.table(table(sum_genes2$combo))
t = t[order(N)]
sum_genes3 = filter(sum_genes2, N == 23)
write.csv(sum_genes3, file="genes_up_down_ALL_gtex_tcga_comparisons.csv", quote=F, row.names=F)

#save res for ss4
write.csv(res, file="all_gtex_lncRNA_comparisons_SS4.csv", quote=F, row.names=F)

