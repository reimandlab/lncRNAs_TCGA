setwd("/Users/kisaev/remote/IC_exp_cands")

library(data.table)
library(dplyr)
library(stringr)
library(metap)
library(ggpubr)
library(activePathways)

####################################################################################################
#LOAD DATA
####################################################################################################

#----RNA-Seq----------------------------------------------------------------------------------------

#SUBID results
#SUBID results using Erik's technique 
#includes results for 28 cancer types from TCGA for those with at least 50 patients 
subid = readRDS("subid_results_KI_IC_march20.rds")
subid$method = "SUBID_erik_data"
colnames(subid)[c(1,6)] = paste(colnames(subid)[c(1,6)], "Eriks_data", sep="_")
colnames(subid)[6] = "Cancer"
colnames(subid)[4] = "ensg"

#MEDIAN DICHOTOMIZED 
medians = readRDS("TCGA_ION_CHANNEL_results_March19.rds")
colnames(medians)[2:ncol(medians)] = paste(colnames(medians)[2:ncol(medians)], "median", sep="_") 
colnames(medians)[7] = "Cancer"
medians = medians[,c("gene", "HR_median", "pval_median", "fdr_pval_median", "name_median", "Cancer")]
colnames(medians)[c(1,5)] = c("ensg", "gene")

#OUTLIER DICHOTOMIZED 
outliers = readRDS("TCGA_ION_CHANNEL_results_Jan3119.rds")
colnames(outliers)[2:ncol(outliers)] = paste(colnames(outliers)[2:ncol(outliers)], "outlier", sep="_") 
colnames(outliers)[7] = "Cancer"
outliers = outliers[,c("gene", "HR_outlier", "pval_outlier", "num_risk_outlier", "perc_risk_outlier", "fdr_pval_outlier", "Cancer")]
outliers$perc_risk_outlier = as.numeric(outliers$perc_risk_outlier)
colnames(outliers)[1] = c("ensg")

#----MICROARRAY-data--------------------------------------------------------------------------------
microarray = readRDS("GBM_ics_microarray_results_feb9.rds") #GBM data

####################################################################################################
#SUMMARIZE DATA
####################################################################################################

#[1] Combine RNA-Seq data 
kidata = merge(medians, outliers, by=c("ensg", "Cancer")) #126 ion channels with both median and outlier available analyiss 
kidata_sub = merge(kidata, subid, by=c("ensg", "Cancer")) #123 ion channels with both median, outlier and subid available 
kidata_sub$HR_median = as.numeric(kidata_sub$HR_median)
kidata_sub$HR_outlier = as.numeric(kidata_sub$HR_outlier)
kidata_sub$pval_outlier = as.numeric(kidata_sub$pval_outlier)

#[2] Combine microarray or add tag because most ICs may not have affy data

#this is for GBM ONLY #####################################################
#z = which(kidata_sub$gene %in% microarray$gene)
#noaffy = kidata_sub[-z,]
#noaffy$HR_microarray = ""
#noaffy$pval_microarray = ""
#noaffy$conc_microarray = ""
#noaffy$fdr_microarray = ""
#alldat = merge(kidata_sub, microarray, by="gene")
#alldat = rbind(alldat, noaffy)

#alldat$HR_microarray = as.numeric(alldat$HR_microarray)
#alldat$pval_microarray = as.numeric(alldat$pval_microarray)
#alldat$conc_microarray = as.numeric(alldat$conc_microarray)
#alldat$fdr_microarray = as.numeric(alldat$fdr_microarray)

#round all digits
#z = which(colnames(alldat) %in% c("HR_median", "pval_median", "fdr_pval_median", "HR_outlier", "pval_outlier", "perc_risk_outlier", "fdr_pval_outlier", "median_p_Eriks_data",
#          "subid_p","subid_fdr",  "median_fdr_Eriks_data", "HR_microarray", "pval_microarray", "conc_microarray", "fdr_microarray"))

alldat = kidata_sub
alldat = as.data.frame(alldat)
#alldat[,z] = apply(alldat[,z], 2, function(x){round(x, digits=4)})
alldat = as.data.table(alldat)
alldat$method = NULL

alldat$HR_match = ""
z = which((alldat$HR_median >1) & (alldat$HR_outlier >1))
alldat$HR_match[z] = "yes"
alldat = alldat[z,]

saveRDS(alldat, file="ion_channels_3_methods_survival_analysis_KI_March20_all_cancers.rds")
write.csv(alldat, file="ion_channels_3_methods_survival_analysis_KI_March20_all_cancers.csv", quote=F, row.names=F)

alldat = as.data.table(alldat)
alldat = alldat[order(pval_outlier, subid_p, pval_median)]
alldat$rank = 1:nrow(alldat)

cands = c("GJB2", "CATSPER1", "AQP9", "SCN9A")
cands_dat = filter(alldat, gene %in% cands)
#write.csv(cands_dat, file="ion_channels_four_exo_cands_KI_Feb22.csv", quote=F, row.names=F)

####################################################################################################
#MAKE FIGURE FOR REPORT 
####################################################################################################

#keep just the cols with the pvalues
#ps = c("gene", "pval_median", "pval_outlier", "subid_p", "pval_microarray")
ps = c("gene", "pval_median", "pval_outlier", "subid_p", "Cancer")
fulldat = alldat
alldat = alldat[,..ps]
#z = which(is.na(alldat$pval_microarray))
#alldat$pval_microarray[z] = 1

alldat = as.data.frame(alldat)
alldat$combo = paste(alldat$gene, alldat$Cancer, sep="_")
rownames(alldat) = alldat$combo
alldat$combo = NULL
alldat$gene = NULL
alldat$Cancer = NULL

#get fisher's merged p-value
#alldat$fish = apply(alldat, 1, function(x){merge_p_values(x)})
#alldat$fish = round(alldat$fish, digits=7)

alldat = as.matrix(alldat)
browns = merge_p_values(alldat, method="Brown")
alldat = as.data.frame(alldat)
alldat$browns = browns

#reorder
alldat$gene = rownames(alldat)

alldat = as.data.table(alldat)
alldat = alldat[order(browns)]

#add cancer type and gene 
alldat$Cancer = sapply(alldat$gene, function(x){unlist(strsplit(x, "_"))[2]})
alldat$gene_name = sapply(alldat$gene, function(x){unlist(strsplit(x, "_"))[1]})

z = which(alldat$gene_name %in% cands)
alldat$cands = ""
alldat$cands[z] = "yes"

#write.csv(alldat, file="ion_channels_merged_pvalues_browns_KI_onlyhazardours_220219.csv", quote=F, row.names=F)

#add fdr
alldat$fdr = p.adjust(alldat$browns, method="fdr")

write.csv(alldat, file="ion_channels_merged_pvalues_browns_KI_onlyhazardours_withFDR_010319_all_cancers.csv", quote=F, row.names=F)

#make summary plot <- for GBM only 
#alldat$fdr_plot = -log10(alldat$fdr)
#ggbarplot(alldat, x = "gene", y = "fdr_plot",
#          fill = "cands",               # change fill color by cyl
#          color = "white",            # Set bar border colors to white
#          palette = "jco",            # jco journal color palett. see ?ggpar
#          sort.val = "desc",          # Sort the value in dscending order
#          sort.by.groups = FALSE,     # Don't sort inside each group
#          x.text.angle = 90           # Rotate vertically x axis texts
#) + labs(x="Ion Channel", y="-log10(FDR)")+
#  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red")

#get HRs and make heatmap 

fulldat$comb = paste(fulldat$gene_name, fulldat$Cancer, sep="_")
fulldat = fulldat[,c("HR_median", "HR_outlier", "comb")]

colnames(alldat)[5] = "comb"
alldat = merge(alldat, fulldat, by="comb")

#make heatmap
alldat = as.data.table(filter(alldat, fdr < 0.1, HR_median < 20, HR_outlier < 20))
table(alldat$Cancer)

alldat$max_hr = sapply(1:nrow(alldat), function(x){max(alldat$HR_median[x], alldat$HR_outlier[x])})
t = as.data.table(table(alldat$Cancer))
t = t[order(-N)]

k = as.data.table(table(alldat$gene_name))
k = k[order(-N)]

alldat$Cancer = factor(alldat$Cancer, levels = t$V1)
alldat$gene_name = factor(alldat$gene_name, levels = k$V1)

g = ggplot(alldat, aes(gene_name, Cancer)) + geom_tile(aes(fill=max_hr)) +
  scale_fill_gradient(low="grey", high="red", na.value = 'transparent') + labs(x = "Cancer", y="Ion channel") #+ coord_flip()

g = ggpar(g, font.xtickslab = c(3,"plain", "black"), font.ytickslab = c(5,"plain", "black"), xtickslab.rt=90)+
  theme(legend.position="none")

xplot = ggplot(alldat, aes(Cancer)) + geom_bar(fill = "black") + theme_bw()+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + coord_flip()

yplot = ggplot(alldat, aes(gene_name)) + geom_bar(fill = "black") + theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

date = Sys.Date()
pdf(paste(date, "summary_figure_ICs.pdf", sep="_"), width=12)
ggarrange(yplot, NULL, g, xplot, 
          ncol = 2, nrow = 2,  align = "hv", 
          widths = c(6, 1), heights = c(1, 4),
          common.legend = FALSE)
dev.off()

#save
colnames(alldat)[9] = "BrownsFDR"
write.csv(alldat, file="merged_pvals_all_cancers_hazardous.csv", quote=F, row.names=F)
  

