set.seed(911)

source("/u/kisaev/lncRNAs_TCGA/NewTCGA_2017scripts/Thesis_final_scripts/Prognostic_lncs_more_detail_paper/final_scripts_2019/revisions_2020/load_data.R")
#setWD
setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNAs_2019_manuscript")

#get candidates files

#Data--------------------------------------------
allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = filter(allCands, data=="TCGA") #179 unique lncRNA-cancer combos, #166 unique lncRNAs 

#which cancer types are the non-unique lncRNAs from?
allCands$Combo = NULL
allCands = allCands[,c("gene", "coef", "HR", "pval", "cancer", "gene_name")]
allCands = allCands[!duplicated(allCands), ]
cands_dups = unique(allCands$gene[which(duplicated(allCands$gene))])

#read PFI and OS results
pfi = fread("SuppTable4_PFI.txt")
colnames(pfi)[7:28] = paste("PFI", colnames(pfi)[7:28], sep="_")

os = fread("SuppTable4.txt")
pfi_cancers = c("LGG", "BRCA", "READ")
merged = merge(os, pfi, by=c("gene"))
colnames(canc_conv)[2]="cancer"
merged=merge(merged, canc_conv, by="cancer")
merged$plot = "others"
merged$plot[merged$type == "LGG"] = "LGG"
merged$plot[merged$type == "READ"] = "READ" 
merged$plot[merged$type == "BRCA"] = "BRCA" 
merged$plot = factor(merged$plot, levels=c("others", "LGG", "BRCA", "READ"))

merged$gene_plot[!(merged$plot == "others")] = merged$gene_name.y[!(merged$plot == "others")]

merged$HR  = log2(merged$HR)
merged$PFI_HR = log2(merged$PFI_HR)
z = which((merged$fdr_pval < 0.05) & (merged$PFI_fdr_pval < 0.05))
merged$significant = ""
merged$significant[z] = "sig"
merged$significant = factor(merged$significant, levels=c("sig", ""))

#1. get correlation between hazard ratios OS vs PFI and highlight correlation for LGG, READ and BRCA 

pdf("/u/kisaev/figure4X_PFI_vs_OS.pdf", width=6, height=6)
ggscatter(merged, x = "HR", y = "PFI_HR", palette=c("gray88", "#1B9E77" ,"#D95F02" ,"#7570B3"), alpha = 0.7, 
          add = "reg.line",  color="plot",  shape="significant",                             # Add regression line
          position="dodge", font.label = c(4, "plain"), 
          conf.int = FALSE,  label="gene_plot", repel = TRUE,                                 # Add confidence interval
          add.params = list(color = "dodgerblue4",
                            fill = "lightgray"), size=1
          )   +
  stat_cor(method = "spearman") + theme_bw() + xlim(c(-4,4)) +  ylim(c(-3,3)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") + 
  geom_vline(xintercept=0, linetype="dashed", color = "black") + xlab("log2HR (OS)")+ylab("log2HR (PFI)") 
dev.off()

merged = merged[,c(1:13,27,30:34, 49, 52)]
write.csv(merged, "/u/kisaev/supp_table_input_XX_PFI_wOStimes.csv", quote=F, row.names=F)


#2. get numbers of lncs that remain significant in PFI for LGG, READ and BRCA

#save supplementary table 