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
library(EnvStats)
library(TCGAbiolinks)

source("check_lnc_exp_cancers.R")

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

#Combined into one dataframe because need to get ranks 
all <- merge(rna, pcg, by = c("patient", "Cancer"))
all = all[,1:25170]

#--------This script ------------------------------------------------

#-correlation of EGFR and ID1 gene expression in GBMs in TCGA. Scatter 
#plus R and P values as outcome.

#-expression levels of ID1 in three subtypes of GBM: classical, proneural, 
#mesenchymal. box plot as outcome, classical vs the two other subtypes.

#--------------------------------------------------------------------
#Clinical files - use TCGAbiolinks
#--------------------------------------------------------------------

clin_subtypes <- TCGAquery_subtype(tumor = "gbm")
#variable descirbing original subtype --> "Original.Subtype"

#subset gene expression matrix to just ID1 and EGFR 
genes = c("ENSG00000125968", "ENSG00000146648")
genes = all[,which(colnames(all) %in% c(genes, "patient"))]
genes[,2:3] = log1p(genes[,2:3])
          
clin_subtypes = merge(clin_subtypes, genes, by = "patient")
clin_subtypes$Original.Subtype = as.character(clin_subtypes$Original.Subtype)
clin_subtypes = as.data.table(clin_subtypes)
clin_subtypes = filter(clin_subtypes, Original.Subtype %in% c("Mesenchymal", "Proneural", "Classical"))

pdf("GBM_ID1_EGFR_analysis_sept26.pdf")

#scatter plot 
sp <- ggscatter(clin_subtypes, x = "ENSG00000125968", y = "ENSG00000146648",
  add = "reg.line",  # Add regressin line
  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
  conf.int = TRUE # Add confidence interval
  )
# Add correlation coefficient
sp = sp + stat_cor(method = "spearman") + theme_bw() + ggtitle("GBM, ID1 vs EGFR Expression")
print(sp)

p <- ggboxplot(clin_subtypes, x = "Original.Subtype", y = "ENSG00000125968",
          color = "Original.Subtype", order = c("Classical", "Mesenchymal", "Proneural"), 
          title = "GBM subtypes vs ID1 expression", 
          add = "jitter", ylab = "ID1 expression",  ggtheme = theme_bw()) +
          stat_compare_means(ref.group = "Classical") + 
          stat_n_text()

        p = ggpar(p,
          font.xtickslab = c(9,"plain", "black"),
          xtickslab.rt = 65, legend="none")
        print(p)

#classical vs other 2 subtypes 
clin_subtypes$subtype = ""
clin_subtypes$subtype[clin_subtypes$Original.Subtype == "Classical"] = "Classical"
clin_subtypes$subtype[!(clin_subtypes$Original.Subtype == "Classical")] = "Other"

p <- ggboxplot(clin_subtypes, x = "subtype", y = "ENSG00000125968",
          color = "subtype", order = c("Classical", "Other"), 
          title = "GBM subtypes vs ID1 expression", 
          add = "jitter", ylab = "ID1 expression",  ggtheme = theme_bw()) +
          stat_compare_means(ref.group = "Classical") + 
          stat_n_text()

        p = ggpar(p,
          font.xtickslab = c(9,"plain", "black"),
          xtickslab.rt = 65, legend="none")
        print(p)

# 1. what if the correlation is restricted to the classical subtype? 
dev.off()

pdf("ID1_EGFR_correlation_classical_subtype_pearson.pdf")
#scatter plot 
sp <- ggscatter(clin_subtypes[which(clin_subtypes$subtype == "Classical"),], x = "ENSG00000125968", y = "ENSG00000146648",
  add = "reg.line",  # Add regressin line
  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
  conf.int = TRUE # Add confidence interval
  )
# Add correlation coefficient
sp = sp + stat_cor(method = "pearson") + theme_bw() + ggtitle("GBM, ID1 vs EGFR Expression, \nClassical only")
print(sp)

dev.off()

# 2. what if the correlation of EGFR is apparent at proteome not transcriptome level? We could check it in the RPPA data. 





