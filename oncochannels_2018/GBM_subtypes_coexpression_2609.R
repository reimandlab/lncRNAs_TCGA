library(survAUC)
#source("source_code_Cox_MonteCarlo_CV_Mar13.R")
require(caTools)
#check if this person is in my analysis: TCGA-61-2095
library(glmnet)
library(survcomp)
library(caret)
library(stringr)

#source("universal_LASSO_survival_script.R")

library(ggpubr)
library(ggrepel)
library(viridis)
library(patchwork)
library(caret)  
library(Rtsne)
library(EnvStats)
library(TCGAbiolinks)
library(stringr)

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

#rna = readRDS("rna_lncRNAs_expression_data_june29.rds")
#pcg = readRDS("rna_pcg_expression_data_june29.rds")

antibody = fread("GBM.antibody_annotation.txt")
antibody = fread("mdanderson.org_GBM.MDA_RPPA_Core.antibody_annotation.txt")

rppa = fread("GBM.rppa.txt")

#-----------------------------------------------------
#process RPPA files 
#-----------------------------------------------------

#1. change patient ids to match rna paitnet ids 
change_id = function(pat){
  l =(unlist(strsplit(pat, "-"))[1:3])
  #make sure it's tumour sample
  check= (unlist(strsplit(pat, "-"))[4])
  z = which(str_detect(check, "01"))
  if(!(length(z)==0)){
    l = paste(l[1],l[2], l[3], sep="-")
  }
  else{
    l = "no"
  }
  return(l)
}

newnames = unlist(llply(colnames(rppa)[2:ncol(rppa)], change_id))
z = which(newnames == "no")
z = z+1
rppa = as.data.frame(rppa)
rppa = rppa[,-z]
newnames = newnames[-which(newnames=="no")]
colnames(rppa)[2:ncol(rppa)] = newnames

#2. element-genes 
get_gene = function(element){
  g = unlist(strsplit(element, "\\|"))[1]
  return(g)
}
rppa$gene_id = unlist(llply(rppa[,1], get_gene))

#------FEATURES-----------------------------------------------------

#Combined into one dataframe because need to get ranks 
#all <- merge(rna, pcg, by = c("patient", "Cancer"))
#all = all[,1:25170]

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
  ) + xlab("ID1 expression") + ylab("EGFR expression")
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
          add = "jitter",ggtheme = theme_bw()) +
          stat_compare_means(ref.group = "Classical") + 
          stat_n_text()+
          ylab("ID1 expression")

        p = ggpar(p,
          font.xtickslab = c(9,"plain", "black"),
          xtickslab.rt = 65, legend="none")
        print(p)

# 1. what if the correlation is restricted to the classical subtype? 
dev.off()

pdf("ID1_EGFR_correlation_classical_subtype_spearman.pdf")
#scatter plot 
sp <- ggscatter(clin_subtypes[which(clin_subtypes$subtype == "Classical"),], x = "ENSG00000125968", y = "ENSG00000146648",
  add = "reg.line",  # Add regressin line
  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
  conf.int = TRUE # Add confidence interval
  )
# Add correlation coefficient
sp = sp + stat_cor(method = "spearman") + theme_bw() + ggtitle("GBM, ID1 vs EGFR Expression, \nClassical only")+
xlab("ID1 expression") + ylab("EGFR expression")
print(sp)

dev.off()

#can you please try sub-type weighted correlation:
#ID1~subtype vs ID1~subtype+EGFR

id1_m1 = lm(clin_subtypes$ENSG00000125968 ~ clin_subtypes$subtype)
id1_m2 = lm(clin_subtypes$ENSG00000125968 ~ clin_subtypes$subtype + clin_subtypes$ENSG00000146648)

anova(id1_m1, id1_m2)

egfr_m1 = lm(clin_subtypes$ENSG00000146648 ~ clin_subtypes$subtype)
egfr_m2 = lm(clin_subtypes$ENSG00000146648 ~ clin_subtypes$subtype + clin_subtypes$ENSG00000125968)

anova(egfr_m1, egfr_m2)

#----------------------------------------------
#####RPPA#######################################
#----------------------------------------------

# 2. what if the correlation of EGFR is apparent at proteome not transcriptome level? We could check it in the RPPA data. 
egfr_rrpa = t(rppa[which(rppa$gene_id == "EGFR"),])
#add rppa data to datafile
colnames(egfr_rrpa) = egfr_rrpa[1,]
egfr_rrpa = egfr_rrpa[-1,]
egfr_rrpa = egfr_rrpa[-(nrow(egfr_rrpa)),]
egfr_rrpa = as.data.frame(egfr_rrpa)
egfr_rrpa$patient = ""
egfr_rrpa$patient = unlist(rownames(egfr_rrpa))

clin_subtypes = merge(clin_subtypes, egfr_rrpa, by="patient")

#plot ID1 expression vs EGFR rppa levels 

#one for each EGFR antibody 
colnames(clin_subtypes)[55] = "EGFR"
colnames(clin_subtypes)[56] = "EGFR_pY1068"
colnames(clin_subtypes)[57] = "EGFR_pY1173"
clin_subtypes$EGFR = as.numeric(clin_subtypes$EGFR)
clin_subtypes$EGFR_pY1068 = as.numeric(clin_subtypes$EGFR_pY1068)
clin_subtypes$EGFR_pY1173 = as.numeric(clin_subtypes$EGFR_pY1173)


pdf("GBM_ID1_EGFR_rppa_analysis_oct1.pdf")
#scatter plot 
sp <- ggscatter(clin_subtypes, x = "ENSG00000125968", y = "EGFR",
  add = "reg.line",  # Add regressin line
  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
  conf.int = TRUE # Add confidence interval
  )
# Add correlation coefficient
sp = sp + stat_cor(method = "spearman") + theme_bw() + ggtitle("GBM, ID1 expression vs EGFR RPPA")+
xlab("ID1 gene expression")
print(sp)

sp <- ggscatter(clin_subtypes, x = "ENSG00000125968", y = "EGFR_pY1068",
  add = "reg.line",  # Add regressin line
  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
  conf.int = TRUE # Add confidence interval
  )
# Add correlation coefficient
sp = sp + stat_cor(method = "spearman") + theme_bw() + ggtitle("GBM, ID1 expression vs EGFR_pY1068 RPPA")+
xlab("ID1 gene expression")

print(sp)

sp <- ggscatter(clin_subtypes, x = "ENSG00000125968", y = "EGFR_pY1173",
  add = "reg.line",  # Add regressin line
  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
  conf.int = TRUE # Add confidence interval
  )
# Add correlation coefficient
sp = sp + stat_cor(method = "spearman") + theme_bw() + ggtitle("GBM, ID1 expression vs EGFR_pY1173 RPPA")+
xlab("ID1 gene expression")

print(sp)

dev.off()


####within classical subtype 


pdf("ID1_EGFR_rppa_correlation_classical_subtype_oct1.pdf")
#scatter plot 
sp <- ggscatter(clin_subtypes[which(clin_subtypes$subtype == "Classical"),], x = "ENSG00000125968", y = "EGFR",
  add = "reg.line",  # Add regressin line
  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
  conf.int = TRUE # Add confidence interval
  )
# Add correlation coefficient
sp = sp + stat_cor(method = "spearman") + theme_bw() + ggtitle("GBM, ID1 expression vs EGFR RPPA, \nClassical only")+
xlab("ID1 gene expression")

print(sp)

#scatter plot 
sp <- ggscatter(clin_subtypes[which(clin_subtypes$subtype == "Classical"),], x = "ENSG00000125968", y = "EGFR_pY1068",
  add = "reg.line",  # Add regressin line
  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
  conf.int = TRUE # Add confidence interval
  )
# Add correlation coefficient
sp = sp + stat_cor(method = "spearman") + theme_bw() + ggtitle("GBM, ENSG00000125968 expression vs EGFR_pY1068 RPPA, \nClassical only")+
xlab("ID1 gene expression")

print(sp)

#scatter plot 
sp <- ggscatter(clin_subtypes[which(clin_subtypes$subtype == "Classical"),], x = "ENSG00000125968", y = "EGFR_pY1173",
  add = "reg.line",  # Add regressin line
  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
  conf.int = TRUE # Add confidence interval
  )
# Add correlation coefficient
sp = sp + stat_cor(method = "spearman") + theme_bw() + ggtitle("GBM, ID1 expression vs EGFR_pY1173 RPPA, \nClassical only")+
xlab("ID1 gene expression")

print(sp)

#scatter plot 
sp <- ggscatter(clin_subtypes[which(clin_subtypes$subtype == "Classical"),], x = "ENSG00000125968", y = "EGFR_pY1173",
  add = "reg.line",  # Add regressin line
  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
  conf.int = TRUE # Add confidence interval
  )
# Add correlation coefficient
sp = sp + stat_cor(method = "spearman") + theme_bw() + ggtitle("GBM, ID1 expression vs EGFR_pY1173 RPPA, \nClassical only")+
xlab("ID1 gene expression")

print(sp)

dev.off()








































