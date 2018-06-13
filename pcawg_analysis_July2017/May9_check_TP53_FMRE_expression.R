
library(data.table)
library(ggplot2)
library(ggpubr)
library(stringr)
library(dplyr)
library(gridExtra)
library(grid)

pcg_rna = readRDS(file="all_rna_may8th.rds")
#subset to TP53 
z = which(colnames(pcg_rna) %in% c("ENSG00000141510", "patient", "canc"))
pcg_rna = pcg_rna[,z]

#[5] all patients in cohort
load("patient2cancertype.rsav")
head(patient2cancer_type) #1844 all together 

fmre_muts = readRDS(file="TP53_FMRE_patients_that_haveit.rds")

#subset RNA-seq to only people in the mutation file
z = which(rownames(pcg_rna) %in% names(patient2cancer_type))
pcg_rna = pcg_rna[z,] #831 

pcg_rna$fmre = ""
z = which(pcg_rna$patient %in% fmre_muts$patient)
pcg_rna$fmre[z] = "FMRE"
pcg_rna$fmre[-z] = "noFMRE"

pdf("TP53_exp_FMRE_vs_no_all_cancers.pdf")
p <- ggboxplot(pcg_rna, x = "fmre", y = "ENSG00000141510",
         color = "fmre", 
         palette = "jco", title = paste("All Cancers - TP53 exp ~ FMRE", table(pcg_rna$fmre)[1], "muts", table(pcg_rna$fmre)[2], "nomuts", sep=" "), 
          add = "jitter")
# Change method
p + stat_compare_means(method = "wilcox.test")

#only pancreatic**************************************************************************
panc = names(patient2cancer_type)[patient2cancer_type == "Panc-AdenoCA"]
#panc patients that have mutation 
fmre_muts[fmre_muts$patient %in% panc,]
panc =subset(pcg_rna, pcg_rna$patient %in% panc)

#pdf("TP53_exp_FMRE_vs_no_pancreatic_cancer_only.pdf")
p <- ggboxplot(panc, x = "fmre", y = "ENSG00000141510",
	color = "fmre",
         palette = "jco", title = paste("Pancreatic Cancer - TP53 exp ~ FMRE", table(panc$fmre)[1], "muts", table(panc$fmre)[2], "nomuts", sep=" "), 
          add = "jitter")
# Change method
p + stat_compare_means(method = "wilcox.test")
#p + stat_compare_means()

#only ovarian**************************************************************************
ov = names(patient2cancer_type)[patient2cancer_type == "Ovary-AdenoCA"]
ov =subset(pcg_rna, pcg_rna$patient %in% ov)

p <- ggboxplot(ov, x = "fmre", y = "ENSG00000141510",
	color = "fmre",
         palette = "jco", title = paste("Ovarian Cancer - TP53 exp ~ FMRE", table(ov$fmre)[1], "muts", table(ov$fmre)[2], "nomuts", sep=" "), 
          add = "jitter")
# Change method
p + stat_compare_means(method = "wilcox.test")

#only uterus**************************************************************************
uterus = names(patient2cancer_type)[patient2cancer_type == "Uterus-AdenoCA"]
uterus =subset(pcg_rna, pcg_rna$patient %in% uterus)

p <- ggboxplot(uterus, x = "fmre", y = "ENSG00000141510",
	color = "fmre",
         palette = "jco", title = paste("Uterus Cancer - TP53 exp ~ FMRE", table(uterus$fmre)[1], "muts", table(uterus$fmre)[2], "nomuts", sep=" "), 
          add = "jitter")
# Change method
p + stat_compare_means(method = "wilcox.test")

#only stomach**************************************************************************
stomach = names(patient2cancer_type)[patient2cancer_type == "Stomach-AdenoCA"]
stomach =subset(pcg_rna, pcg_rna$patient %in% stomach)

p <- ggboxplot(stomach, x = "fmre", y = "ENSG00000141510",
	color = "fmre",
         palette = "jco", title = paste("Stomach Cancer - TP53 exp ~ FMRE", table(stomach$fmre)[1], "muts", table(stomach$fmre)[2], "nomuts", sep=" "), 
          add = "jitter")
# Change method
p + stat_compare_means(method = "wilcox.test")

#only lung**************************************************************************
lungSCC = names(patient2cancer_type)[patient2cancer_type == "Lung-SCC"]
lungSCC =subset(pcg_rna, pcg_rna$patient %in% lungSCC)

p <- ggboxplot(lungSCC, x = "fmre", y = "ENSG00000141510",
	color = "fmre",
         palette = "jco", title = paste("LungSCC Cancer - TP53 exp ~ FMRE", table(lungSCC$fmre)[1], "muts", table(lungSCC$fmre)[2], "nomuts", sep=" "), 
          add = "jitter")
# Change method
p + stat_compare_means(method = "wilcox.test")

#only breast**************************************************************************
breast = names(patient2cancer_type)[patient2cancer_type == "Breast-AdenoCa"]
breast =subset(pcg_rna, pcg_rna$patient %in% breast)

p <- ggboxplot(breast, x = "fmre", y = "ENSG00000141510",
	color = "fmre",
         palette = "jco", title = paste("Lung Cancer - TP53 exp ~ FMRE", table(breast$fmre)[1], "muts", table(breast$fmre)[2], "nomuts", sep=" "), 
          add = "jitter")
# Change method
p + stat_compare_means(method = "wilcox.test")

#only lungadeno**************************************************************************
lungAdeno = names(patient2cancer_type)[patient2cancer_type == "Lung-AdenoCA"]
lungAdeno =subset(pcg_rna, pcg_rna$patient %in% lungAdeno)

p <- ggboxplot(lungAdeno, x = "fmre", y = "ENSG00000141510",
	color = "fmre",
         palette = "jco", title = paste("LungAdeno Cancer - TP53 exp ~ FMRE", table(lungAdeno$fmre)[1], "muts", table(lungAdeno$fmre)[2], "nomuts", sep=" "), 
          add = "jitter")
# Change method
p + stat_compare_means(method = "wilcox.test")

#only headSCC**************************************************************************
headSCC = names(patient2cancer_type)[patient2cancer_type == "Head-SCC"]
headSCC =subset(pcg_rna, pcg_rna$patient %in% headSCC)

p <- ggboxplot(headSCC, x = "fmre", y = "ENSG00000141510",
	color = "fmre",
         palette = "jco", title = paste("HeadSCC Cancer - TP53 exp ~ FMRE", table(headSCC$fmre)[1], "muts", table(headSCC$fmre)[2], "nomuts", sep=" "), 
          add = "jitter")
# Change method
p + stat_compare_means(method = "wilcox.test")

#only Colorectal**************************************************************************
colorectal = names(patient2cancer_type)[patient2cancer_type == "ColoRect-AdenoCA"]
colorectal =subset(pcg_rna, pcg_rna$patient %in% colorectal)

p <- ggboxplot(colorectal, x = "fmre", y = "ENSG00000141510",
	color = "fmre",
         palette = "jco", title = paste("Colorectal Cancer - TP53 exp ~ FMRE", table(colorectal$fmre)[1], "muts", table(colorectal$fmre)[2], "nomuts", sep=" "), 
          add = "jitter")
# Change method
p + stat_compare_means(method = "wilcox.test")

#only bladder**************************************************************************
bladder = names(patient2cancer_type)[patient2cancer_type == "Bladder-TCC"]
bladder =subset(pcg_rna, pcg_rna$patient %in% bladder)

p <- ggboxplot(bladder, x = "fmre", y = "ENSG00000141510",
	color = "fmre",
         palette = "jco", title = paste("Bladder Cancer - TP53 exp ~ FMRE", table(bladder$fmre)[1], "muts", table(bladder$fmre)[2], "nomuts", sep=" "), 
          add = "jitter")
# Change method
p + stat_compare_means(method = "wilcox.test")

#only prosate*************************************************************************
prostate = names(patient2cancer_type)[patient2cancer_type == "Prost-AdenoCA"]
prostate =subset(pcg_rna, pcg_rna$patient %in% prostate)

p <- ggboxplot(prostate, x = "fmre", y = "ENSG00000141510",
  color = "fmre",
         palette = "jco", title = paste("Prostate Cancer - TP53 exp ~ FMRE", table(prostate$fmre)[1], "muts", table(prostate$fmre)[2], "nomuts", sep=" "), 
          add = "jitter")
# Change method
p + stat_compare_means(method = "wilcox.test")

dev.off()












