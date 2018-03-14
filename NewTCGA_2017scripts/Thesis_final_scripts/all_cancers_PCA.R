###---------------------------------------------------------------
###Load libraries and data 
###---------------------------------------------------------------

source("source_code_Cox_MonteCarlo_CV_Mar13.R")
require(caTools)

#start with only lncRNA_intergenic
lincs = subset(fantom, CAT_geneClass == "lncRNA_intergenic")
z = which(colnames(rna) %in% lincs$gene)
rna = as.data.frame(rna)
rna = rna[,c(z, (ncol(rna)-5):ncol(rna))]

###[2.] Data splitting 

###---------------------------------------------------------------
###PCA using lncRNA expression 
#can then compare how using all genes compared to just using
#the ones chosen by LASSO at the end 
###---------------------------------------------------------------

#check if this person is in my analysis: TCGA-61-2095
library(glmnet)
library(survcomp)
library(caret)
library(stringr)
library(factoextra)

rna = subset(rna, canc %in% c("Kidney renal clear cell carcinoma", "Liver hepatocellular carcinoma", 
	"Ovarian serous cystadenocarcinoma", "Pancreatic adenocarcinoma"))

#remove 0 sums
sums = apply(rna[,1:(ncol(rna)-5)], 2, sum)
z = which(sums==0)
rna = rna[,-(which(colnames(rna) == names(z)))]

groups <- as.factor(rna$canc)
res.pca <- prcomp(rna[,1:(ncol(rna)-5)],  scale = TRUE)

pdf("4cancers_PCA_using_intergenic_lncRNAs.pdf")

# Change title and axis labels
p = fviz_pca_ind(res.pca, geom="point", label="none", habillage=rna$canc, addEllipses=TRUE, ellipse.level=0.95) +
  labs(title ="PCA using expression of all intergenic lncRNAs") + 
   xlim(-20, 20) + ylim (-30, 20)

# Change group colors using RColorBrewer color palettes
p + scale_color_brewer(palette="Dark2") +
     theme_minimal()

dev.off()

