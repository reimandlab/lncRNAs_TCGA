source("source_code_Cox_MonteCarlo_CV_April12.R")
require(caTools)

#start with only lncRNA_intergenic
#lincs = subset(fantom, CAT_geneClass == "lncRNA_intergenic")
#z = which(colnames(rna) %in% lincs$gene)
#rna = as.data.frame(rna)
#rna = rna[,c(z, (ncol(rna)-5):ncol(rna))]

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

rna = as.data.frame(rna)

ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
#z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
#ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,6]))
ucsc <- ucsc[-z,]

#add location to lncRNAs 

#how many lncRNAs  --> 5,785

rownames(rna) = rna$patient
z1 = which(str_detect(colnames(rna), "ENSG"))	
z2 = which(colnames(rna) %in% "Cancer")
rna = rna[,c(z1,z2)]

meds_lncs = apply(rna[,1:(ncol(rna)-1)], 2, median)
meds_lncs = as.data.frame(meds_lncs)
meds_lncs$type = "lncRNA"
colnames(meds_lncs)[1] = "medianFPKM"
meds_lncs$gene = rownames(meds_lncs)
meds_lncs = merge(meds_lncs, fantom, by="gene")
meds_lncs = meds_lncs[,c(1,2,3,5)]
colnames(meds_lncs)[4] = "genetype"

#get all pcg median expression across all cancers
rownames(pcg) = pcg$patient
meds_pcgs = apply(pcg[,2:19351], 2, median)
meds_pcgs = as.data.frame(meds_pcgs)
meds_pcgs$gene = rownames(meds_pcgs)
meds_pcgs$type = "pcg"
colnames(meds_pcgs)[1] = "medianFPKM"
meds_pcgs$genetype = "pcg"
meds_pcgs = meds_pcgs[,c("gene", "medianFPKM", "type", "genetype")]

all_meds = rbind(meds_lncs, meds_pcgs)
rownames(all_meds) = 1:nrow(all_meds)

all_meds$medianFPKM = as.numeric(all_meds$medianFPKM)
#all_meds$medianFPKM = log1p(all_meds$medianFPKM)
all_meds = as.data.frame(all_meds)

logged = all_meds
logged$medianFPKM = log1p(logged$medianFPKM)

pdf("figure1supllement.pdf", width=8, height=5)
ggdensity(logged, x = "medianFPKM", color="genetype", palette = mypal)
dev.off()

#remove outliers 
z = which(all_meds$medianFPKM >= 311653)
all_meds = all_meds[-z,]

#how much bigger is the pcg median 
meds = as.data.table(aggregate(all_meds[,2], list(all_meds$type), median))	

wilcox.test(all_meds$medianFPKM[all_meds$type=="pcg"], all_meds$medianFPKM[all_meds$type=="lncRNA"])
wilcox.test(all_meds$medianFPKM[all_meds$type=="lncRNA"], all_meds$medianFPKM[all_meds$type=="pcg"])

#not logged --> remvoe outliers 
mypal = wes_palette("Royal2")

pdf("figure1supllementB_notlogged.pdf", width=6, height=7)
g = ggboxplot(all_meds, x = "type", y = "medianFPKM", color="black", fill="type", palette = mypal, notch = TRUE, size=0.5, width=0.5) +
stat_compare_means(method = "wilcox.test") + theme_minimal() + geom_jitter(position=position_jitter(width=0.05,height=0),
         alpha=0.1,
         size=0.75) + ylab("FPKM") 
print(g)
dev.off()


all_meds$logged = log1p(all_meds$medianFPKM)

pdf("figure1supllementB.pdf", width=6, height=7)
g = ggboxplot(all_meds, x = "type", y = "logged", color="black", fill="type", palette = mypal, notch = TRUE, size=0.5, width=0.5)+
stat_compare_means(method = "wilcox.test") + theme_minimal() + geom_jitter(position=position_jitter(width=0.05,height=0),
         alpha=0.1,
         size=0.75) + ylab("log1p(FPKM)") + 
  scale_y_continuous(breaks = round(seq(min(all_meds$logged), max(all_meds$logged), by = 5),1))

print(g)
a = g
dev.off()


pdf("figure1supllementB_alltypes.pdf", width=6, height=7)
g = ggboxplot(all_meds, x = "genetype", y = "logged", color="black", fill="genetype", order =c("others", "e_lncRNA", "p_lncRNA_intergenic", "p_lncRNA_divergent", "pcg"), palette = mypal, size=0.5, width=0.5, notch = TRUE)+
stat_compare_means() + theme_minimal() + geom_jitter(position=position_jitter(width=0.05,height=0),
         alpha=0.1,
         size=0.75) + ylab("log1p(FPKM)") + 
  scale_y_continuous(breaks = round(seq(min(all_meds$logged), max(all_meds$logged), by = 5),1))

g = ggpar(g,
 font.tickslab = c(7,"plain", "black"),
 xtickslab.rt = 45)

print(g)
b = g
dev.off()


#combine plots into one 
library(patchwork)
pdf("combined_boxplot.pdf", width=11)
a + b
dev.off()





