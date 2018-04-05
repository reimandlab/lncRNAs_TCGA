#Preamble#-------------------------------------------------
options(stringsAsFactors=F)

#Later on - incorporate FANTOM5 and CRISPRi lncRNAs 

#Libraries#------------------------------------------------
library(data.table)
library(survival)
library(ggplot2)
library(ggfortify)
library(cluster)
library(qqman)
library(dplyr)
library(rafalib)
library(RColorBrewer) 
library(genefilter)
library(gplots) ##Available from CRAN
library(survminer)
library(MASS)
library(Hmisc)
library(gProfileR)
library(wesanderson)
library(ggsci)
library(gridExtra)
library(ggpubr)
library(ggthemes)
library(plyr)

mypal = pal_npg("nrc", alpha = 0.7)(10)

#Data-------------------------------------------------------

#List of canddidates and cox results
#allCands <- fread("7tier1_35tier2_lncRNA_candidates_August28th.txt", sep=";")
#List of canddidates and cox results
allCands <- readRDS("chosen_features_wFANTOM_data_Mar22_1000CVs_8020splits.rds")

#UCSC gene info
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
ucsc <- ucsc[z,]
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

#Clinical file 
#clin <- fread("pcawg_specimen_histology_August2016_v6.tsv", data.table=F)
#conversion <- fread("pcawgConversion.tsv", data.table=F)

#RNA-Seq file 
#lncRNA
lnc_rna <- readRDS("5919_lncRNAs_tcga_all_cancers_March13_wclinical_data.rds")
lnc_rna <- as.data.frame(lnc_rna)
lnc_rna$patient <- rownames(lnc_rna)

#PCGs
pcg_rna <- readRDS("19438_lncRNAs_tcga_all_cancers_March13_wclinical_data.rds")
pcg_rna <- as.data.frame(pcg_rna)
pcg_rna$patient <- rownames(pcg_rna)

#remove duplicated column names 
dups <- colnames(pcg_rna)[which(duplicated(colnames(pcg_rna)))]   

#Combined into one dataframe because need to get ranks 
all <- merge(lnc_rna, pcg_rna, by = c("canc", "time", "status", "sex", "patient"))

genes <- fread("all_genes_used_inRankingAnalysisPCAWG_Mar26.txt", sep=";")
z = which(colnames(all) %in% genes$x)
all = all[,c(1:5, z)]
all = all[,-c(6:10)]

#####-----------------------------------------------------------------------------------------------
####Plot distribution plot comparing to PCGs--------------------------------------------------------
#####-----------------------------------------------------------------------------------------------

#1. remove all genes that are not expressed at all
lnc_rna = subset(lnc_rna, canc %in% c("Ovarian serous cystadenocarcinoma", "Liver hepatocellular carcinoma", 
	"Kidney renal clear cell carcinoma", "Pancreatic adenocarcinoma"))
pcg_rna = pcg_rna[which(rownames(pcg_rna) %in% rownames(lnc_rna)),]

sums = apply(lnc_rna[,1:(ncol(lnc_rna)-5)], 2, sum)
z = which(sums==0)
lnc_rna = lnc_rna[,-(which(colnames(lnc_rna) %in% names(z)))]

sums = apply(pcg_rna[,1:(ncol(pcg_rna)-5)], 2, sum)
z = which(sums==0)
pcg_rna = pcg_rna[,-(which(colnames(pcg_rna) %in% names(z)))]


####################################################
#2. get median FPKM for each gene 
#####PLOTTING#######################################

pdf("dist_lncRNA_medians_TCGA_ind_cancers.pdf", width=5, height=5)
for(i in 1:length(unique(lnc_rna$canc))){
canc = subset(lnc_rna, canc == unique(lnc_rna$canc)[i])
meds_lncs = apply(canc[,1:(ncol(canc)-5)], 2, median)
meds_lncs = as.data.frame(meds_lncs)
meds_lncs$gene = rownames(meds_lncs)
meds_lncs$type = ""
for(i in 1:nrow(meds_lncs)){
	z = which(fantom$CAT_geneID == meds_lncs$gene[i])
	meds_lncs$type[i] = fantom$CAT_geneCategory[z]
}
colnames(meds_lncs)[1] = "median"
meds_lncs = as.data.table(meds_lncs)

pcg_canc = pcg_rna[which(rownames(pcg_rna) %in% rownames(canc)),]
meds_pcgs = apply(pcg_canc[,1:(ncol(pcg_canc)-5)], 2, median)
meds_pcgs = as.data.frame(meds_pcgs)
meds_pcgs$gene = rownames(meds_pcgs)
meds_pcgs$type = "pcg"
colnames(meds_pcgs)[1] = "median"
meds_pcgs = as.data.table(meds_pcgs)

all_meds = rbind(meds_lncs, meds_pcgs)
#all_meds$median = floor(all_meds$median)
all_meds$type = as.factor(all_meds$type)

#remove super highly expressed outliers
z = which(all_meds$median >= 188860)
if(!(length(z)==0)){
all_meds = all_meds[-z,]
}

#zoomed in 
all_meds$gene = NULL
gg <- ggplot(all_meds) + labs(title=canc$canc[1], y = "density", x="median FPKM")
gg <- gg + geom_density(aes(x=median, y=..scaled.., colour=type, fill=type), alpha=0.01)
gg <- gg + theme_bw()
print(ggpar(gg, xlim=c(0,100000)))

#zoomed out
gg <- ggplot(all_meds) + labs(title=canc$canc[1], y = "density", x="median FPKM")
gg <- gg + geom_density(aes(x=median, y=..scaled.., colour=type, fill=type), alpha=0.01)
gg <- gg + theme_bw()
#ggpar(p, xlim=c(0,250))
print(gg)

#3. log medians first 
all_meds$median = log1p(all_meds$median)
#zoomed in 
gg <- ggplot(all_meds) + labs(title=canc$canc[1], y = "density", x="log1p(median FPKM)")
gg <- gg + geom_density(aes(x=median, y=..scaled.., colour=type, fill=type), alpha=0.01)
gg <- gg + theme_bw()
print(gg)
#print(ggpar(gg, xlim=c(0,100)))
}
dev.off()
