#Preamble#-------------------------------------------------
options(stringsAsFactors=F)

#Later on - incorporate FANTOM5 and CRISPRi lncRNAs 

#Libraries#------------------------------------------------
library(colorout)
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
library(plyr)

mypal = pal_npg("nrc", alpha = 0.7)(10)

#Data#---------------------------------------------------

#UCSC gene info
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,6]))
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
colnames(fantom)[1] = "gene"

#lincs = subset(fantom, (CAT_geneClass == "lncRNA_intergenic") & (CAT_geneCategory %in% c("e_lncRNA", "p_lncRNA_intergenic")))
lincs = fantom

#Clinical file 
#clin <- fread("pcawg_specimen_histology_August2016_v6.tsv", data.table=F)
#conversion <- fread("pcawgConversion.tsv", data.table=F)

#RNA-Seq file 
#lncRNA
lnc_rna <- readRDS("6028_pcawg_lncRNAs_RNASeq_data.rds")
lnc_rna <- as.data.frame(lnc_rna)
lnc_rna$patient <- rownames(lnc_rna)

#PCGs
pcg_rna <- readRDS("20166_pcawg_PCGs_RNASeq_data.rds")
pcg_rna <- as.data.frame(pcg_rna)
pcg_rna$patient <- rownames(pcg_rna)

#remove duplicated column names 
dups <- colnames(pcg_rna)[which(duplicated(colnames(pcg_rna)))]   
#save them in a list for future reference 
pcg_rna <- pcg_rna[,-(which(colnames(pcg_rna) %in% dups))]

#Clinical file - available only for 485/497 patients 
clin <- readRDS("Jan26_PCAWG_clinical")
z <- which(clin$icgc_donor_id %in% rownames(lnc_rna))
clin <- clin[z,]

lnc_rna <- lnc_rna[which(rownames(lnc_rna) %in% clin$icgc_donor_id),] #485 patients remain
pcg_rna <- pcg_rna[which(rownames(pcg_rna) %in% clin$icgc_donor_id),] #485 patients remain 

#only look at lncRNAs included in fantom
z = which(colnames(lnc_rna) %in% lincs$gene)
lnc_rna = lnc_rna[,z]
#which are detectable in all cancers?
meds = apply(lnc_rna, 2, median)
det = meds[meds>=1]

#For each patient add survival status and days since last seen 
lnc_rna$canc = ""
lnc_rna$status = ""
lnc_rna$time = ""
lnc_rna$sex = ""

#lncs
for(i in 1:nrow(lnc_rna)){
  pat <- rownames(lnc_rna)[i]
  z <- which(clin$icgc_donor_id %in% pat)
  lnc_rna$canc[i] <- clin$histology_abbreviation[z]
  lnc_rna$status[i] <- clin$donor_vital_status[z]
  lnc_rna$sex[i] <- clin$donor_sex[z]
  t <- clin$donor_survival_time[z]
  if(is.na(t)){
        t <- clin$donor_interval_of_last_followup[z]
        }
        lnc_rna$time[i] <- t
}

#############################################################################
#Get list of lncRNAs detectable in each cancer-------------------------------
#############################################################################
lnc_rna = subset(lnc_rna, canc %in% c("Kidney-RCC", "Ovary-AdenoCA", "Panc-AdenoCA", "Liver-HCC"))
cancers = as.list(unique(lnc_rna$canc))

get_detectable_genes = function(cancer){
	dat = subset(lnc_rna, canc == cancer)
	meds = apply(dat[,1:(ncol(dat)-4)], 2 ,median)	
	det = meds[meds>=1]
	det = as.data.frame(det)
	det$canc = cancer
	det$gene = rownames(det)
	return(det)
}

detectable = llply(cancers, get_detectable_genes)
detectable = ldply(detectable, data.frame)
colnames(detectable)[1] = "medianFPKM"
detectable$name = ""
for(i in 1:nrow(detectable)){
	z = which(lincs$gene == detectable$gene[i])
	detectable$name[i] = lincs$CAT_geneName[z]
}

saveRDS(detectable, file="PCAWG_detectable_genes_4cancers_March20.rds")





























