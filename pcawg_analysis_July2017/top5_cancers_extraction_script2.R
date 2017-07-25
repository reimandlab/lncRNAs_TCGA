#top5_cancers_extraction_script1.R

#Karina Isaev
#July 24th, 2017

#Purpose: Using PCAWG, extract the top 5 cancer types with the 
#most patients and good distribution of samples across
#multiple histological subtypes 

#For each cancer type, identify list of 100-200 candidate 
#lncRNAs that wiill be used for further survival and co-expression
#analysis 

#Script2 - using the top 5 cancer types chosen
#PLOTS:

#1. For each cancer type, identify list of candidate lncRNAs 
#whose expression is specific to that cancer 

#2. Plot their violin plots with increasing medians 

#3. Plot to compare how cancer specific lncRNAs are expressed 
#in other cancer types 

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

mypal = pal_npg("nrc", alpha = 0.7)(10)

#---------------------------------------------------------
#Data
#---------------------------------------------------------

#UCSC gene info
ucsc <- fread("UCSC_hg19_gene_annotations_downlJune12byKI.txt", data.table=F)
z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,6]))
ucsc <- ucsc[-z,]

#Clinical file 
clin <- fread("pcawg_specimen_histology_August2016_v6.tsv", data.table=F)
conversion <- fread("pcawgConversion.tsv", data.table=F)

#RNA-Seq file 
rna <- fread("joint_fpkm_uq.tsv", data.table=F)

#Cancers to use 
tum_types <- fread("top5_cancers_andHISTO_to_keepJuly20.txt", data.table=F)

#---------------------------------------------------------
#Processing
#---------------------------------------------------------

###"NORMAL SAMPLES"
z <- which(conversion$normal_rna_seq_aliquot_id %in% colnames(rna))
norm_pats <- conversion$icgc_donor_id[z]

###"TUMOUR SAMPLES" - PROCESSING RNA FILE 
z <- which(conversion$tumor_rna_seq_aliquot_id %in% colnames(rna))
tum_pats <- conversion$icgc_donor_id[z]

for(i in 1:ncol(rna)){
	z <- which(conversion$tumor_rna_seq_aliquot_id %in% colnames(rna)[i])
	if(!(length(z)==0)){
		colnames(rna)[i] <- conversion$icgc_donor_id[z]
	}
}

extract <- function(row){
	gene <- as.character(row[[1]])
	ens <- unlist(strsplit(gene, "::"))[3]
	return(ens)
}
rna[,1] <- apply(rna[,1:2], 1, extract) 

#seperate first by "_"
extract2 <- function(row){
	gene <- as.character(row[[1]])
	ens <- unlist(strsplit(gene, "_"))[2]
	return(ens)
}
rna[,1] <- apply(rna[,1:2], 1, extract2) 

#now only keep data for ensembl id genes 
ensg <- function(row){
	gene <- as.character(row[[1]])
	ens <- grepl("ENSG", gene)
	return(ens)
}
check <- apply(rna[,1:2], 1, ensg)
z <- which(check==TRUE)
rna <- rna[z,]

#2. Check how many IDs match the lncRNA IDs
#none match while trancript number is present 
#remove ie, ENSG00000201285.1 --> ENSG00000201285
#in both rna file and lncs file

extract3 <- function(row){
	gene <- as.character(row[[1]])
	ens <- gsub("\\..*","",gene)
	return(ens)
}
rna[,1] <- apply(rna[,1:2], 1, extract3) ; 

#Remove duplicate genes 
z <- which(duplicated(rna[,1]))
genes <- rna[z,1]
z <- which(rna[,1] %in% genes)
rna <- rna[-z,]

#Using UCSC keep only antisense, lincRNA and protein-coding genes
z <- which(rna[,1] %in% ucsc$hg19.ensGene.name2)
rna <- rna[z,]
 
#rows
rownames(rna) <- rna[,1]
rna <- rna[,-1]

##Keep only patient samples that also have clinical data
z <- which(colnames(rna) %in% clin$icgc_donor_id)
rna <- rna[,z]

#Divide RNA into lnc_RNA and pcg_RNA
z <- which(ucsc$hg19.ensemblSource.source == "protein_coding")
zz <- which(rownames(rna) %in% ucsc$hg19.ensGene.name2[z])

###################
###PCG_RNA#########
###################
pcg_rna <- rna[zz,]


###################
###LNC_RNA#########
###################
lnc_rna <- rna[-zz,]


###SUBSET CLINICAL AND EXPRESSION FILE TO ONLY THE TOP 5 CANCERS/HISTOS

##CLIN
clin_top5 <- subset(clin, clin$icgc_donor_id %in% tum_pats)
z <- which(duplicated(clin_top5$icgc_donor_id))
clin_top5 <- clin_top5[-z,]
clin_top5 <- subset(clin_top5, (clin_top5$histology_tier2 %in% tum_types$V1) & (clin_top5$histology_tier4 %in% tum_types$V2))
clin_top5$combined_tum_histo <- ""
clin_top5$combined_tum_histo <- paste(clin_top5[,15], clin_top5[,17])

#EXPRESSION - pcgs
z <- which(colnames(pcg_rna) %in% clin_top5$icgc_donor_id)
pcg_rna_top5 <- pcg_rna[,z] #20166
sums <- apply(pcg_rna_top5, 1, sum) 
s <- which(sums==0)
z <- which(rownames(pcg_rna_top5) %in% names(s))
pcg_rna_top5 <- pcg_rna_top5[-z,] #20022
pcg_rna_top5 <- t(pcg_rna_top5)
pcg_rna_top5 <- as.data.frame(pcg_rna_top5)
pcg_rna_top5$canc <- ""
#add patient tum type
for(i in 1:nrow(pcg_rna_top5)){
	pat <- rownames(pcg_rna_top5)[i]
	z <- which(clin_top5$icgc_donor_id %in% pat)
	hist <- clin_top5$combined_tum_histo[z]
	pcg_rna_top5$canc[i] <- hist
}


#EXPRESSION - lncs
z <- which(colnames(lnc_rna) %in% clin_top5$icgc_donor_id)
lnc_rna_top5 <- lnc_rna[,z] #12598
sums <- apply(lnc_rna_top5, 1, sum) 
s <- which(sums==0)
z <- which(rownames(lnc_rna_top5) %in% names(s))
lnc_rna_top5 <- lnc_rna_top5[-z,] #12543
lnc_rna_top5 <- t(lnc_rna_top5)
lnc_rna_top5 <- as.data.frame(lnc_rna_top5)
lnc_rna_top5$canc <- ""
#add patient tum type
for(i in 1:nrow(lnc_rna_top5)){
	pat <- rownames(lnc_rna_top5)[i]
	z <- which(clin_top5$icgc_donor_id %in% pat)
	hist <- clin_top5$combined_tum_histo[z]
	lnc_rna_top5$canc[i] <- hist
}


#---------------------------------------------------------
#Analysis - plot medians cutoffs
#---------------------------------------------------------

#For each subtype:
#Generate (1) plot showing how many lncRNAs left after different 
#median cutoffs --> then decide on a cutoff based on how many lncRNAs 
#there are

plots <- list()

for(i in 1:length(unique(lnc_rna_top5$canc))){
	tis <- unique(lnc_rna_top5$canc)[i]
	z <- which(lnc_rna_top5$canc %in% tis)
	tis_exp <- lnc_rna_top5[z,]
	#measure medians of genes and save them 
	meds <- as.data.frame(matrix(nrow=dim(tis_exp)[2]-1,ncol=2))
	colnames(meds) <- c("Gene", "MedianE")
	meds$Gene <- colnames(tis_exp[,1:dim(tis_exp)[2]-1])
	meds$MedianE <- apply(tis_exp[,1:dim(tis_exp)[2]-1], 2, median)
	meds$check1 <- ""
	meds$check2 <- ""
	meds$check3 <- ""
	meds$check4 <- ""
	meds$check5 <- ""
	
	for(y in 1:nrow(meds)){
		m <- meds$MedianE[y]
		if(m >=5){
			meds[y,3:7] <- 1
		}
		if(m >=4){
			meds[y,3:6] <- 1
		}
		if(m >=3){
			meds[y,3:5] <- 1
		}
		if(m >=2){
			meds[y,3:4] <- 1
		}
		if(m >=1){
			meds[y,3] <- 1
		}
		if(m < 1){
			meds[y,3] <- 0
		}
	}
	#plot how many lncRNAs meet each median cutoff
	plot_meds <- as.data.frame(matrix(nrow=5,ncol=2))
	colnames(plot_meds) <- c("Median", "Number_lncRNAs")
	plot_meds[,1] <- 1:5
	plot_meds[5,2] <- length(which(meds$check5==1))
	plot_meds[4,2] <- length(which(meds$check4==1 & (!(meds$check5==1))))
	plot_meds[3,2] <- length(which(meds$check3==1 & (!(meds$check4==1))))
	plot_meds[2,2] <- length(which(meds$check2==1 & (!(meds$check3==1))))
	plot_meds[1,2] <- length(which(meds$check1==1 & (!(meds$check2==1))))

	#save plot
	g <- ggbarplot(plot_meds, x="Median", y="Number_lncRNAs", palette=mypal, col="Median", fill="Median", label = TRUE, lab.pos = "in", lab.size = 2.6)
	g <- ggpar(g, legend="none")
	g <- g + labs(title = tis, y="Number of lncRNAs") + 
     theme(plot.title = element_text(hjust = 0.5))
    plots[[i]] <- g

    #save list of high expression lncRNAs 
    name_file <- paste(tis, "list_great5med_lncs.txt")
    save <- as.data.frame(meds[meds$MedianE >=5,1])
    colnames(save)[1] <- "gene"
    save$canc <- tis
    write.table(save, name_file, quote=F, row.names=F, sep="_")

} #end loop

g1 <- plots[[1]]
g2 <- plots[[2]]
g3 <- plots[[3]]
g4 <- plots[[4]]
g5 <- plots[[5]]
g6 <- plots[[6]]
g7 <- plots[[7]]

pdf("top5_cancers_lncRNAs_above_diffMedians.pdf", pointsize=5, height=13, width=13)
plot_grid(g1,g2,g3,g4,g5,g6,g7, labels = "AUTO", ncol = 2, align = 'v', label_size = 10, scale = 0.9)
dev.off()


#---------------------------------------------------------
#Analysis - how many lncRNAs in common?
#---------------------------------------------------------

f <- fread("all_lncs_cancers.txt", data.table=F, sep="_")
#remove individual file headers
z <- which(f[,1] %in% "gene")
f <- f[-z,]

#how many times does each gene appear? if 7 == all cancers
counts <- as.data.table(table(f[,1]))
counts <- counts[order(N)]
#173/305 lncRNAs appear only once
#21/305 lncRNAs appear in all cancer types 

all <- counts$V1[counts$N ==7]




























