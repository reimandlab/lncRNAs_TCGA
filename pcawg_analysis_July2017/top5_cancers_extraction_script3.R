#top5_cancers_extraction_script3.R

#Karina Isaev
#July 27th, 2017

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
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,6]))
ucsc <- ucsc[-z,] #total unique genes left 32843

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
fantom[,1] <- apply(fantom[,1:2], 1, extract3)


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

#keep only FANTOM lncRNAs
z <- which(colnames(lnc_rna_top5) %in% fantom$CAT_geneID)
lnc_rna_top5 <- lnc_rna_top5[,c(z,12544)]

#change to hugo ids 
for(i in 1:5607){
	g <- colnames(lnc_rna_top5)[i]
	colnames(lnc_rna_top5)[i] <- fantom$CAT_geneName[which(fantom$CAT_geneID %in% g)]
}

saveRDS(lnc_rna_top5, "5607_pcawg_lncRNAs_RNASeq_data.rds")

#change PCGs to Hugo Ids
for(i in 1:20166){
	g <- colnames(pcg_rna_top5)[i]
	colnames(pcg_rna_top5)[i] <- ucsc$hg19.ensemblToGeneName.value[which(ucsc$hg19.ensGene.name2 %in% g)]
}

saveRDS(pcg_rna_top5, "20166_pcawg_PCGs_RNASeq_data.rds")

#---------------------------------------------------------
#Analysis - how many/which lncRNAs are specific to each
#cancer type?
#---------------------------------------------------------

f <- fread("all_lncs_cancers.txt", data.table=F, sep="_")
#remove individual file headers
z <- which(f[,1] %in% "gene")
f <- f[-z,]

#subset to include only lncs covered by FANTOM
z <- which(f[,1] %in% fantom[,1])
f <- f[z,]

#remove 7SK gene from the list as its not cancer unique 
z <- which(f[,1] %in% fantom[which(fantom$CAT_geneName=="7SK"),1])
f <- f[-z,] #end up with 86 unique genes in the list

#how many times does each gene appear? if 7 == all cancers
counts <- as.data.table(table(f[,1], f[,2]))
counts <- counts[order(N)]

#remove 0s 
counts <- counts[!(counts$N==0)]

#want to get list of genes that occur in only one cancer type what are they?
once <- as.data.table(table(f[,1]))
once <- once[order(N)]
once <- once[once$N==1] #50 unique genes that appear in one cancer 

counts <- counts[which(counts$V1 %in% once$V1), ]
counts <- as.data.frame(counts)

z <- which(ucsc$hg19.ensGene.name2 %in% counts[,1])
#counts$gene <- ""
#counts$gene <- ucsc$hg19.ensemblToGeneName.value[z]

num_cancers <-  as.data.table(table(counts$V2))
colnames(num_cancers)[1] <- "new"

tum_types$new <- ""
for(i in 1:nrow(tum_types)){
	t1 <- tum_types$V1[i]
	t2 <- tum_types$V2[i]
	n <- paste(t1, t2, sep=" ")
	tum_types$new[i] <- n
}

num_cancers <- merge(tum_types, num_cancers, by="new")
num_cancers <- num_cancers[,-(2:3)]
num_cancers$rows <- ""
for(i in 1:nrow(num_cancers)){
	t1 <- num_cancers$V3[i]
	t2 <- num_cancers$N[i]
	n <- t1*t2
	num_cancers$rows[i] <- n
}

counts$num_patients <- ""
for(i in 1:nrow(counts)){
	t <- counts[i,2]
	z <- which(num_cancers$new %in% t)
	counts$num_patients[i] <- num_cancers$V3[z]
}

rows <- sum(as.numeric(num_cancers$rows))
to_plot <- as.data.frame(matrix(ncol=4, nrow=rows))
colnames(to_plot) <- c("Gene", "Cancer", "Patient", "GeneE")

for(i in 1:nrow(counts)){
	gene2 <- counts[i,1]
	tis <- counts$V2[i]
	if(i == 1){
		a <- as.numeric(counts$num_patients[i])
		to_plot[1:a,1] <- gene2
		z <- which(colnames(lnc_rna_top5) %in% gene2)
		dat <- lnc_rna_top5[,c(z,12544)]
		z <- which(dat$canc %in% tis)	
		dat <- dat[z,]
		to_plot[1:a,2] <- dat$canc
		to_plot[1:a,3] <- rownames(dat)
		to_plot[1:a,4] <- dat[,1]
	}
	if(!(i==1)){
		a <- as.numeric(counts$num_patients[i])-1
		z <- which(is.na(to_plot[,1]))[1]
		to_plot[z:(z+a), 1] <- gene2

		z2 <- which(colnames(lnc_rna_top5) %in% gene2)
		dat <- lnc_rna_top5[,c(z2,12544)]	
		
		z3 <- which(dat$canc %in% tis)	
		dat <- dat[z3,]

		to_plot[z:(z+a),2] <- dat$canc
		to_plot[z:(z+a),3] <- rownames(dat)
		to_plot[z:(z+a),4] <- dat[,1]
	}
}

to_plot$Cancer <- as.factor(to_plot$Cancer)

##Plot violin
to_plot$gene2 <- ""
for(i in 1:nrow(to_plot)){
	g <- to_plot$Gene[i]	 
	z <- which(fantom$CAT_geneID %in% g)
	to_plot$gene2[i] <- fantom$CAT_geneName[z]
}
to_plot_logged <- to_plot  
to_plot_logged$GeneE <- log1p(to_plot_logged$GeneE)


pdf("50_high_SPECIFIC_bycancer_tissues_lncs.pdf", pointsize=5, width=16, height=14)
g <- ggviolin(to_plot_logged, x="gene2", y="GeneE", color="Cancer", fill="Cancer", palette=mypal,draw_quantiles = 0.5)
g <- ggpar(g, font.legend = c(8, "plain", "black"), x.text.angle=75, font.x=c(6, "plain", "black")) 
g <- g + labs(title = "50 lncRNAs Highly Expressed in Individual Cancers", y="log1p(FPKM)", x="Gene") + 
     theme(plot.title = element_text(hjust = 0.5))
   g <- g  + geom_hline(aes(yintercept=log1p(10)), colour="#990000", linetype="dashed")
g
dev.off()




















