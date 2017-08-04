#---------------------------------------------------------------------------
#processing_pancanatlas_files.R
#---------------------------------------------------------------------------

#Data: August 4th

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
library(factoextra)

mypal = pal_npg("nrc", alpha = 0.7)(10)

#---------------------------------------------------------------------------
#Data
#---------------------------------------------------------------------------

#Clinical file
clin <- fread("clinical_PANCANatlas_patient_with_followup.tsv")

#Expression file
exp <- fread("unc.edu_PANCAN_IlluminaHiSeq_RNASeqV2.geneExp_whitelisted.tsv")

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

#TSS codes
tss_codes <- fread("tcga_tss_codes.txt")
tss_codes <- as.data.frame(tss_codes)
colnames(tss_codes) <- c("tss", "centre", "cancer", "loc")

#---------------------------------------------------------------------------
#Processing
#---------------------------------------------------------------------------

#1. Gene IDs 

genes <- exp[,1]

#seperate first by "|"
extract2 <- function(row){
	gene <- as.character(row[[1]])
	ens <- unlist(strsplit(gene, "\\|"))[1]
	return(ens)
}
exp[,1] <- apply(exp[,1], 1, extract2) 
exp <- as.data.frame(exp)
#remove all the rows with '?'
exp <- exp[-c(1:29),]
#remove duplicated patient IDs - none 

#2. Patient IDs
clin <- as.data.frame(clin)
#from exp file need to shorten column names to this format: bcr_patient_barcode: TCGA-OR-A5J2
ids <- colnames(exp)[2:10328]

extract3 <- function(row){
	gene <- as.character(row[[1]])
	ens <- unlist(strsplit(gene, "-"))[1:4]
	p1 <- ens[1]
	p2 <- ens[2]
	p3 <- ens[3]
	p4 <- ens[4]
	ens <- paste(p1, p2, p3,p4, sep="-")
	return(ens)
}

colnames(exp)[2:10328] <-  unlist(lapply(ids, extract3))
#remove duplicated patient IDs - none 
z <- which(duplicated(colnames(exp)))
dups <- colnames(exp)[z]
z <- which(colnames(exp) %in% dups)
exp <- exp[,-z]

##Top 5 PCAWG cancers
z <- which(tss_codes[,3] %in% c("Breast invasive carcinoma", "Kidney Chromophobe", "Kidney renal clear cell carcinoma", 
	"Kidney renal papillary cell carcinoma", "Liver hepatocellular carcinoma", "Ovarian serous cystadenocarcinoma", "Pancreatic adenocarcinoma"))
tss_codes <- tss_codes[z,]

#the second segment of patient id is the TSS
patients <- as.data.frame(matrix(ncol=3, nrow=10323)) ; colnames(patients) <- c("id", "tss", "source")
patients$id <- colnames(exp)[2:10324]
for(i in 1:nrow(patients)){
	pat <- patients$id[i]
	tss <- unlist(strsplit(pat, "-"))[2]
	patients$tss[i] <- tss
}

#add cancer type
patients <- merge(patients, tss_codes, by="tss") #3094 patients remain 

#add source ie, cancer vs normal
#Tumor types range from 01 - 09, normal types from 10 - 19 and control samples from 20 - 29. See Code Tables Report for a complete list of sample codes 
for(i in 1:nrow(patients)){
	pat <- patients$id[i]
	source <- unlist(strsplit(pat, "-"))[4]
	source <-  substr(source, 1, nchar(source)-1)
	if(as.numeric(source) <= 9){
		patients$source[i] <- "tumour"
	}
	if(as.numeric(source) <=19 & as.numeric(source) >9){
    	patients$source[i] <- "normal"
	}
}

#280 normal 
#2814 tumour

#Subset Exp file to include only those tumour samples from relevant 5 PCAWG cancers 

	#1. remove duplciated gene entries 
	g <- exp[which(duplicated(exp[,1])),1] ; 
	exp <- exp[-(which(exp[,1] %in% g)),]

	rownames(exp) <- exp[,1] ; exp <- exp[,-1]
	#2. keep only genes in UCSC file and FANTOM5 file
	exp <- exp[(which(rownames(exp) %in% ucsc$hg19.ensemblToGeneName.value)), ]

#only tumour samples 
patients <- patients[patients$source=="tumour",]
z <- which(colnames(exp) %in% patients$id)
exp <- exp[,z] #2814 total patients 

saveRDS(patients, file="TOP5_PCAWGcancertypes_tcga_rnaseqfile_patients_cancertypes_conversion.rds")

#save RNA-Seq seperatley
saveRDS(exp, file="TOP5_PCAWGcancertypestcga_pancanatlas_rnaseq_processesAug3.rds")


#how many of the 50 PCAWG lncRNAs are sequenced in this file?
pcawg_lncs <- fread("ids_of_50unique_lncRNAs_pcawg.txt", data.table=F)
#convert to hugo ids 
colnames(pcawg_lncs)[1] <- "hg19.ensGene.name2"
pcawg_lncs <- merge(pcawg_lncs, ucsc, by="hg19.ensGene.name2")

#all 86 lncs, some appear more than one in multiple cancers
#how many sequences in pancan?
all_lncs <- fread("86_unique_lncs_fromPCAWG.txt", data.table=F,sep=";")
#convert to hugo ids 
colnames(all_lncs)[1] <- "hg19.ensGene.name2"
all_lncs <- merge(all_lncs, ucsc, by="hg19.ensGene.name2")
#10/86 lncRNAs 

#subset all_lncs to include only the 5 cancer specific genes
lncs <- rownames(exp)[which(rownames(exp) %in% pcawg_lncs$hg19.ensemblToGeneName.value)]
all_lncs <- subset(all_lncs, all_lncs$hg19.ensemblToGeneName.value %in% lncs)

#---------------------------------------------------------------------------
#Make file for survival analysis - for the cancer specific lncRNAs
#---------------------------------------------------------------------------

#Divide each gene into rows where each row is each patient's expression of that 
#gene plus that patients clinical info including survival status and time 

mini_exp <- exp[which(rownames(exp) %in% lncs),]

lncs <- as.data.frame(lncs)

#Make matrix, nrow= #genes(5) *#patients(2814) == 14070
#canc specific - 5 lncRNAs 

specific_genes <- as.data.frame(matrix(nrow=14070, ncol=5))
colnames(specific_genes) <- c("Gene", "Cancer", "Patient", "GeneE", "ref")

for(i in 1:nrow(lncs)){
	gene <- as.character(lncs[i,1])
	if(i == 1){
		specific_genes[1:2814,1] <- gene
		z <- which(rownames(mini_exp) %in% gene)
		dat <- as.data.frame(t(mini_exp[z,]))	
		#add cancer type to dat
		dat$canc <- ""
		for(j in 1:nrow(dat)){
			pat <- rownames(dat)[j]
			h <- which(patients$id == pat)
			dat$canc[j] <- patients$cancer[h]
		}
		specific_genes[1:2814,2] <- dat$canc
		specific_genes[1:2814,3] <- rownames(dat)
		specific_genes[1:2814,4] <- dat[,1]
		ref_cord <- which(all_lncs$hg19.ensemblToGeneName.value == gene)
		specific_genes[1:2814,5] <- all_lncs$V2[ref_cord]
	}
	if(!(i==1)){
		z <- which(is.na(specific_genes[,1]))[1]
		specific_genes[z:(z+2813), 1] <- gene
		
		z2 <- which(rownames(mini_exp) %in% gene)
		dat <- as.data.frame(t(mini_exp[z2,]))	
		#add cancer type to dat
		dat$canc <- ""
		for(j in 1:nrow(dat)){
			pat <- rownames(dat)[j]
			h <- which(patients$id == pat)
			dat$canc[j] <- patients$cancer[h]
		}
		specific_genes[z:(z+2813),2] <- dat$canc
		specific_genes[z:(z+2813),3] <- rownames(dat)
		specific_genes[z:(z+2813),4] <- dat[,1]
		ref_cord <- which(all_lncs$hg19.ensemblToGeneName.value == gene)
		specific_genes[z:(z+2813),5] <- all_lncs$V2[ref_cord]
	}
}

##Plot violin
specific_genes_logged <- specific_genes  
specific_genes_logged$GeneE <- log1p(specific_genes_logged$GeneE)

##***********************************
#Plot each cancer type seperatley ... 

for(i in 1:length(unique(specific_genes_logged$ref))){
	data <- subset(specific_genes_logged, specific_genes_logged$ref %in% unique(specific_genes_logged$ref)[i])	
	#make plot - reference group is the cancer type 
	ref <- unique(specific_genes_logged$ref)[i]

	data$Gene <- as.factor(data$Gene)
		file <- paste(ref, "cancer_specific_expressionamongOthers.pdf", sep="_")
		pdf(file, pointsize=8, width=14, height=13)
	
		for (y in seq(1, length(unique(data$Gene)), 4)) {
    
    	g <- ggviolin(data[data$Gene %in% levels(data$Gene)[y:(y+3)], ], 
                  x="Cancer", y="GeneE", color="Cancer", fill="Cancer", palette=mypal,draw_quantiles = 0.5, facet.by="Gene",
                  order = unique(data$Cancer)[c(1,3,7,2,4,5,6)])
   	    g <- ggpar(g, font.legend = c(6, "plain", "black")) 
		g <- g + labs(title = ref, y="log1p(FPKM)") + 
     	 theme(plot.title = element_text(hjust = 0.5))
		g <- g + rremove("x.text")
		g <- g + rremove("x.ticks")  
		g <- g  + geom_hline(aes(yintercept=log1p(20)), colour="#990000", linetype="dashed")
		print(g)
		}
	
	    dev.off()

}

###Using the data for the 50 lncRNA expressed in a tissue specific manner 
###save it and use it for survival analysis and co-expression analysis ... 

write.table(specific_genes_logged, file="logged_geneExpression_5unique_specific_lncRNAs_pancan.txt", quote=F, row.names=F, sep=";")

extract4 <- function(row){
	gene <- as.character(row[[3]])
	ens <- unlist(strsplit(gene, "-"))[1:4]
	p1 <- ens[1]
	p2 <- ens[2]
	p3 <- ens[3]
	ens <- paste(p1, p2, p3, sep="-")
	return(ens)
}

specific_genes_logged$Patient <- apply(specific_genes_logged[,1:3], 1, extract4)














