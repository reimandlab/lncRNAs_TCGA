#GTEX_data0.R

###---------------------------------------------------------------
###Load Libraries
###---------------------------------------------------------------

source("source_file.R")

###---------------------------------------------------------------
###Load Data 
###---------------------------------------------------------------

#1. Gene Reads
gene_reads = fread("GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct")

#2. Transcript expected count 
#transcript = fread("GTEx_Analysis_2016-01-15_v7_RSEMv1.2.22_transcript_expected_count.txt")

#3. Sample attributes 
#sample_atts = fread("GTEx_Analysis_v7_Annotations_SampleAttributesDD.xlsx")
atts = fread("GTEx_v7_Annotations_SampleAttributesDS.txt")

atts = filter(atts, SAMPID %in% colnames(gene_reads))
atts = as.data.table(atts)
#for now just need liver
#atts = filter(atts, SMTS == "Liver")

cols = colnames(gene_reads)[which(colnames(gene_reads) %in% atts$SAMPID)]
cols = c(c("Name", "Description"), cols)

#RNASeq data for all tissue types 
gene_reads = gene_reads[, ..cols]
#1. Gene IDs 
extract3 <- function(row){
	gene <- as.character(row[[1]])
	ens <- gsub("\\..*","",gene)
	return(ens)
}

gene_reads[,1] <- apply(gene_reads[,1:2], 1, extract3) ; 
#4. Subject phenotypes
#sub_pheno = fread("GTEx_Analysis_v7_Annotations_SubjectPhenotypesDD.xlsx")
pheno = fread("GTEx_v7_Annotations_SubjectPhenotypesDS.txt")

#5. lncRNAs
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
gene_reads = as.data.frame(gene_reads)
z <- which(fantom$CAT_geneID %in% gene_reads[,1])
fantom = fantom[z,]

#6. UCSC Hg19
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
#z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
#ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,6]))
ucsc <- ucsc[-z,]

colnames(ucsc)[6] = "CAT_geneID"
fantom = merge(fantom, ucsc, by = "CAT_geneID")
fantom = as.data.table(fantom)

#filter(fantom, functional_evidence == 4)
z <- which(gene_reads[,1] %in% fantom$CAT_geneID)
#lncRNA = gene_reads[z,]
#saveRDS(lncRNA, file="lncRNAs_GTEX_V7.rds")
#pcg = gene_reads[-z,]
#saveRDS(pcg, file = "PCGs_GTEX_V7.rds")

###---------------------------------------------------------------
###Process Data
###---------------------------------------------------------------

#tcga_genes = fread("genes_used_TCGA.txt")
tcga_genes = fread("all_genes_used_inRankingAnalysisTCGA_May4th.txt")
gene_reads = filter(gene_reads, Name %in% tcga_genes$x)
gene_reads = as.data.table(gene_reads)
genes_gtex = gene_reads$Name
write.table(genes_gtex, file = "genes_used_GTExApril17.txt")

tcga_genes$type = ""
tcga_genes$type[tcga_genes$x %in% fantom$CAT_geneID] = "lncRNA"
tcga_genes$type[is.na(tcga_genes$type)] = "pcg"

###---------------------------------------------------------------
###Add scores 
###---------------------------------------------------------------
rna = gene_reads
rna$Description = NULL
rna = as.data.frame(rna)
rownames(rna) = rna$Name
rna$Name = NULL
rna = t(rna)

#1. log1p
rna = log1p(rna)

#2. Get lncRNA - median within each tissue type
tissues <- unique(atts$SMTS)
tissues = tissues[c(1, 5, 6, 8, 9, 11, 12, 13, 14, 16, 26, 17, 27, 28, 30, 21)]

#Function 1
#input: tissue 
#output: list of dataframes by tissue
get_tissue_specific <- function(tissue){
	z = which(atts$SMTS == tissue)
	pats = atts$SAMPID[z]
	z = which(rownames(rna) %in% pats)
	tis = rna[z,]
	tis = as.data.frame(tis)
	tis$tis = tissue
	return(tis)
}
tissues_data <- lapply(tissues, get_tissue_specific)

#Function 2
#input: dataframe with lncRNA/pcg-RNAseq data 
#output: new row added to dataframe indicating gene's score within 
#each patient 

addScores <- function(dtt){
	names <- as.list(rownames(dtt))
	
	getScores <- function(patient){
	z = which(rownames(dtt) %in% patient)
	row = dtt[z,]
	score=""
	expression <- data.frame(exp=as.numeric(row[1:(length(row)-1)]), gene=names(row)[1:(length(row)-1)]) #check if characters 
	expression$score <- score
	
	expression <- as.data.table(expression)
	expression <- expression[order(exp)]
	expression$score <- as.numeric(rownames(expression))/length(rownames(expression))
	
	#subset to just lncrnas
	lncs = tcga_genes$x[tcga_genes$type == "lncRNA"]
	z <- which(expression$gene %in% lncs)
	expression <- expression[z, ]
	expression = as.data.frame(expression)
	expression$patient = ""
	expression$patient = patient 
	return(expression)
		}

	patients <- llply(names, getScores, .progress="text") #list of dataframes, need to coerce together
	#names <- rownames(dataframe)
	patients1 <- rbindlist(patients)
	patients1 <- as.data.frame(patients1)
	patients1$data <- "GTEX"
	patients1$tis = dtt$tis[1]
	return(patients1)
}	

scored <- llply(tissues_data, addScores, .progress="text") #list of dataframes
all_tissues_scored <-  rbindlist(scored)
all_tissues_scored$exp = as.numeric(all_tissues_scored$exp)

#saveRDS(all_tissues_scored, "allGTEX_lncRNAs_scored_May23.rds")
saveRDS(all_tissues_scored, "allGTEX_lncRNAs_scored_Aug21.rds")































