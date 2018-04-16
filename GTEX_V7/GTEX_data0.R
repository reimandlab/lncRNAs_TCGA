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

#RNASeq data for 175 liver samples 
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
lncRNA = gene_reads[z,]
saveRDS(lncRNA, file="lncRNAs_GTEX_V7.rds")

pcg = gene_reads[-z,]
saveRDS(pcg, file = "PCGs_GTEX_V7.rds")

###---------------------------------------------------------------
###Process Data
###---------------------------------------------------------------

#add mean tpm to each gene and sort 
gene_reads = as.data.frame(gene_reads)
gene_reads$mean = apply(gene_reads[,3:177], 1, mean)
gene_reads$type[gene_reads$Name %in% lncRNA$Name] = "lncRNA"

pro = ucsc[ucsc$hg19.ensemblSource.source == "protein_coding",]
gene_reads$type[gene_reads$Name %in% pro$CAT_geneID] = "pcg"
gene_reads = as.data.table(gene_reads)

gene_reads = filter(gene_reads, type %in% c("lncRNA", "pcg"))
gene_reads = as.data.table(gene_reads)

gene_reads = gene_reads[order(mean)]

tcga_genes = fread("genes_used_TCGA.txt")
gene_reads = filter(gene_reads, Name %in% tcga_genes$x)
gene_reads = as.data.table(gene_reads)
genes_gtex = gene_reads$Name
write.table(genes_gtex, file = "genes_used_GTEx.txt")

###---------------------------------------------------------------
###Add scores 
###---------------------------------------------------------------

get_score = function(patient){
	z = which(colnames(gene_reads) %in% patient)
	z =  c(z, 1:2, 178:179)
	data = gene_reads[,..z]
	data = as.data.frame(data)
	data$order = data[,1]
	data = as.data.table(data)
	data = data[order(order)]
	data$rank = 1:nrow(data)
	data$score = data$rank/nrow(data)
	data$patient = patient
	data$canc = "liver_gtex"
	data = as.data.frame(data)
	data = data[,-1]
	return(data)
}

patients = as.list(colnames(gene_reads)[3:(ncol(gene_reads)-2)])
ranked_genes = llply(patients, get_score)
ranked_genes = ldply(ranked_genes, data.frame)

saveRDS(ranked_genes, file="GTEX_liver_ranked_genes_Dec14.rds")

































