###---------------------------------------------------------------
###TCGA_cancers_data3.R
###---------------------------------------------------------------

###October 10th, 2017
###Goal: Get distribution of patients in each cancer type with
#RNA data 

###---------------------------------------------------------------
###Load Libraries
###---------------------------------------------------------------

source("source_file.R")

###---------------------------------------------------------------
###Load Data 
###---------------------------------------------------------------

#1. RNA data 
df_liv = readRDS("liver_lncRNA_data.rds")
dfgenes_liv = readRDS("liver_lncRNA_transposed_data_additional.rds")
dfgenes_pcg_liv = readRDS("liver_PCG_data.rds")

#2. Clinical data
clin = read.csv("all_clin_XML_tcgaSept2017.csv")
clin = clin[,1:90]
clin = as.data.table(clin)
clin = filter(clin, bcr_patient_barcode %in% rownames(df_liv))
z <- which(duplicated(clin$bcr_patient_barcode))
clin = clin[-z,]

#3. Fantom data 
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
z <- which(fantom$CAT_geneID %in% colnames(df_liv))
fantom = fantom[z,]

#4. UCSC Hg19
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,6]))
ucsc <- ucsc[-z,]

colnames(ucsc)[6] = "CAT_geneID"
fantom = merge(fantom, ucsc, by = "CAT_geneID")
fantom = as.data.table(fantom)

filter(fantom, functional_evidence == 4)

###---------------------------------------------------------------
###Convert counts to TPMs
###---------------------------------------------------------------

tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

genes <- fread("geneLengths.txt")
genes[,1] <- apply(genes[,1:2], 1, extract3)
genes = as.data.frame(genes)
rownames(genes) = genes$gene_id
genes$aggregate_length = genes$aggregate_length/1000

rownames(dfgenes_liv) = dfgenes_liv$gene
counts1 = dfgenes_liv[,1:319]
counts2 = dfgenes_pcg_liv
counts = rbind(counts1, counts2)

#make rownames the same for genes and counts dataframes
genes = genes[match(rownames(counts), rownames(genes)), ]

tpms <- apply(counts, 2, function(x) tpm(x, genes$aggregate_length))

#add mean tpm to each gene and sort 
tpms = as.data.frame(tpms)
tpms$mean = apply(tpms, 1, mean)
tpms$gene = rownames(tpms)
tpms$type[tpms$gene %in% rownames(counts1)] = "lncRNA"

pro = ucsc[ucsc$hg19.ensemblSource.source == "protein_coding",]
tpms$type[tpms$gene %in% pro$CAT_geneID] = "pcg"
tpms = as.data.table(tpms)

tpms = filter(tpms, type %in% c("lncRNA", "pcg"))
tpms = as.data.table(tpms)

tpms = tpms[order(mean)]

genes_used_TCGA = tpms$gene
write.table(genes_used_TCGA, file="genes_used_TCGA.txt", quote=F, row.names=F)

gtex_genes = fread("genes_used_GTEx.txt", header=F)
tpms = filter(tpms, gene %in% gtex_genes$V2)
tpms = as.data.table(tpms)

###---------------------------------------------------------------
###Rank all genes within each patient 
###---------------------------------------------------------------

get_score = function(patient){
	z = which(colnames(tpms) %in% patient)
	z =  c(z, 321, 322)
	data = tpms[,..z]
	data = as.data.frame(data)
	data$order = data[,1]
	data = as.data.table(data)
	data = data[order(order)]
	data$rank = 1:nrow(data)
	data$score = data$rank/nrow(data)
	data$patient = patient
	data$canc = "liver_tcga"
	data = as.data.frame(data)
	data = data[,-1]
	return(data)
}

patients = as.list(colnames(tpms)[1:(ncol(tpms)-3)])
ranked_genes = llply(patients, get_score)
ranked_genes = ldply(ranked_genes, data.frame)

saveRDS(ranked_genes, file="TCGA_liver_ranked_genes_Dec19.rds")



















