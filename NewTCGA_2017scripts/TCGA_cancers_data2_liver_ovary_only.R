###---------------------------------------------------------------
###TCGA_cancers_data2.R
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
rna = readRDS("lnc_rna_ovary_liver_plus_clinical.rds")
rna = as.data.table(rna)

#2. Clinical data
clin = read.csv("all_clin_XML_tcgaSept2017.csv")
clin = clin[,1:90]
clin = as.data.table(clin)
clin = filter(clin, bcr_patient_barcode %in% rownames(rna))
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

###---------------------------------------------------------------
###Processing data 
###---------------------------------------------------------------

liver = filter(rna, canc == "Liver hepatocellular carcinoma")
ovary = filter(rna, canc == "Ovarian serous cystadenocarcinoma")

datasets = list(liver, ovary)

#1. Remove any genes with 0 counts in all people
check_low = function(df){
	sums = apply(df[,1:(ncol(df)-4)], 2, sum)
	z <- which(sums == 0)
	print(names(sums[z]))
	df = df[,-z]
	return(df)
}

datasets = llply(datasets, check_low)

#2. Plot histogram of counts/gene within each cancer type 
make_histo = function(gene){
	z <- which(colnames(df) %in% gene)
	df = df[,c((ncol(rna)-3),z)]
	gene = fantom$CAT_geneName[which(fantom$CAT_geneID %in% gene)]
	canc = df[1,1]
	df[,2] = as.numeric(df[,2])
	colnames(df)[2] = "gene"
	g1 = gghistogram(df, x = "gene", fill = "lightgray",
   add = "mean", rug = TRUE, title = paste(canc, gene, "raw counts"))
	df[,2] = log1p(df[,2])
	g2 = gghistogram(df, x = "gene", fill = "lightgray",
   add = "mean", rug = TRUE, title = paste(canc, gene, "log1p counts"))
	plot = plot_grid(g1, g2, labels = c("A", "B"), align = "h")
	print(plot)
}

pdf("Liver_distributions.pdf", width = 20, height=14)
df = liver
genes = colnames(df)[1:(ncol(df)-4)]
sapply(genes, make_histo)
dev.off()

pdf("Ovary_distributions.pdf", width = 20, height=14)
df = ovary
genes = colnames(df)[1:(ncol(df)-4)]
sapply(genes, make_histo)
dev.off()



