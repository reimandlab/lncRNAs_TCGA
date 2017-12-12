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
rna$patient = rownames(rna)
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
z <- which(fantom$CAT_geneID %in% colnames(rna))
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

z <- which(colnames(rna) %in% fantom$CAT_geneID)
rna = as.data.frame(rna)
rna = rna[,c(z, 5920:5924)]
rna = as.data.table(rna)

###---------------------------------------------------------------
###Processing data 
###---------------------------------------------------------------

liver = filter(rna, canc == "Liver hepatocellular carcinoma")
ovary = filter(rna, canc == "Ovarian serous cystadenocarcinoma")

datasets = list(liver, ovary)

#1. Remove any genes with 0 counts in all people
check_low = function(df){
	sums = apply(df[,1:(ncol(df)-5)], 2, sum)
	z <- which(sums == 0)
	print(names(sums[z]))
	df = df[,-z]
	rownames(df) = df$patient
	return(df)
}

datasets = llply(datasets, check_low)
df = datasets[[1]] #liver
dfgenes = t(df[,1:(ncol(df)-5)])
dfgenes = as.data.frame(dfgenes)
dfgenes$mean = apply(dfgenes, 1, mean)
dfgenes$gene = rownames(dfgenes)
dfgenes = as.data.table(dfgenes)
dfgenes = dfgenes[order(mean)]

dfgenes = as.data.frame(dfgenes)
dfgenes$length = ""
dfgenes$score = ""

fantom = as.data.frame(fantom)
fantom$length = fantom$hg19.ensGene.txEnd - fantom$hg19.ensGene.txStart

for(i in 1:nrow(dfgenes)){
	z <- which(fantom$CAT_geneID %in% dfgenes$gene[i])
	dfgenes$length[i] = fantom$length[z]
	dfgenes$score[i] = fantom$functional_evidence[z]
}

#2. How does expression differ based on length of lncRNA? 
#ie - is expression associated with length? 

#Summary stats 
df = subset(dfgenes, select = c("mean", "length", "gene", "score"))
df$length = as.numeric(df$length)
df$score = as.numeric(df$score)
summary(df)
superhigh = df$gene[which(df$mean > 5000)]
z <- which(df$mean > 5000)
df = df[-z,]
pdf("length_lncRNA_versusNumReads.pdf")
plot(df[,1:2])
dev.off()

#more zoomed in 
z <- which(df$mean > 1000)
df = df[-z,]
pdf("length_lncRNA_versusNumReadsZoomed.pdf")
plot(df[,1:2])
dev.off()

#doesn't seem to be a correlation there
# Fit our regression model
lengthmod <- lm(mean ~ length + score, # regression formula
              data=df) # data set
# Summarize and print the results
summary(lengthmod) # show regression coefficients table
df = datasets[[1]] #liver
saveRDS(df, file="liver_lncRNA_data.rds")
saveRDS(dfgenes, file="liver_lncRNA_transposed_data_additional.rds")


#3. Plot histogram of counts/gene within each cancer type 
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
df = datasets[[1]]
genes = colnames(df)[1:(ncol(df)-4)]
sapply(genes, make_histo)
dev.off()

#pdf("Ovary_distributions.pdf", width = 20, height=14)
#df = ovary
#genes = colnames(df)[1:(ncol(df)-4)]
#sapply(genes, make_histo)
#dev.off()




