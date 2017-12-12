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
###Establish filter for survival analysis 
###---------------------------------------------------------------

