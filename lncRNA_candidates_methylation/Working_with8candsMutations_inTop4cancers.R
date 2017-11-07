###Working_with8candsMutations_inTop4cancers.R
#+++++++++++++++++++++++++++++++++++++++++++++

#Purpose: this script will iteriterate over the mini
#methylation files to extract probes for 
#ovarian cancer patients, overlapping lncRNA regions 

###Load Libraries
#+++++++++++++++++++++++++++++++++++++++++++++

source("source_file.R")
library(GenomicRanges)
library(stringr)

#---------------------------------------------------------
#Data
#---------------------------------------------------------

#UCSC gene info
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
ucsc <- ucsc[z,]
#z <- which(duplicated(ucsc[,6]))
#ucsc <- ucsc[-z,]

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

#Candidate lncRNAs 
lncs = fread("results_October12_42candsFromPCAWG.txt")
lncs = filter(lncs, pval < 0.05)

#ucsc only keep those lncrnas in fantom
lincs = ucsc[which(ucsc$hg19.ensGene.name2 %in% lncs$gene),]
colnames(lincs)[6] = "gene"
lncs = merge(lncs, lincs, by = "gene")
lncs = lncs[,c(1:6, 12:13)]
colnames(lncs)[8] = "hugo"

#patients
canc_pats = fread("485_patient_IDs_top5CancersPCAWG.txt")

#Intersected mutations with 8 candidates 
muts = fread("muts_2500bpUpstream_200downstream_8candsLNCRNAs.bed")
colnames(muts) = c("lnc_chr", "lnc_statrt", "lnc_end", "lnc_ensg", "lnc_gene", "type", "transcript", "mut_chr", "mut_start", "mut_end", 
	"patient", "canc")

#people that went through mutation analysis 

#Right now if gene is affected by mutation, each trancript will be recorded, ie same gene mutation repeated
#multiple times in file
mut_ppl = fread("mutated_samples.txt", header=F)

#---------------------------------------------------------
#How many patients does each mutation affect? which cancers?
#predicted cancer?
#---------------------------------------------------------

genes = unique(muts$lnc_gene)
get_mut_info = function(gene){
	gene_mut = filter(muts, lnc_gene == gene)
	z <- which(duplicated(gene_mut$patient))
	if(!(length(z)==0)){gene_mut = gene_mut[-z,]}
	#which is the predicted cancer to be mutated 
	z <- which(lncs$hugo == gene)
	canc = lncs$canc[z[[1]]]
	gene_mut$canc = str_sub(gene_mut$canc, 1, 4)
	#mutation frequency
	cancs = as.data.table(table(gene_mut$canc))
	cancs = cancs[order(-N)]
	result_list = list(gene, canc, cancs, gene_mut)
	return(result_list)
}

genes_muts = llply(genes, get_mut_info)


#---------------------------------------------------------
#For lncRNAs with mutations in resepctive cancer 
#Is there a difference in expression between high and low
#lncRNA expressing patients? 
#---------------------------------------------------------

genes = genes[c(4, 8)]
lnc_rna <- readRDS("5607_pcawg_lncRNAs_RNASeq_data.rds")
lnc_rna <- as.data.frame(lnc_rna)
lnc_rna$patient <- rownames(lnc_rna)
z <- which(lnc_rna$patient %in% canc_pats$patient)
lnc_rna = lnc_rna[z,]

#canc = c("Kidney Adenocarcinoma, clear cell type", "Ovary Serous cystadenocarcinoma", "Liver Hepatocellular carcinoma", "Pancreas Pancreatic ductal carcinoma")
z <- which(lnc_rna$canc %in% "Liver Hepatocellular carcinoma")
lnc_rna = lnc_rna[z,]

#Keep only patients that have gone through mutation testing 
z <- which(colnames(lnc_rna) %in% genes)
lnc_rna = lnc_rna[,c(z, 5608:5609)]
z <- which(mut_ppl$V1 %in% lnc_rna$patient)
lnc_rna$Neat1mut = ""
neat1 = genes_muts[[4]][[4]]
neat1pats = neat1$patient[which(neat1$patient %in% lnc_rna$patient)]
lnc_rna$Neat1mut[which(lnc_rna$patient %in% neat1pats)] = "Mut"
lnc_rna$Neat1mut[lnc_rna$Neat1mut == ""] = "NotMut"

lnc_rna$Adormut = ""
ador = genes_muts[[8]][[4]]
adorpatients = ador$patient[which(ador$patient %in% lnc_rna$patient)]
lnc_rna$Adormut[which(lnc_rna$patient %in% adorpatients)] = "Mut"
lnc_rna$Adormut[lnc_rna$Adormut == ""] = "NotMut"

colnames(lnc_rna)[1] = "ADORA2A"

#Add high or low expression tag 
med_neat1 = median(as.numeric(lnc_rna$NEAT1))
med_ador = median(as.numeric(lnc_rna$ADORA2A))

lnc_rna$neat1TAG[lnc_rna$NEAT1 >= med_neat1] = "High"
lnc_rna$neat1TAG[lnc_rna$NEAT1 < med_neat1] = "Low"

lnc_rna$adorTAG[lnc_rna$ADORA2A >= med_ador] = "High"
lnc_rna$adorTAG[lnc_rna$ADORA2A < med_ador] = "Low"

#Compare expression 
pdf("neat1MutationvsExpression_84Liver.pdf", width=9)
g  = ggboxplot(lnc_rna, x="Neat1mut", y="NEAT1", add="jitter", palette=mypal[c(2,1)], col="Adormut", title="Neat1 expression versus mutation status in 84 Liver Patients")
g + stat_compare_means()
dev.off()

#Compare expression 
pdf("ADOR2AS1_MutationvsExpression_84Liver.pdf", width=9)
g  = ggboxplot(lnc_rna, x="Adormut", y="ADORA2A", add="jitter", palette=mypal[c(2,1)], col="Neat1mut", title="ADORA2A-AS1 expression versus mutation status in 84 Liver Patients")
g + stat_compare_means()
dev.off()


















