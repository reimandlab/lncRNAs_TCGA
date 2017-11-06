###compilingMiniMethylationFiles.R
#+++++++++++++++++++++++++++++++++

#Purpose: this script will iteriterate over the mini
#methylation files to extract probes for 
#ovarian cancer patients, overlapping lncRNA regions 

###Load Libraries
#+++++++++++++++++++++++++++++++++

source("source_file.R")
library(GenomicRanges)

###Data
#+++++++++++++++++++++++++++++++++

canc_pats = fread("485_patient_IDs_top5CancersPCAWG.txt")

lnc_probes = read.table("cand_lincs_wMethylationProbes.txt")
colnames(lnc_probes) = c("Chr_lnc", "start_lnc", "end_lnc", "ensg", "hugo", "type", 
	"transcript", "chr_probe", "start_probe", "end_probe", "probe")

###Cat methylation results into one file 
results = fread("merged_methylation_files_for8candidates_top_cancers.txt", fill=TRUE)
z <- which(results$V2 == "V2")
results = results[-z,]

#z <- which(results$V2 == "project_code")
#results = results[-z,]

###which lncRNAs are covered by probes 
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

###Merge methylation file with lnc-RNA probes so we know 
###which probe corresponds to which gene 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++

cols = c(2, 3, 9, 10, 11)
results = as.data.frame(results)
results = results[,cols]
colnames(results) = c("patient", "country", "probe", "beta_value", "measure")
results = merge(results, lnc_probes, by = "probe")

###Get lncRNA expression data so we can assign 
###high or low lncRNA candidate expression to each patient 
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++

lnc_rna <- readRDS("5607_pcawg_lncRNAs_RNASeq_data.rds")
lnc_rna <- as.data.frame(lnc_rna)
lnc_rna$patient <- rownames(lnc_rna)
z <- which(lnc_rna$patient %in% canc_pats$patient)
lnc_rna = lnc_rna[z,]

canc = c("Kidney Adenocarcinoma, clear cell type", "Ovary Serous cystadenocarcinoma", "Liver Hepatocellular carcinoma", "Pancreas Pancreatic ductal carcinoma")
z <- which(lnc_rna$canc %in% canc)
lnc_rna = lnc_rna[z,]

cancers = unique(lnc_rna$canc)
get_canc_data = function(canc){
	z <- which(lnc_rna$canc == canc)
	canc_data = lnc_rna[z,]
	#canc specific lncrnas 
	canc_data$canc = str_sub(canc_data$canc, 1, 4)
	genes = lncs$gene[which(lncs$canc == canc_data$canc[1])]
	for(i in 1:length(genes)){
		z <- which(ucsc$hg19.ensGene.name2 == genes[[i]])
		genes[[i]] = ucsc$hg19.ensemblToGeneName.value[z]
	}
	z <- which(colnames(canc_data) %in% genes)
	canc_data = canc_data[,c(z, 5608, 5609)]
	return(canc_data)	
}

canc_data_list = llply(cancers, get_canc_data)
















