###compilingMiniMethylationFiles.R
#+++++++++++++++++++++++++++++++++

#Purpose: this script will iteriterate over the mini
#methylation files to extract probes for 
#ovarian cancer patients, overlapping lncRNA regions 

###Load Libraries
#+++++++++++++++++++++++++++++++++

source("source_file.R")
library(GenomicRanges)
library(stringr)

###Data
#+++++++++++++++++++++++++++++++++

canc_pats = fread("485_patient_IDs_top5CancersPCAWG.txt")

lnc_probes = read.table("cand_lincs_wMethylationProbes.txt")
colnames(lnc_probes) = c("Chr_lnc", "start_lnc", "end_lnc", "ensg", "hugo", "type", 
	"transcript", "chr_probe", "start_probe", "end_probe", "probe")

###Cat methylation results into one file 
results = fread("merged_KiRC_lncRNA_files.txt", fill=TRUE)
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
z <- which(results$V2 == "project_code")
if(!length(z)==0){results = results[-z,]}
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

#Liver candidates 
Kidn = filter(lncs, canc == "Kidn")
for(i in 1:nrow(Kidn)){
	z <- which(lincs$hg19.ensGene.name2 == Kidn$gene[i])
	Kidn$gene[i] = lincs$hg19.ensemblToGeneName.value[z]	
}

z <- which(results$hugo %in% Kidn$gene)
results = results[z,]

#most patients have two entried for each probe? why...not sure 

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Add tag to each patient for each lncRNA based on high or low
#expression 

#canc = c("Kidney Adenocarcinoma, clear cell type", "Ovary Serous cystadenocarcinoma", "Liver Hepatocellular carcinoma", "Pancreas Pancreatic ductal carcinoma")
canc = "Kidney Adenocarcinoma, clear cell type"
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

canc_data = get_canc_data(canc)

check = function(vector_value, median){
	check = vector_value >=median
		if(check){
			return("High")
		}
		if(!(check)){
			return("Low")
		}
	}

#Only 42 patients have methylation data 
tag_exp = function(canc_data){
	canc_data_tagged = canc_data
	#how many columns/genes are there to evaluate 
	cols = ncol(canc_data_tagged)-2
	for(i in 1:cols){
		gene = colnames(canc_data_tagged)[i]
		median = median(as.numeric(canc_data_tagged[,i]))
		tags = unlist(llply(canc_data_tagged[,i], check, median=median))
		canc_data_tagged[,i] = tags
	}
	return(canc_data_tagged)
}

tagged_canc_data = tag_exp(canc_data)

#Methlyation merge with expression data 
merged_methylation_exp = merge(results, tagged_canc_data, by = "patient")
merged_methylation_exp$beta_value = as.numeric(merged_methylation_exp$beta_value)

pdf("AC006126.4_expressionVSmethylation_KIRC.pdf", width=9, height=9)
g = ggboxplot(merged_methylation_exp, x = "AC006126.4", y="beta_value", add ="jitter", palette=mypal, fill="AC006126.4")
g = facet(g, facet.by = "probe")
g = g + stat_compare_means()
print(g)
dev.off()









