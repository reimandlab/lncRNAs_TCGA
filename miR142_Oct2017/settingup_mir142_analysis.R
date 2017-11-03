###settingup_mir142_analysis.R

###Libraries
source("source_file.R")
library(GenomicRanges)
library(stringr)

###Data 

#1. ppl with mutations in miR142
load("mirna142_muts_for_Lina_250117.rsav")
#20 mutations, 17 people in total 
muts_ppl = mirna142_muts$mirna_pre 

#2. Gene Expression data
exp = readRDS("180lymphomaPatients_RNA-Seqfile.rds")

#3. linker file 
link = read.csv("PCAWGMay2016DataReleasev1_v1.1_v1.2_v1.3_v1.4release_may2016v1.4.csv")

#4. cancer subtype 
canc = fread("pcawg_specimen_histology_August2016_v6.tsv")

lymph = canc$icgc_donor_id[canc$histology_tier2=="Lymphoid"]

#5. all samples that have gone through mutation analysis 
#*you should remove any samples not found in this list from your analysis, so that downstream statistical testing is fair
mutanalysis = fread("mutated_samples.txt", header=F)

#6. identified mir142 targets comrpise set of genes tests with annotations 
genes = fread("file_for_network_w_target_type_and_cancer_genes.txt")

#subset exp 
z <- which(exp$Hugo %in% genes$all.targets)
exp = exp[z,]

#7. copy number data for lymphoma patients (196/208) from JR's file obtained on October 30th 
cna = readRDS("196pcawg_lymphoma_cna_matrix_fromJR_oct30.rds") 

#8. UCSC coordinates 
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)

###Converting RNA-Seq specimen ids to donor ids 
#RNA seq specimen IDs to donor ID 

###CNA figure out how segments match to gene 
segments = rownames(cna) 
segments = str_split_fixed(segments, "::", n=2)

#make sure it is still atttached to patients so that later doesn't get lost, order of cnas...
cna = cbind(segments, cna)

#some rows are segments and some are gene IDs, let's divide those into seperate dataframes
segments1 = cna[which(cna[,1] == "ENCODEmerge"),]
segments2 =  cna[which(!(cna[,1] == "ENCODEmerge")),]

#segments 1
chr = str_split_fixed(segments1[,2], ":", n=4)
coordinates = str_split_fixed(chr[,2], "-", n=2)
chr = chr[,1]
segments1 = cbind(segments1, chr, coordinates)
segments1 = as.data.frame(segments1)
segments1$strand = "*"

#segments 2
chr = str_split_fixed(segments2[,2], ":", n=4)
coordinates = str_split_fixed(chr[,2], "-", n=2)
chr = chr[,1]
segments2 = cbind(segments2, chr, coordinates)
segments2 = as.data.frame(segments2)
segments2$strand = "*"

#seperate gene ID from the rest 
gene = str_split_fixed(segments2[,2], "::", n=3)
gene = cbind(segments2, gene)
gene = gene[which(gene[,203] %in% c("gencode")),]

#basically ignoring segments 1 from now as they are just segments 
#in final CNA file, gene -- 25,687 unique Hugo Gene Symbols
#25765 unique Ensembl gene ids 

##subset to genes we are interested in 
gene <- gene[which(gene[,204] %in% genes$all.targets), ] #960/964 genes from list with CNA data 
#remove duplicated rows
z <- which(duplicated(gene[,204]))
gene = gene[-z,]

#which patients in this CNA file have mutations/not have mutations
#add gene expression values to genes 
#merge clinical, expression and cna 

##-----CNA-matrix--------------------------------------------------------------------------
gene = gene[,-(c(1,2,200:203))]
gene$chr = NULL
colnames(gene)[197:198] = c("HUGO", "ensg")

#common patients to CNA and exp 
z <- which(colnames(exp) %in% colnames(gene))
exp = exp[,c(1,z,182)]
gene = gene[,c(which(colnames(gene) %in% colnames(exp)), 197:198)]

###add patient info - mutations status-----------------------------------------------------
mutationStatus = c(1:172)
names(mutationStatus) = colnames(gene)[1:172]

for(i in 1:172){
	pat = colnames(gene)[i]
	z <- which(muts_ppl$mut_patient %in% pat)
	if(!(length(z)==0)){
		mutationStatus[i] = "mutation"
	}
	if(length(z)==0){
		z = which(mutanalysis$V1 %in% pat)
		if(!(length(z)==0)){
		mutationStatus[i] = "Nomutation"
	}
		if(length(z)==0){
			mutationStatus[i] == "NoTested"
		}
	}
}

#all of these patients were analyzed for mutations and have lymphoma
#and have CNA data for 959/964 genes 

#Subset RNA-Seq to these 960 genes
z = which(duplicated(exp$Hugo))
z2 = which(exp$Hugo %in% exp$Hugo[z])
exp = exp[-z2,]

z = which(exp$Hugo %in% gene$HUGO)
exp = exp[z,]

z = which(gene$HUGO %in% exp$Hugo)
gene = gene[z,]

#add histological subtype 
###add patient info - mutations status-----------------------------------------------------
lymphsubtype = c(1:172)
names(lymphsubtype) = colnames(gene)[1:172]

for(i in 1:172){
	pat = colnames(gene)[i]
	z = which(canc$icgc_donor_id %in% pat)
	if(length(z) > 1){
		z = z[[1]]
	}
	histo = canc$histology_tier3[z]
	lymphsubtype[i] = histo
}

rownames(gene) = gene$HUGO
gene = gene[,1:172]

gene = rbind(gene, mutationStatus, lymphsubtype)

rownames(exp) = exp$Hugo
exp = exp[,2:173]

##########################SAVE ALL DATA------------------------------------------------------

#CNA 172 patients
saveRDS(gene, "172_Lymph_CNA_plus_MutStatus_Histo.rds")

#RNA-Seq 172 patients 
saveRDS(exp, "172_Lymph_RNASeq.rds")



























