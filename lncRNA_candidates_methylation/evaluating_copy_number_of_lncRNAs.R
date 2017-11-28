###evaluating_copy_number_of_lncRNAs.R
#++++++++++++++++++++++++++++++++++++++++

#Purpose: this script will iteriterate over the mini
#methylation files to extract probes for 
#ovarian cancer patients, overlapping lncRNA regions 

###Load Libraries
#++++++++++++++++++++++++++++++++++++++++

source("source_file.R")
library(GenomicRanges)
library(stringr)

###Data
#++++++++++++++++++++++++++++++++++++++++
ov_pats = fread("ovarianPatientsRNAseqPCAWGn=70.txt")

##which lncRNAs are covered by probes 
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

#RNA-Seq file 
#lncRNA
lnc_rna <- readRDS("5607_pcawg_lncRNAs_RNASeq_data.rds")
lnc_rna <- as.data.frame(lnc_rna)
lnc_rna$patient <- rownames(lnc_rna)

#PCGs
pcg_rna <- readRDS("20166_pcawg_PCGs_RNASeq_data.rds")
pcg_rna <- as.data.frame(pcg_rna)
pcg_rna$patient <- rownames(pcg_rna)

#remove duplicated column names 
dups <- colnames(pcg_rna)[which(duplicated(colnames(pcg_rna)))]   
#save them in a list for future reference 
pcg_rna <- pcg_rna[,-(which(colnames(pcg_rna) %in% dups))]

#Clinical file - available only for 485/497 patients 
clin_pats <- readRDS("Jan26_PCAWG_clinical")
z <- which(clin_pats$icgc_donor_id %in% ov_pats$V2)
clin_pats <- clin_pats[z,]

lnc_rna <- lnc_rna[which(rownames(lnc_rna) %in% clin_pats$icgc_donor_id),] 
pcg_rna <- pcg_rna[which(rownames(pcg_rna) %in% clin_pats$icgc_donor_id),] 

#make sure there aren't any lncRNAs in protein coding file 
table(ucsc[,7][ucsc[,8] %in% colnames(pcg_rna)])
#Remove 
z <- which(colnames(pcg_rna) %in% fantom[,2])
pcg_rna <- pcg_rna[,-z]

#Candidate lncRNAs 
#List of canddidates and cox results
allCands <- fread("lncRNAs_sig_FDR_0.1_Nov23.txt")
allCands = filter(allCands, gene %in% c("NEAT1", "RP11-622A1.2", "GS1-251I9.4", "ZNF503-AS2", "AC009336.24"))
allCands = allCands[-6,]
lncs = filter(allCands, canc == "Ovary Serous cystadenocarcinoma")

#ucsc only keep those lncrnas in fantom
ucsc = ucsc[which(ucsc$hg19.ensemblToGeneName.value %in% lncs$gene),]

#lncRNA candidate expression status 
medianexp = read.table("medianScoresOvarianCancerTop3_lncRNAs.txt", sep=";", header=T)
#mutations from MAF file 
maf = fread("merged_mutations_ovarian_cancer_patients.txt", fill=TRUE)
z = which(maf$V43 == "V43")
maf = maf[-z,]

colnames(maf)[15] = "patient"
maf = filter(maf, patient %in% ov_pats$V2)

#which patients US which patients AUS
us = fread("sample.OV-AU.tsv")
us = unique(us$icgc_donor_id)
aus = fread("sample.OV-US.tsv")
aus = unique(aus$icgc_donor_id)

maf$country = ""
maf$country[which(maf$patient %in% us)] ="US"
maf$country[which(maf$patient %in% aus)] = "AUS"
maf = maf[,-1]

colnames(maf) = c("Hugo", "Chr", "Start", "End", "Strand", "Variant_Classification", "Variant_Type", "Reference_allele", "Reference_allele2", "alt_allele", "dbSNP_RS", "gc_content",
	"Project_Code", "Donor_ID", "country")

#copy number data for ovarian patients from JR's file obtained on October 30th 
cna = readRDS("ovarianPatients_CNA_data.rds") 

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
gene = gene[which(gene[,73] %in% c("gencode")),]

#are lncRNAs in here?
g = function(row){
	gg = row[[75]]
	gg = gsub("\\..*","",gg)
	return(gg)
}
gene$ensg = apply(gene, 1, g)
z <- which(gene$ensg %in% ucsc$hg19.ensGene.name2)
lncwCNA = gene[z,]

rownames(lncwCNA) = lncwCNA$ensg
lncwCNA = lncwCNA[,-(c(1:2, 69:76))]

###add patient info - expression status-----------------------------------------------------
exp = c(1:66)
names(exp) = colnames(lncwCNA)[1:66]

for(i in 1:66){
	pat = colnames(lncwCNA)[i]
	z <- which(medianexp$patient %in% pat)
	if(!(length(z)==0)){
		exp[i] = medianexp$GS1[z]
	}
}

lncwCNA = rbind(lncwCNA, exp)
lncwCNA = t(lncwCNA)
colnames(lncwCNA) = c("CNA_GS1", "Expression_GS1")
lncwCNA = as.data.frame(lncwCNA)
lncwCNA$patient = rownames(lncwCNA)
lncwCNA$CNA_GS1 = as.numeric(lncwCNA$CNA_GS1)

pdf("GS1-251I9.4_Boxplot_CNA_expressionOct31.pdf", width=9, height=8)
g  = ggboxplot(lncwCNA, x="Expression_GS1", y="CNA_GS1", add="jitter", palette=mypal[c(2,1)], col="Expression_GS1", title="GS1-251I9.4 Expression vs CNA in 66 Ovarian Cancer Patients")
g + stat_compare_means(method = "t.test")
dev.off()

pdf("GS1-251I9.4_Boxplot_CNA_expressionOct31Wilcoxon.pdf", width=9, height=8)
g  = ggboxplot(lncwCNA, x="Expression_GS1", y="CNA_GS1", add="median", palette=mypal[c(2,1)], col="Expression_GS1", title="GS1-251I9.4 Expression vs CNA in 66 Ovarian Cancer Patients")
g + stat_compare_means()
dev.off()


###CNA level of PCG antisense to GS1 lncRNA 
z <- which(gene$ensg %in% "ENSG00000155100")
pcg = gene[z[[1]],]

rownames(pcg) = pcg$ensg
pcg = pcg[,-(c(1:2, 69:76))]

###add patient info - expression status PCG-----------------------------------------------------
exp = c(1:66)
names(exp) = colnames(pcg)[1:66]

for(i in 1:66){
	pat = colnames(pcg)[i]
	z <- which(medianexp$patient %in% pat)
	if(!(length(z)==0)){
		exp[i] = medianexp$GS1[z]
	}
}

pcg = rbind(pcg, exp)
pcg = t(pcg)
colnames(pcg) = c("CNA_OTUD6B", "Expression_GS1")
pcg = as.data.frame(pcg)
pcg$patient = rownames(pcg)
pcg$CNA_OTUD6B = as.numeric(pcg$CNA_OTUD6B)

pdf("OTUD6B_Boxplot_CNA_expressionOct31.pdf", width=9, height=8)
g  = ggboxplot(pcg, x="Expression_GS1", y="CNA_OTUD6B", add="jitter", palette=mypal[c(2,1)], col="Expression_GS1", title="OTUD6B Expression vs CNA in 66 Ovarian Cancer Patients")
g + stat_compare_means(method = "t.test")
dev.off()

pdf("OTUD6B_Boxplot_CNA_expressionOct31Wilcoxon.pdf", width=9, height=8)
g  = ggboxplot(pcg, x="Expression_GS1", y="CNA_OTUD6B", add="jitter", palette=mypal[c(2,1)], col="Expression_GS1", title="OTUD6B Expression vs CNA in 66 Ovarian Cancer Patients")
g + stat_compare_means()
dev.off()























