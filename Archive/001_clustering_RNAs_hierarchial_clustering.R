library(survAUC)
#source("source_code_Cox_MonteCarlo_CV_Mar13.R")
require(caTools)
#check if this person is in my analysis: TCGA-61-2095
library(glmnet)
library(survcomp)
library(caret)
library(stringr)
library(EnvStats)
library(patchwork)

source("universal_LASSO_survival_script.R")

library(ggpubr)
library(ggrepel)
library(viridis)
library(patchwork)
library(caret)  
library(Rtsne)
library(data.table)
library(GenomicRanges)

#------DATA---------------------------------------------------------
#UCSC gene info
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
#z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
#ucsc <- ucsc[z,]
#z <- which(duplicated(ucsc[,8]))
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

#save RNA and PCG files locally
#saveRDS(rna, file="rna_lncRNAs_expression_data_june29.rds")
#saveRDS(pcg, file="rna_pcg_expression_data_june29.rds")

setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNAs_2019_manuscript")

#------FEATURES-----------------------------------------------------

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
#allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_Aug8.rds")
allCands = subset(allCands, data == "TCGA") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 
allCands$combo = unique(paste(allCands$gene, allCands$cancer, sep="_"))


mypal = c("#E5DFD9","#EAD286" ,"#D1EB7B", "#96897F" ,"#E5C0A6" ,
  "#72A93B", "#74DAE3" ,"#49B98D" ,"#D97B8F" ,"#70A2A4", "#64709B" ,"#DFBF38" ,"#61EA4F" ,
  "#C7CBE7", "#786DDA",
"#CFA0E0" ,"#67E9D0" ,"#7C9BE1", "#D94753" ,
"#AAE6B0", "#D13BDF" ,"#DEAEC7" ,"#BBE6DF" ,"#B2B47A" ,"#E6ECBA", "#C86ED7",
 "#7BEE95" ,"#6F46E6" ,"#65B9E0", "#C0EC3E",
"#DE8D54" ,"#DF4FA6")

#-------------------------------------------------------------------
#-----------------PCA using just all expressed lncRNA --------------
#-------------------------------------------------------------------

#remove cancer types with less than 50 patients 
pats_num = as.data.table(table(rna$Cancer))
pats_num = filter(pats_num, N <50)
canc_rm = pats_num$V1

#remove those ones
cancers = unlist(cancers[which(!(cancers %in% canc_rm))])

#rna = pcg
rna = subset(rna, Cancer %in% cancers)

dim(rna)
rownames(rna) = rna$patient
z1 = which(str_detect(colnames(rna), "ENSG"))
z2 = which(colnames(rna) %in% "type")
rna = as.data.frame(rna)
rna = rna[,c(z1, z2)]

#1. using medians of each lncRNA in each cancer type 
#cluster cancer types by lncRNA expression 

get_meds = function(canc){
	dat = as.data.table(filter(rna, type == canc))
	#first scale values 	
	z1 = which(str_detect(colnames(dat), "ENSG"))
	meds = apply(dat[,..z1], 2, median)
	meds_dat = as.data.frame(meds)
	meds_dat$lnc = names(meds)
	meds_dat$type = canc
	return(meds_dat)
}

cancers = unique(rna$type)
all_meds = llply(cancers, get_meds)
all_meds = ldply(all_meds)
rownames(all_meds) = all_meds$type

library(tidyr)
data_wide <- spread(all_meds, lnc, meds)
rownames(data_wide) = data_wide$type
data_wide$type = NULL

data_wide <- scale(data_wide)
head(data_wide)

# Dissimilarity matrix
d <- dist(data_wide, method = "manhattan")

# Hierarchical clustering using Complete Linkage
hc1 <- hclust(d, method = "ward.D2" )

# Plot the obtained dendrogram
pdf("lncRNAs_figure1a_hierarchial_clustering.pdf")
plot(hc1, hang = -1)
dev.off()

#mRNA 

get_meds = function(canc){
	dat = as.data.table(filter(pcg, type == canc))
	#first scale values 	
	z1 = which(str_detect(colnames(dat), "ENSG"))
	meds = apply(dat[,..z1], 2, median)
	meds_dat = as.data.frame(meds)
	meds_dat$lnc = names(meds)
	meds_dat$type = canc
	return(meds_dat)
}

cancers = unique(rna$type)
all_meds = llply(cancers, get_meds)
all_meds = ldply(all_meds)
rownames(all_meds) = all_meds$type

library(tidyr)
data_wide <- spread(all_meds, lnc, meds)
rownames(data_wide) = data_wide$type
data_wide$type = NULL

data_wide <- scale(data_wide)
head(data_wide)

# Dissimilarity matrix
d <- dist(data_wide, method = "manhattan")

# Hierarchical clustering using Complete Linkage
hc1 <- hclust(d, method = "ward.D2" )

# Plot the obtained dendrogram
pdf("pcgs_figure1a_hierarchial_clustering.pdf")
plot(hc1, hang = -1)
dev.off()






