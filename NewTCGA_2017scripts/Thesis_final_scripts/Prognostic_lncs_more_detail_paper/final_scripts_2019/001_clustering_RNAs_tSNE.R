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

#1. remove those not expressed at all
z1 = which(str_detect(colnames(rna), "ENSG"))
sums = apply(rna[,z1], 2, var) #get 500 most variables genes 
sums = sums[order(-sums)]
keep = sums[1:1000] 
keep = names(keep)
z = which(colnames(rna) %in% c("type", keep))
rna = rna[,z]

#2. remove those with MAD < 0? 
#1. remove those not expressed at all
z1 = which(str_detect(colnames(rna), "ENSG"))
#sums = apply(rna[,z1], 2, mad)
#z = which(sums <= 0)
#rna = rna[,-z]

#2. cluster cancer types, ie -> each dot should be cancer type 
library("FactoMineR")
z1 = which(str_detect(colnames(rna), "ENSG"))
logged_rna = rna
logged_rna[,z1] = log1p(logged_rna[,z1])
#z = which(colnames(logged_rna) %in% rm)
#logged_rna = logged_rna[,-z]
z1 = which(str_detect(colnames(logged_rna), "ENSG"))

#logged - DONE! 
#pdf("logged_nonMAD0lncRNAs_cands_PCA_plots_Aug27.pdf", width=12)
#autoplot(prcomp(logged_rna[,z1]), data = logged_rna, colour = 'type')+ 
#scale_colour_manual(values = mypal)
#dev.off()

#try t-SNE 
library(Rtsne)
#iris_unique <- unique(iris) # Remove duplicates
#iris_matrix <- as.matrix(iris_unique[,1:4])
set.seed(42) # Set a seed if you want reproducible results
tsne_out <- Rtsne(logged_rna[,z1], dims = 2, perplexity=30, verbose=TRUE, max_iter = 500) # Run TSNE
saveRDS(tsne_out, file="tsne_out_28_cancers_Oct1_lncRNAs.rds")

tsne_plot = readRDS("tsne_out_28_cancers_Oct1_lncRNAs.rds")
# Show the objects in the 2D tsne representation
tsne_plot <- data.frame(x = tsne_out$Y[,1], y = tsne_out$Y[,2], col = logged_rna$type)

z = which(duplicated(tsne_plot$col))
tsne_plot$label[z] ="no"
tsne_plot$label[-z] = "yes"

z = (which(tsne_plot$label == "yes"))
for(i in 1:length(z)){
	canc = tsne_plot$col[z[i]]
	tsne_plot$label[z[i]] = canc
}

#assign color to variable ******
pdf("tSNE_28_cancers_april10_top_var_1000_lncRNAs.pdf", width=9)
ggplot(tsne_plot,aes(x, y, label = label)) + geom_point(aes(x=x, y=y, color=col)) + scale_colour_manual(values=mypal)+
geom_text_repel(data = subset(tsne_plot, !(label == "no")))
dev.off()
