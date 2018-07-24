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

#------DATA---------------------------------------------------------
#UCSC gene info
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
#z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
#ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,8]))
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

#save RNA and PCG files locally
#saveRDS(rna, file="rna_lncRNAs_expression_data_june29.rds")
#saveRDS(pcg, file="rna_pcg_expression_data_june29.rds")

rna = readRDS("rna_lncRNAs_expression_data_june29.rds")
pcg = readRDS("rna_pcg_expression_data_june29.rds")

#------FEATURES-----------------------------------------------------

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = subset(allCands, data == "TCGA") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 
allCands$combo = unique(paste(allCands$gene, allCands$cancer, sep="_"))
cands_dups = unique(allCands$gene[which(duplicated(allCands$gene))])

val_cands = read.csv("175_lncRNA_cancers_combos_23_cancer_types_july5.csv")
val_cands = as.data.table(val_cands)
val_cands = subset(val_cands, data == "PCAWG") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 
val_cands$combo = unique(paste(val_cands$gene, val_cands$cancer, sep="_"))
val_cands = subset(val_cands, top_pcawg_val == "YES") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 

#------PCGs for pathway analysis--------------------------------------

pcgs_risks = readRDS("pcgs_enriched_in_risk_groups_non_lncRNA_risk_groups_pcg_analysis_july24.rds")
risk = as.data.table(filter(pcgs_risks, V2 == "Risk"))
nonrisk = as.data.table(filter(pcgs_risks, V2 == "NonRisk"))

##---------pathways enriched by PCGs that appear in at least 2 cancers and in risk group
#remove those that are in both risk and non-risk ...? 

risk = as.data.table(filter(risk, both_risk_groups == "", N >=5))
genes = risk$V1
combined_paths <- gprofiler(genes, organism = "hsapiens", exclude_iea=TRUE, ordered_query= TRUE, min_set_size=5, max_set_size = 300, min_isect_size=10, correction_method="fdr")
print(dim(combined_paths)[1])

if(!(dim(combined_paths)[1]==0)){
#only keep GO or REACTOME
reac <- grep("REAC", combined_paths$term.id)
#go <- grep("GO", combined_paths$term.id)
combined_paths <- combined_paths[c(reac), ]
combined_paths <- combined_paths[,c(9,12, 3, 3, 1, 14)]
colnames(combined_paths) <- c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Genes")
write.table(combined_paths, sep= "\t", file="pathways_enriched_in_PCGS_enriched_in_lncRNA_risk_groups_july13.txt", quote=F, row.names=F)
}

#enriched in non-risk
nonrisk = as.data.table(filter(nonrisk, both_risk_groups == "", N >=5))
genes = nonrisk$V1
combined_paths <- gprofiler(genes, organism = "hsapiens", exclude_iea=TRUE, ordered_query= TRUE, min_set_size=5, max_set_size = 300, min_isect_size=20, correction_method="fdr")
print(dim(combined_paths)[1])
if(!(dim(combined_paths)[1]==0)){
#only keep GO or REACTOME
reac <- grep("REAC", combined_paths$term.id)
#go <- grep("GO", combined_paths$term.id)
combined_paths <- combined_paths[c(reac), ]
combined_paths <- combined_paths[,c(9,12, 3, 3, 1, 14)]
colnames(combined_paths) <- c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Genes")
write.table(combined_paths, sep= "\t", file="pathways_enriched_in_PCGS_enriched_in_lncRNA_NONrisk_groups_july13.txt", quote=F, row.names=F)
}



