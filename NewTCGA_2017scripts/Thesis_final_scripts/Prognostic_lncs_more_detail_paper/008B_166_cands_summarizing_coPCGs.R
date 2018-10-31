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
library(reshape2)

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

pcgs_risks_both = readRDS("pcgs_enriched_in_risk_groups_non_lncRNA_risk_groups_pcg_analysis_july24.rds")
risk = as.data.table(filter(pcgs_risks_both, V2 == "Risk"))
nonrisk = as.data.table(filter(pcgs_risks_both, V2 == "NonRisk"))

all_de_results = readRDS("coexpression_results_processed_july24.rds")
all_de_results = as.data.table(all_de_results)

#using count data
all_de_results = readRDS("diff_expressed_PCGs_lncRNA_risk_groups_Aug21.rds")
all_de_results = as.data.table(all_de_results)


#------make matrix of lncRNA candidates within each cancer type
#and their associated PCGs

cancers = unique(all_de_results$cancer)

make_matrix_for_ap = function(canc){
	dat = dplyr::filter(all_de_results, cancer %in% canc)
	dat$gene_name = as.character(dat$gene_name)

	#columns -> lnc
	#row -> gene_name
	#cell -> P.value 
	#upregulated genes 
	dat_up = filter(dat, pcg_risk == "upregulated_in_risk")
	dat_up_matrix = acast(dat_up, ID~lnc, value.var="P.Value")
	dat_up_matrix[is.na(dat_up_matrix)] = 1
	#file = paste("Aug22_DE_genes_fActivePathways/", canc, "upregulated_risk_genes_matrix_for_activepathways_Aug22.rds", sep="_")
	#saveRDS(dat_up_matrix, file)

	#downregulated genes
	dat_down = filter(dat, pcg_risk == "downregulated_in_risk")
	dat_down_matrix = acast(dat_down, ID~lnc, value.var="P.Value")
	dat_down_matrix[is.na(dat_down_matrix)] = 1
	#file = paste("Aug22_DE_genes_fActivePathways/", canc, "downregulated_risk_genes_matrix_for_activepathways_Aug22.rds", sep="_")
	#saveRDS(dat_down_matrix, file)

	#make lists of genes for pathway enrichment analysis 
	#upregulated
	upreg_genes = as.data.table(filter(dat_up, adj.P.Val <= 0.05))
	lncs = unique(dat$lnc)
	if(length(lncs) > 1){
		sum = as.data.table(table(upreg_genes$gene_name))
		sum = as.data.table(filter(sum, N >= 1))
		upreg_gene_list = unique(sum$V1)
	}
	if(!(length(lncs) > 1)){
		upreg_gene_list = unique(upreg_genes$gene_name)
	}

	if(!(length(upreg_gene_list)==0)){
	#pathways 
	combined_paths <- gprofiler(upreg_gene_list, organism = "hsapiens", exclude_iea=TRUE, ordered_query= TRUE, min_set_size=5, max_set_size = 500, min_isect_size=5, correction_method="fdr")
	print(dim(combined_paths)[1])

		if(!(dim(combined_paths)[1]==0)){
		#only keep GO or REACTOME
		reac <- grep("REAC", combined_paths$term.id)
		#go <- grep("GO", combined_paths$term.id)
		combined_paths <- combined_paths[c(reac), ]
		combined_paths <- combined_paths[,c(9,12, 3, 3, 1, 14)]
		colnames(combined_paths) <- c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Genes")
		file = paste("Aug22_DE_genes_gProfiler_results/", canc, "upregulated_risk_genes_pathways.txt", sep="_")
		write.table(combined_paths, sep= "\t", file, quote=F, row.names=F)
		up_path = length(unique(combined_paths$Description))

		}
	}

	#downregulated
	downreg_genes = as.data.table(filter(dat_down, adj.P.Val <= 0.05))
	if(length(lncs) > 1){
		sum = as.data.table(table(downreg_genes$ID))
		sum = as.data.table(filter(sum, N >= 1))
		downreg_gene_list = unique(sum$V1)
	}
	if(!(length(lncs) > 1)){
		downreg_gene_list = unique(downreg_genes$ID)
	}

	if(!(length(downreg_gene_list)==0)){
	#pathways 
	combined_paths <- gprofiler(downreg_gene_list, organism = "hsapiens", exclude_iea=TRUE, ordered_query= TRUE, min_set_size=5, max_set_size = 500, min_isect_size=5, correction_method="fdr")
	print(dim(combined_paths)[1])

		if(!(dim(combined_paths)[1]==0)){
		#only keep GO or REACTOME
		reac <- grep("REAC", combined_paths$term.id)
		#go <- grep("GO", combined_paths$term.id)
		combined_paths <- combined_paths[c(reac), ]
		combined_paths <- combined_paths[,c(9,12, 3, 3, 1, 14)]
		colnames(combined_paths) <- c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Genes")
		file = paste("Aug22_DE_genes_gProfiler_results/", canc, "downregulated_risk_genes_pathways.txt", sep="_")
		write.table(combined_paths, sep= "\t", file, quote=F, row.names=F)
		down_path = length(unique(combined_paths$Description))

		}
	}

	print(paste(canc, "done"))
}

llply(cancers, make_matrix_for_ap, .progress="text")

#done gave Marta matriced produced by this code for ActivePathways 



##---------pathways enriched by PCGs that appear in at least 2 cancers and in risk group
#remove those that are in both risk and non-risk ...? 

risk = as.data.table(filter(risk, both_risk_groups == "", N >=2))
genes = risk$V1
combined_paths <- gprofiler(genes, organism = "hsapiens", exclude_iea=TRUE, ordered_query= TRUE, min_set_size=5, max_set_size = 200, min_isect_size=5, correction_method="fdr")
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

#barplot: x-axis pathway, y-axis FDR + rotate plot
combined_paths$num_genes =  sapply(combined_paths$Genes, function(x){length(unlist(strsplit(x, ",")))})
combined_paths = as.data.table(combined_paths) 
#order by first num of genes then FDR
combined_paths = combined_paths[order(num_genes)]
combined_paths$Description = factor(combined_paths$Description, levels=combined_paths$Description)

pdf("FINAL_figure6E_risk_pathways.pdf")
g = ggplot(combined_paths, aes(num_genes, Description, color=FDR)) + xlab("Number of Risk Genes") + ylab("Reactome Pathway") + 
        geom_point() + scale_colour_gradientn(colours = terrain.colors(10))
ggpar(g, font.ytickslab = c(6, "plain", "black"))
dev.off()

#enriched in non-risk
nonrisk = as.data.table(filter(nonrisk, both_risk_groups == "", N >=2))
genes = nonrisk$V1
combined_paths <- gprofiler(genes, organism = "hsapiens", exclude_iea=TRUE, ordered_query= TRUE, min_set_size=5, max_set_size = 200, min_isect_size=5, correction_method="fdr")
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

#barplot: x-axis pathway, y-axis FDR + rotate plot
combined_paths$num_genes =  sapply(combined_paths$Genes, function(x){length(unlist(strsplit(x, ",")))})
combined_paths = as.data.table(combined_paths) 
#order by first num of genes then FDR
combined_paths = combined_paths[order(num_genes)]
combined_paths$Description = factor(combined_paths$Description, levels=combined_paths$Description)

pdf("FINAL_figure6E_NONrisk_pathways.pdf", height=3)
g = ggplot(combined_paths, aes(num_genes, Description, color=FDR)) + xlab("Number of Non-Risk Genes") + ylab("Reactome Pathway") + 
        geom_point() + scale_colour_gradientn(colours = terrain.colors(10))
ggpar(g, font.ytickslab = c(6, "plain", "black"))
dev.off()

