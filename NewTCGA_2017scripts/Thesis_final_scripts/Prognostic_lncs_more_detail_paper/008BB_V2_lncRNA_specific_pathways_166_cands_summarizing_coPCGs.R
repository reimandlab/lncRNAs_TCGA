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
val_cands = subset(val_cands, as.numeric(pval) < 0.05)

#------PCGs for pathway analysis--------------------------------------

pcgs_risks_both = readRDS("pcgs_enriched_in_risk_groups_non_lncRNA_risk_groups_pcg_analysis_july24.rds")
risk = as.data.table(filter(pcgs_risks_both, V2 == "Risk"))
nonrisk = as.data.table(filter(pcgs_risks_both, V2 == "NonRisk"))

all_de_results = readRDS("coexpression_results_processed_july24.rds")
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
	
	#all genes - both up and downregulated genes 

	#upregulated genes 
	dat_all = dat
	dat_all_matrix = acast(dat_all, ID~lnc, value.var="P.Value")
	dat_all_matrix[is.na(dat_all_matrix)] = 1

	file = paste("Aug22_DE_genes_fActivePathways/", canc, "all_up_down_genes_fActivepathways_Sep14.rds", sep="_")
	saveRDS(dat_all_matrix, file)

	lncs = unique(dat$lnc)
	
	lnc_spef_pe = function(ln){	

		#make lists of genes for pathway enrichment analysis 
		#upregulated
		upreg_genes = as.data.table(filter(dat_all, lnc == ln, adj.P.Val <= 0.05))
		upreg_gene_list = unique(upreg_genes$gene)

		if(!(length(upreg_gene_list)==0)){
		#pathways 
		combined_paths <- gprofiler(upreg_gene_list, organism = "hsapiens", exclude_iea=TRUE, ordered_query= TRUE, min_set_size=5, max_set_size = 500, min_isect_size=5, correction_method="fdr")
		print(dim(combined_paths)[1])

		if(!(dim(combined_paths)[1]==0)){
		#only keep GO or REACTOME
		reac <- grep("REAC", combined_paths$term.id)
		go <- grep("GO", combined_paths$term.id)
		combined_paths <- combined_paths[c(reac, go), ]
		combined_paths <- combined_paths[,c(9,12, 3, 3, 1, 14)]
		colnames(combined_paths) <- c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Genes")
		file = paste("Aug22_DE_genes_gProfiler_results/", ln, canc, "all_Sept14_pathways.txt", sep="_")
		write.table(combined_paths, sep= "\t", file, quote=F, row.names=F)
		up_path = length(unique(combined_paths$Description))
		if(!(dim(combined_paths)[1]==0)){
		combined_paths$lnc = ln
		return(combined_paths)

		}
		} #if(!(dim(combined_paths)[1]==0)){
		} #if(!(length(upreg_gene_list)==0)){
		}

	all_lncs_paths = llply(lncs, lnc_spef_pe, .progress="text")
	all_lncs_paths = ldply(all_lncs_paths)
	
	if(!(dim(all_lncs_paths)[1]==0)){
	all_lncs_paths$canc = canc
	print(paste(canc, "done"))
	return(all_lncs_paths)
	}
}

all_lnc_pathways = llply(cancers, make_matrix_for_ap, .progress="text")
all_lnc_pathways_df = ldply(all_lnc_pathways)
#done gave Marta matriced produced by this code for ActivePathways 


###---------------Summary figure-------------------------------###

saveRDS(all_lnc_pathways_df, file="pathways_for_each_lncRNA_Sept14.rds")

#number of PCGs/lncRNA vs number of Pathways/lncRNAs 

#1. DE PCGs/lncRNA
sig_des = as.data.table(filter(all_lnc_pathways_df, adj.P.Val <= 0.05))
sig_des_sum = as.data.table(table(sig_des$combo))
sig_des_sum = sig_des_sum[order(N)]

#keep those wtih at least 20pcgs
sig_des_sum = as.data.table(filter(sig_des_sum, N > 20))
colnames(sig_des_sum) = c("combo", "num_sig_des")

#2. pathways/lncRNA
all_lnc_pathways_df$combo = paste(all_lnc_pathways_df$lnc, all_lnc_pathways_df$canc, sep="_")
all_lnc_pathways_df = as.data.table(filter(all_lnc_pathways_df, FDR <= 0.05))
sig_paths_sum = as.data.table(table(all_lnc_pathways_df$combo))
sig_paths_sum = sig_paths_sum[order(N)]

#keep those with least 5 pathways 
sig_paths_sum = as.data.table(filter(sig_paths_sum, N > 5))
colnames(sig_paths_sum) = c("combo", "num_sig_pathways")
sig_paths_sum = unique(sig_paths_sum)

#3. combine 
sig_paths_sum = merge(sig_paths_sum, sig_des_sum, by="combo")
sig_paths_sum$canc = sapply(sig_paths_sum$combo, function(x){unlist(strsplit(x, "_"))[2]})

canc_conv = rna[,c("type", "Cancer")]
canc_conv = unique(canc_conv)
colnames(canc_conv)[2] = "canc"
sig_paths_sum = merge(sig_paths_sum, canc_conv, by="canc")


pdf("summary_#DE_pcgs_vs_pathways_pathways_figure_sep14.pdf", width=9, height=5)
g = ggscatter(sig_paths_sum, x = "num_sig_des", y = "num_sig_pathways",
   fill = "type", shape = 21, size = 3, # Points color, shape and size
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   cor.coeff.args = list(method = "spearman"), xlab="# of Sig DE PCGs", ylab="# of Sig Pathways")
ggpar(g, legend="right")
dev.off()


t = as.data.table(table(sig_paths_sum$type))
t=t[order(N)]

ggbarplot(t, "V1", "N",
   fill = "V1")
dev.off()







