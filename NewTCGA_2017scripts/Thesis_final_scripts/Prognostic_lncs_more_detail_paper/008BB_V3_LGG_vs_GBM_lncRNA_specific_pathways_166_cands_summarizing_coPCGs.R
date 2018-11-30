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

#Combined into one dataframe because need to get ranks 
all <- merge(rna, pcg, by = c("patient", "Cancer"))
all = all[,1:25170]

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

all_de_results = readRDS("diff_expressed_PCGs_lncRNA_risk_groups_Aug21.rds")
all_de_results = as.data.table(all_de_results)

all_de_results = readRDS("diff_expressed_PCGs_lncRNA_risk_groups_lgg_nov30.rds")
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

	file = paste("Aug22_DE_genes_fActivePathways/", canc, "all_up_down_genes_fActivepathways_Oct30.rds", sep="_")
	saveRDS(dat_all_matrix, file)

	lncs = unique(dat$lnc)
	
	lnc_spef_pe = function(ln){	

		#make lists of genes for pathway enrichment analysis 
		#upregulated
		upreg_genes = as.data.table(filter(dat_all, lnc == ln, adj.P.Val <= 0.05))
		upreg_gene_list = unique(upreg_genes$gene)

		if(!(length(upreg_gene_list)==0)){
		#pathways 
		combined_paths <- gprofiler(upreg_gene_list, organism = "hsapiens", exclude_iea=TRUE, ordered_query= TRUE, min_set_size=10, max_set_size = 250, min_isect_size=5, correction_method="fdr")
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

#saveRDS(all_lnc_pathways_df, file="pathways_for_each_lncRNA_Oct30.rds")

#all_lnc_pathways_df = readRDS("pathways_for_each_lncRNA_Oct30.rds")

#number of PCGs/lncRNA vs number of Pathways/lncRNAs 

#1. DE PCGs/lncRNA
t = as.data.table(filter(all_de_results, adj.P.Val <= 0.05))
sig_des=t
sig_des$combo = paste(sig_des$lnc, sig_des$canc, sep="_")
sig_des_sum = as.data.table(table(sig_des$combo))
sig_des_sum = sig_des_sum[order(N)]

#keep those wtih at least 20pcgs
#sig_des_sum = as.data.table(filter(sig_des_sum, N > 20))
colnames(sig_des_sum) = c("combo", "num_sig_des")

#2. pathways/lncRNA
all_lnc_pathways_df$combo = paste(all_lnc_pathways_df$lnc, all_lnc_pathways_df$canc, sep="_")
all_lnc_pathways_df = as.data.table(filter(all_lnc_pathways_df, FDR <= 0.01))
sig_paths_sum = as.data.table(table(all_lnc_pathways_df$combo))
sig_paths_sum = sig_paths_sum[order(N)]

#keep those with least 5 pathways 
#sig_paths_sum = as.data.table(filter(sig_paths_sum, N > 5))
colnames(sig_paths_sum) = c("combo", "num_sig_pathways")
sig_paths_sum = unique(sig_paths_sum)

#3. combine 
sig_paths_sum = merge(sig_paths_sum, sig_des_sum, by="combo")
sig_paths_sum$canc = sapply(sig_paths_sum$combo, function(x){unlist(strsplit(x, "_"))[2]})

canc_conv = rna[,c("type", "Cancer")]
canc_conv = unique(canc_conv)
colnames(canc_conv)[1] = "canc"
sig_paths_sum = merge(sig_paths_sum, canc_conv, by="canc")


pdf("summary_#DE_pcgs_vs_pathways_pathways_figure_sep14.pdf", width=9, height=5)
g = ggscatter(sig_paths_sum, x = "num_sig_des", y = "num_sig_pathways",
   fill = "canc", shape = 21, size = 3, # Points color, shape and size
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   cor.coeff.args = list(method = "spearman"), xlab="# of Sig DE PCGs", ylab="# of Sig Pathways")
ggpar(g, legend="right")
dev.off()


t = as.data.table(table(sig_paths_sum$canc))
t=t[order(N)]

ggbarplot(t, "V1", "N",
   fill = "V1")
dev.off()


paths = sig_paths_sum[,c("canc", "combo", "num_sig_pathways", "Cancer")]
paths$type = "pathways"
genes = sig_paths_sum[,c("canc", "combo", "num_sig_des", "Cancer")]
genes$type = "genes"

colnames(paths)[3] = "num_sig_des"
all_res = rbind(paths, genes)
sig_paths_sum = sig_paths_sum[order(-num_sig_des)]
all_res$combo = factor(all_res$combo, levels=unique(sig_paths_sum$combo))

write.csv(sig_paths_sum, file="summary_num_diff_exp_genes_pathways_nov16.csv", quote=F, row.names=F)

#network part A 
pdf("summary_barplot_#DE_pcgs_vs_pathways_pathways_figure_sep14.pdf", width=16, height=5)
g = ggplot(all_res, aes(combo, num_sig_des, group=type)) + theme_classic() + 
   geom_col(position = "dodge", aes(fill=type), color="grey29")+xlab("lncRNA-cancer pair") + ylab("Number of significant \ngenes/pathways")+
   theme(legend.title=element_blank(), legend.position="top", axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
print(g)
dev.off()

#add cancer covariate 
#all_res = merge(all_res, canc_conv, by="canc")
colnames(all_res)[c(1, 4)] = c("type", "canc_name")
all_res$canc_name = factor(all_res$canc_name, levels=unique(sig_paths_sum$Cancer))

library(patchwork)

#20 col palette
pall =  c("#72E1E0" ,"#7CE147", "#73E590", "#E180E3" ,"#D39794", "#D4DC8D", "#A6E6BB", "#A4C9E9",
 "#DC5FA5" ,"#D8DB50" ,"#6891DB" ,"#633CE1" ,"#DDDDCF" ,"#DAAADA" ,"#D849DE",
"#7C9F6A" ,"#E6565F" ,"#8C68D2" ,"#63859B" ,"#DD9955")

pdf("summary_barplot_w_covariate_#DE_pcgs_vs_pathways_pathways_figure_sep14.pdf", width=16, height=5)

#covariate for cancer
cov = ggplot(all_res, aes(combo, 1)) + geom_tile(aes(fill = type)) +
theme_void() + theme(legend.position="none")+
scale_fill_manual(values=pall)

g + cov + plot_layout(ncol = 1, heights = c(15,1))
dev.off()

pdf("canc_cov_legend.pdf")
ggplot(all_res, aes(combo, 1)) + geom_tile(aes(fill = type)) +
theme_void() + coord_flip() + 
scale_fill_manual(values=pall)
dev.off()

#instead of cancer type do gene name
all_res$gene = as.character(all_res$combo)
all_res$gene = sapply(all_res$gene, function(x){unlist(strsplit(x, "_"))[1]})

get_name = function(ensg){
    z = which(fantom$CAT_geneID == ensg)
    return(fantom$CAT_geneName[z][1])
}

all_res$genename = unlist(llply(all_res$gene, get_name))
colnames(all_res)[1] = "canc"
all_res$combo2 = paste(all_res$genename, all_res$canc, sep=" ")

sig_paths_sum$gene = ""
sig_paths_sum$gene = sig_paths_sum$combo
sig_paths_sum$gene = sapply(sig_paths_sum$gene, function(x){unlist(strsplit(x, "_"))[1]})
sig_paths_sum$genename = unlist(llply(sig_paths_sum$gene, get_name))
sig_paths_sum$combo2 = paste(sig_paths_sum$genename, sig_paths_sum$canc, sep=" ")
sig_paths_sum = sig_paths_sum[order(-num_sig_des)]

all_res$combo2 = factor(all_res$combo2, levels = unique(sig_paths_sum$combo2))
all_res$font_col = ""

all_res$font_col[all_res$canc == "LGG"] = "1"

sig_paths_sum$font_col = ""
sig_paths_sum$font_col[sig_paths_sum$canc == "LGG"] = "1"
a <- ifelse(sig_paths_sum$font_col == 1, "red", "black")

all_res = as.data.table(filter(all_res, num_sig_des >=75))

#117 unique lncRNA-cancer pairs 
#114 unique lncRNAs 
pdf("summary_barplot_genenames_#DE_pcgs_vs_pathways_pathways_figure_sep14.pdf", width=12, height=4)
g = ggplot(all_res, aes(combo2, num_sig_des, group=type)) + theme_classic() + 
   geom_col(position = "dodge", aes(fill=type), color="grey29")+xlab("lncRNA-cancer pair") + ylab("Number of sig genes/pathways")+
   theme(legend.title=element_blank(), legend.position="top", axis.title.x=element_blank(), 
   	axis.text.x = element_text(angle = 65, hjust = 1, colour = a, size=6)) + scale_fill_manual(values=c("snow3", "#E69F00", "#56B4E9"))
#ggpar(g,
# font.tickslab = c(6,"plain", "black"),
# xtickslab.rt = 45)
print(g)
dev.off()


###---------------Make heatmap-------------------------------###

#liver gluco related pathways 

load("_LGG_all_up_down_genes_.2018-10-31.rdata")
colnames(res)
#canc = subset(res, res$term.name %in% canc_paths_paths)
canc_paths_paths = res[which(str_detect(res$term.name, "brain")),]$term.name
canc = subset(res, term.name %in% canc_paths_paths)

#1. get all PCGs that are in these pathways 
pcgs_brain = unique(unlist(canc$overlap)) #185 unique PCGs
lncs_brain = unique(unlist(canc$evidence)) #8 unique lncRNAs to be used as covariates high vs low 

##2-----------------label patients by risk------------------------------

dat = subset(all, Cancer=="Brain Lower Grade Glioma")
lgg_nearby_genes = readRDS("lgg_nearby_genes.rds")

res_all_brain = res[,c("term.id", "term.name", "adjusted.p.val", "overlap", "evidence")]
#get count for overlap and for evidence 
for(i in 1:nrow(res_all_brain)){
  over = length(unique(unlist(res_all_brain$overlap[i])))
  res_all_brain$overlap[i] = over
  evidence = length(unique(unlist(res_all_brain$evidence[i])))
  if(evidence>=3){
    res_all_brain$evidence[i] = evidence
  }
}
res_all_brain$overlap = as.character(res_all_brain$overlap)
res_all_brain$evidence = as.character(res_all_brain$evidence)

write.csv(res_all_brain, file="all_brain_pathways_lgg_nov16.csv", quote=F, row.names=F)

#get counts of pathway per lcnRNA
res_all_brain = res[,c("term.id", "term.name", "adjusted.p.val", "overlap", "evidence")]
for(i in 1:nrow(res_all_brain)){
  over = length(unique(unlist(res_all_brain$overlap[i])))
  res_all_brain$overlap[i] = over
  evidence = length(unique(unlist(res_all_brain$evidence[i])))
  res_all_brain$evidence[i] = evidence
}
res_all_brain$overlap = as.character(res_all_brain$overlap)
res_all_brain$evidence = as.character(res_all_brain$evidence)

  dat_keep = dat[,which(colnames(dat) %in% c("patient", lncs_brain, pcgs_brain))]
  rownames(dat_keep) = dat_keep$patient
  dat_keep$patient = NULL

  #heatmap 
  heat = dat_keep
  heat = log1p(heat)
  heat = t(heat)

  #change pcg names
  for(i in 1:nrow(heat)){
    print(i)
    pcg = rownames(heat)[i]
    newname = ucsc$hg19.ensemblToGeneName.value[which(ucsc$hg19.ensGene.name2 == pcg)][1]
    rownames(heat)[i] = newname
  }
  z = which(is.na(rownames(heat)))
  if(!(length(z)==0)){
    heat=heat[-z,]
  }

  #which genes are cancer gene census genes?
  pcgs = unique(rownames(heat))
  #COSMIC cancer gene census
  census = read.csv("Census_allFri_Jul_13_16_55_59_2018.csv")
  #get ensg
  get_census_ensg = function(genes){
  glist = unlist(strsplit(genes, ","))
  z = which(str_detect(glist, "ENSG"))
  ensg = glist[z]
  return(ensg)
  }
  census$ensg = sapply(census$Synonyms, get_census_ensg)

heat = t(heat)

full_heat = heat 

hox_genes = colnames(heat)[which(str_detect(colnames(heat), "HOX"))]
random_genes = c()#which(rownames(heat) %in% c("G6PC", "ALDH5A1", "IGF1", "PDK4", "SLC2A2", "DHDH", "IGF2", "G6PD"))
z1 = which(colnames(heat) %in% c(census$Gene.Symbol))
z2 = which(colnames(heat) %in% ucsc$hg19.ensemblToGeneName.value[which(ucsc$hg19.ensGene.name2 %in% lncs_brain)])
z3 = which(colnames(heat) %in% hox_genes)
heat = heat[,unique(c(z1,z2,z3))]

#get expression of these genes and compare between low lGG, high LGG and GBM 

#1. get list of genes 
development_genes = colnames(heat)
get_ensg_pcg = function(pcg){
  z = which(ucsc$hg19.ensemblToGeneName.value == pcg)
  if(length(z)>1){
    z = z[1]
  }
  return(ucsc$hg19.ensGene.name2[z])
}

development_genes = sapply(development_genes, get_ensg_pcg)

#2. get LGG 
lgg = heat 
lgg = as.data.frame(lgg)
rownames(lgg) = paste(rownames(lgg), "lgg")

#3. get GBM 
z = which(all$type == "GBM")
gbm = all[z,]
z = which(colnames(gbm) %in% c("patient", development_genes))

gbm = gbm[,z]
rownames(gbm) = gbm$patient
gbm$patient = NULL
gbm = t(gbm)

#change pcg names
  for(i in 1:nrow(gbm)){
    print(i)
    pcg = rownames(gbm)[i]
    newname = ucsc$hg19.ensemblToGeneName.value[which(ucsc$hg19.ensGene.name2 == pcg)][1]
    rownames(gbm)[i] = newname
  }
  z = which(is.na(rownames(gbm)))
  if(!(length(z)==0)){
    gbm=gbm[-z,]
  }
gbm = t(gbm)
gbm = log1p(gbm)
gbm = as.data.frame(gbm)
rownames(gbm) = paste(rownames(gbm), "gbm")

#compare expression of these genes between LGG and GBM
gbm <- gbm[names(lgg)]
gbm$type = "gbm"
lgg$type = "lgg"
all_canc = rbind(lgg, gbm)

all_canc$patient = rownames(all_canc)
all_canc$patient = sapply(all_canc$patient, function(x){unlist(strsplit(x, " "))[1]})
all_canc$full_patient = rownames(all_canc)

#lgg info 
lgg_idh = readRDS("TCGA_lgg_wsubtype_info_biolinks.rds")
lgg_idh = lgg_idh[,c("IDH.status", "patient")]

gb_idh = readRDS("TCGA_gbm_wsubtype_info_biolinks.rds")
gb_idh = gb_idh[,c("IDH.status", "patient")]
all_brain = rbind(lgg_idh, gb_idh)

all_canc = merge(all_canc, all_brain, by="patient")
rownames(all_canc) = all_canc$full_patient
all_canc$full_patient = NULL
all_canc$patient = NULL

get_ensg = function(lnc){
  z = which(fantom$CAT_geneName == lnc)
  return(fantom$CAT_geneID[z])
}
z = which(colnames(all_canc) %in% (sapply(lncs_brain, get_name)))
for(i in 1:length(z)){
  colnames(all_canc)[z[i]] = get_ensg(colnames(all_canc)[z[i]])
}
z = which(colnames(all_canc) == "HOXA10-AS")
colnames(all_canc)[z] = "ENSG00000253187"

mat = all_canc[,c(which(!(colnames(all_canc) %in% c("type", "IDH.status"))))]
mat = scale(mat)
mat=t(mat)

df = as.data.frame(all_canc[,c(which(colnames(all_canc) %in% c("type", "IDH.status", "ENSG00000239552", "ENSG00000224950", 
  "ENSG00000250360", "ENSG00000253187", "ENSG00000254635", "ENSG00000255020", "ENSG00000256482", "ENSG00000257261")))])
df$patient = rownames(df)

lnc_med = median(df$ENSG00000254635[df$type=="lgg"])
z1=which((df$ENSG00000254635 >= lnc_med) & (df$type=="lgg"))
z2=which((df$ENSG00000254635 < lnc_med) & (df$type=="lgg"))
df$ENSG00000254635[z1] = 1
df$ENSG00000254635[z2] = 0

lnc_med = median(df$ENSG00000253187[df$type=="lgg"])
z1=which((df$ENSG00000253187 > lnc_med) & (df$type=="lgg"))
z2=which((df$ENSG00000253187 <= lnc_med) & (df$type=="lgg"))
df$ENSG00000253187[z1] = 1
df$ENSG00000253187[z2] = 0

lnc_med = median(df$ENSG00000256482[df$type=="lgg"])
z1=which((df$ENSG00000256482 > lnc_med) & (df$type=="lgg"))
z2=which((df$ENSG00000256482 <= lnc_med) & (df$type=="lgg"))
df$ENSG00000256482[z1] = 1
df$ENSG00000256482[z2] = 0

lnc_med = median(df$ENSG00000224950[df$type=="lgg"])
z1=which((df$ENSG00000224950 >= lnc_med) & (df$type=="lgg"))
z2=which((df$ENSG00000224950 < lnc_med) & (df$type=="lgg"))
df$ENSG00000224950[z1] = 1
df$ENSG00000224950[z2] = 0

lnc_med = median(df$ENSG00000255020[df$type=="lgg"])
z1=which((df$ENSG00000255020 >= lnc_med) & (df$type=="lgg"))
z2=which((df$ENSG00000255020 < lnc_med) & (df$type=="lgg"))
df$ENSG00000255020[z1] = 1
df$ENSG00000255020[z2] = 0

lnc_med = median(df$ENSG00000257261[df$type=="lgg"])
z1=which((df$ENSG00000257261 >= lnc_med) & (df$type=="lgg"))
z2=which((df$ENSG00000257261 < lnc_med) & (df$type=="lgg"))
df$ENSG00000257261[z1] = 1
df$ENSG00000257261[z2] = 0

lnc_med = median(df$ENSG00000239552[df$type=="lgg"])
z1=which((df$ENSG00000239552 > lnc_med) & (df$type=="lgg"))
z2=which((df$ENSG00000239552 <= lnc_med) & (df$type=="lgg"))
df$ENSG00000239552[z1] = 1
df$ENSG00000239552[z2] = 0

lnc_med = median(df$ENSG00000250360[df$type=="lgg"])
z1=which((df$ENSG00000250360 > lnc_med) & (df$type=="lgg"))
z2=which((df$ENSG00000250360 <= lnc_med) & (df$type=="lgg"))
df$ENSG00000250360[z1] = 1
df$ENSG00000250360[z2] = 0

df$ENSG00000254635[df$type=="gbm"]=NA
df$ENSG00000253187[df$type=="gbm"]=NA
df$ENSG00000256482[df$type=="gbm"]=NA
df$ENSG00000224950[df$type=="gbm"]=NA
df$ENSG00000255020[df$type=="gbm"]=NA
df$ENSG00000257261[df$type=="gbm"]=NA
df$ENSG00000239552[df$type=="gbm"]=NA
df$ENSG00000250360[df$type=="gbm"]=NA

#which are favourable lncs? (for those ones flip 0s and 1s)
cands_lncs = as.data.table(filter(allCands, cancer== "Brain Lower Grade Glioma"))
cands_lncs = as.data.table(filter(cands_lncs, HR <1))
cands_lncs = cands_lncs$gene

z = which(colnames(df) %in% cands_lncs)
for(i in 1:length(z)){
  k = which(df[,z[i]] == 0)
  df[k,z[i]] = 1
  df[-k,z[i]] = 0
}

df$risk = apply(df[,1:8], 1, sum)

#df$patient = NULL
#another version of ha

#df = as.data.frame(all_canc[,c(which(colnames(all_canc) %in% c("type","ENSG00000239552", "ENSG00000224950", 
#  "ENSG00000250360", "ENSG00000253187", "ENSG00000254635", "ENSG00000255020", "ENSG00000256482", "ENSG00000257261")))])
#df$patient = rownames(df)
#df_c = df[,c("type", "patient")]
#df = df[,-which(colnames(df) %in% c("type", "patient"))]
#df = t(df)
#z = which(str_detect(rownames(df) ,"lgg"))
#df[z,] = apply(df[z,], 2, scale)
#z = which(str_detect(rownames(df) ,"gbm"))
#df[z,] = 0
#df = cbind(df, df_c)
#df$patient = NULL

#ha = HeatmapAnnotation(df = df , col = list(type = c("lgg" = "orange", "gbm" = "purple"), 
#  ENSG00000254635 = colorRamp2(c(-10, 0, 5), c("chartreuse4", "white", "tomato4")),
#    ENSG00000253187 = colorRamp2(c(-10, 0, 5), c("chartreuse4", "white", "tomato4")),
#ENSG00000256482 = colorRamp2(c(-10, 0, 5), c("chartreuse4", "white", "tomato4")),
#ENSG00000224950 = colorRamp2(c(-10, 0, 5), c("chartreuse4", "white", "tomato4")),
#ENSG00000255020 = colorRamp2(c(-10, 0, 5), c("chartreuse4", "white", "tomato4")),
#ENSG00000257261 = colorRamp2(c(-10, 0, 5), c("chartreuse4", "white", "tomato4")),
#ENSG00000239552 = colorRamp2(c(-10, 0, 5), c("chartreuse4", "white", "tomato4")),
#ENSG00000250360 = colorRamp2(c(-10, 0, 5), c("chartreuse4", "white", "tomato4"))))

#make sure order of patients in df matches order of patients in amtrix
df = df[,c("patient", "risk", "type", "IDH.status")]
df$patient = NULL
identical(rownames(df), colnames(mat))

order = colnames(mat)
library(ComplexHeatmap)
library(circlize)
values = df$risk
df$risk = NULL
typevals = df$type
idh = as.character(df$IDH.status)

ha = HeatmapAnnotation(points = anno_points(values, gp = gpar(size=0.6), axis = TRUE),
  type = typevals, idh_m = idh,
  col = list(type = c("lgg" = "orange", "gbm" = "purple"), idh_m = c("Mutant" = "black", "WT" = "white")),
    show_annotation_name = TRUE,
    annotation_name_offset = unit(2, "mm"),
    annotation_name_rot = c(0, 0, 90))

pdf("developmental_genes_lgg_gbm_all_small_version_genes_new.pdf", width=12, height=10)
Heatmap(mat, column_names_gp = gpar(fontsize = 1), top_annotation = ha, show_column_names = FALSE,
  heatmap_legend_param = list(legend_height = unit(3, "cm"), legend_width = unit(3, "cm")),
  top_annotation_height = unit(3, "cm"), clustering_distance_rows = "pearson")
dev.off()






