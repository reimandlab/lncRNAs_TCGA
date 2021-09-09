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

census = read.csv("Census_allFri_Jul_13_16_55_59_2018.csv")
  #get ensg
  get_census_ensg = function(genes){
  glist = unlist(strsplit(genes, ","))
  z = which(str_detect(glist, "ENSG"))
  ensg = glist[z]
  return(ensg)
  }
census$ensg = sapply(census$Synonyms, get_census_ensg)

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

#which diff exp genes in CGC?
all_de_results = as.data.table(filter(all_de_results, adj.P.Val < 0.05, abs(logFC) >=1))
z = which(census$Gene.Symbol %in% all_de_results$gene_name)
length(unique(census$Gene.Symbol[z]))

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

#all_lnc_pathways = llply(cancers, make_matrix_for_ap, .progress="text")
#all_lnc_pathways_df = ldply(all_lnc_pathways)
#done gave Marta matriced produced by this code for ActivePathways 


###---------------Summary figure-------------------------------###

#saveRDS(all_lnc_pathways_df, file="pathways_for_each_lncRNA_Oct30.rds")

all_lnc_pathways_df = readRDS("pathways_for_each_lncRNA_Oct30.rds")

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
z = which(!(duplicated(all_res$combo2)))
add = all_res[z,]
add$type = "pathways"
add$num_sig_des = 0
all_res = rbind(all_res, add)

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
canc = subset(res, res$term.name %in% canc_paths_paths)
canc_paths_paths = res[which(str_detect(res$term.name, "brain")),]$term.name

#1. get all PCGs that are in these pathways 
pcgs_brain = unique(unlist(canc$overlap)) #185 unique PCGs
lncs_brain = unique(unlist(canc$evidence)) #8 unique lncRNAs to be used as covariates high vs low 

##2-----------------label patients by risk------------------------------

dat = subset(all, Cancer=="Brain Lower Grade Glioma")

  dat_keep = dat[,which(colnames(dat) %in% c("patient", lncs_brain, pcgs_brain))]
  rownames(dat_keep) = dat_keep$patient
  dat_keep$patient = NULL
  #figure out which patients are high risk and which patients low risk
  for(i in 1:length(lncs_brain)){	

  	   z = which(colnames(dat_keep) %in% lncs_brain[i])
  	   median2 <- quantile(as.numeric(dat_keep[,z]), 0.5)

       if(median2 ==0){
        #if median = 0 then anyone greater than zero is 1 
        l1 = which(dat_keep[,z] > 0)
        l2 = which(dat_keep[,z] ==0)
        dat_keep[l1,z] = 1
        dat_keep[l2,z] = 0
        }

       if(!(median2 ==0)){
        l1 = which(dat_keep[,z] >= median2)
        l2 = which(dat_keep[,z] < median2)
        dat_keep[l1,z] = 1
        dat_keep[l2,z] = 0
        }
  }

  #subset gene expression to those pcgs
  #label patients by either high/low lncRNA expression 
  z = which(colnames(dat_keep) %in% c(pcgs_brain, "patient"))
  #heatmap 
  heat = dat_keep[,z]
  heat = log1p(heat)

  library(ComplexHeatmap)

  #tags <- dat_keep$risk
  #color.map <- function(tags) { if (tags=="RISK") "#FF0000" else "#0000FF" }
  #patientcolors <- unlist(lapply(tags, color.map))

  df2 = dat_keep[,c(2,7)] #lncRNA data

  col = list(
    ENSG00000254635 = c("1" = "red", "0" = "blue"),
    ENSG00000253187 = c("1" = "red", "0" = "blue"),
    ENSG00000256482 = c("1" = "red", "0" = "blue"),
    ENSG00000224950 = c("1" = "red", "0" = "blue"),
    ENSG00000255020 = c("1" = "red", "0" = "blue"),
    ENSG00000257261 = c("1" = "red", "0" = "blue"),
    ENSG00000239552 = c("1" = "red", "0" = "blue"),
    ENSG00000250360 = c("1" = "red", "0" = "blue"))

  # Create the heatmap annotation
  #ha <- HeatmapAnnotation(df2, col = col, annotation_height = unit.c(unit(0.25, "cm"), unit(0.25, "cm"), unit(0.25, "cm"), unit(0.25, "cm"),
  #	unit(0.25, "cm"), unit(0.25, "cm")))

  ha <- HeatmapAnnotation(df2, col = col, annotation_height = unit.c(unit(2, "cm"), unit(2, "cm")))

  # cluster on correlation
  heat = scale(heat)
  heat = t(heat)

  #change pcg names
  for(i in 1:nrow(heat)){
    pcg = rownames(heat)[i]
    newname = ucsc$hg19.ensemblToGeneName.value[which(ucsc$hg19.ensGene.name2 == pcg)][1]
    rownames(heat)[i] = newname
  }

  library(circlize)
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

  hox_genes = which(str_detect(rownames(heat), "HOX"))
  random_genes = c()#which(rownames(heat) %in% c("G6PC", "ALDH5A1", "IGF1", "PDK4", "SLC2A2", "DHDH", "IGF2", "G6PD"))
  z = which(rownames(heat) %in% census$Gene.Symbol)
  rownames(heat)[z] = paste(rownames(heat)[z], "*")

  subset = unique(c(hox_genes, random_genes, z))
  labels = rownames(heat)[subset]
  #labels = rownames(heat)


  #geom_tile plot 
  #show high vs low exp for each lncRNA
  #row --> patients 
  #column --> lncRNA high vs low 
  df2 = as.data.frame(df2)
  df2$patient = rownames(df2)
  df2 = as.data.table(melt(df2)) 
  df2$value = factor(df2$value, levels=c(0,1))
  df2 = df2[order(value, variable)]
  df2$variable = factor(df2$variable, levels=unique(df2$variable))

  pdf("lgg_brain_development_eight_lncRNAs.pdf", height=4, width=8)
  ggplot(df2, aes(variable, patient)) + geom_tile(aes(fill = value)) +
  #theme_void() + theme(legend.position="none")+
  scale_fill_manual(values=c("blue", "orange")) + coord_flip()
  dev.off()


  #pdf("heatmap_brain_development_pathway_LGG_sep19.pdf", width=12, height=4)
  pdf("heatmap_liver_brain_development_oct31.pdf", width=11, height=20)

  Heatmap(heat, clustering_distance_columns = "pearson", show_row_names = FALSE, 
  clustering_distance_rows = "pearson", col = colorRamp2(c(-3, 0, 3), c("steelblue1", "white", "orange")), 
  cluster_rows = TRUE, cluster_columns = TRUE, 
  clustering_method_rows = "complete", top_annotation = ha, 
  clustering_method_columns = "complete",show_column_names = FALSE, row_names_gp = gpar(fontsize = 1))+
  rowAnnotation(link = row_anno_link(at = subset, labels = labels),
  width = unit(0.75, "cm") + max_text_width(labels))
	
  dev.off()

#try another plot 

z = which(colnames(dat_keep) %in% c(pcgs_brain, "patient", lncs_brain))

#heatmap 
heat = dat_keep[,z]
heat = log1p(heat)
heat = scale(heat)
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

hox_genes = which(str_detect(rownames(heat), "HOX"))
random_genes = c()#which(rownames(heat) %in% c("G6PC", "ALDH5A1", "IGF1", "PDK4", "SLC2A2", "DHDH", "IGF2", "G6PD"))
z = which(rownames(heat) %in% census$Gene.Symbol)
rownames(heat)[z] = paste(rownames(heat)[z], "*")


pdf("brain_development_pathways.pdf", width=11, height=7)
res.dist <- get_dist(heat, stand = TRUE, method = "pearson")
fviz_dist(res.dist, 
   gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
dev.off()


#try corrplot
z = which(colnames(dat_keep) %in% c(pcgs_brain, "patient", lncs_brain))
#heatmap 
heat = dat_keep[,z]
heat = log1p(heat)
#heat = scale(heat)
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

hox_genes = colnames(heat)[which(str_detect(colnames(heat), "HOX"))]
random_genes = c()#which(rownames(heat) %in% c("G6PC", "ALDH5A1", "IGF1", "PDK4", "SLC2A2", "DHDH", "IGF2", "G6PD"))
z1 = which(colnames(heat) %in% c(census$Gene.Symbol))
z2 = which(colnames(heat) %in% ucsc$hg19.ensemblToGeneName.value[which(ucsc$hg19.ensGene.name2 %in% lncs_brain)])
colnames(heat)[z1] = paste(colnames(heat)[z1], "*")
z3 = which(colnames(heat) %in% hox_genes)
heat=heat[,unique(c(z1,z2, z3))]

M<-cor(heat)
library(corrplot)
res1 <- cor.mtest(heat, conf.level = .95)
corrplot(M, method="circle")
dev.off()
col<- colorRampPalette(c("blue", "white", "red"))(20)
pdf("brain_development_pathway.pdf", width=8, height=8)
corrplot(M, type="upper", order="hclust", col=col, tl.col="black", tl.srt=90, tl.cex=0.9, p.mat = res1$p, insig = "blank")
dev.off()


#get expression of these genes and compare between low lGG, high LGG and GBM 

#1. get list of genes 

#2. get LGG 

#3. get GBM 











#get PCA
pca_dat = heat

res.pca = prcomp(pca_dat, scale=TRUE)
pca_dat = as.data.frame(pca_dat)
pca_dat$lnc_combo = ""

get_mean = function(row){
  means = sum(as.numeric(row[[1]], row[[2]]))
  return(means)
}
pca_dat$lnc_combo = apply(pca_dat, 1, get_mean)
pca_dat = as.data.table(pca_dat)
pca_dat = pca_dat[order(lnc_combo)]

groups <- as.factor(dat$risk)

print(fviz_pca_ind(res.pca,
             #col.ind = groups, # color by groups
             #palette = c("#00AFBB",  "#FC4E07"),
             addEllipses = TRUE, # Concentration ellipses
             #legend.title = "Groups",
             repel = TRUE, label="none", ellipse.level=0.5
             ) #+ labs(title = paste("PCA", lnc, canc)))










