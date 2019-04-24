library(survAUC)
require(caTools)
#check if this person is in my analysis: TCGA-61-2095
library(glmnet)
library(survcomp)
library(caret)
library(stringr)
library(ggpubr)
library(ggrepel)
library(viridis)
library(patchwork)
library(caret)  
library(Rtsne)
library(EnvStats)
library(reshape2)

source("check_lnc_exp_cancers.R")

#------FEATURES-----------------------------------------------------

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = subset(allCands, data == "TCGA") #173 unique lncRNA-cancer combos, #166 unique lncRNAs 
allCands$combo = unique(paste(allCands$gene, allCands$cancer, sep="_"))
cands_dups = unique(allCands$gene[which(duplicated(allCands$gene))])
colnames(allCands)[7] = "Cancer"
allCands = merge(allCands, canc_conv, by="Cancer")
allCands$combo = paste(allCands$gene, allCands$type)

#------PCGs for pathway analysis--------------------------------------

all_de_results = readRDS("diff_expressed_PCGs_lncRNA_risk_groups_Aug21.rds")
all_de_results = as.data.table(all_de_results)

z = which(all_de_results$cancer == "LGG")
all_de_results = all_de_results[-z,]

all_de_results_lgg = readRDS("diff_expressed_PCGs_lncRNA_risk_groups_lgg_nov30.rds")
all_de_results_lgg = as.data.table(all_de_results_lgg)
all_de_results = rbind(all_de_results, all_de_results_lgg)

all_de_results$combo = paste(all_de_results$lnc, all_de_results$cancer) #148/158 unique combos across 22 unique cancer types 
all_de_results = as.data.table(filter(all_de_results, combo %in% allCands$combo))

full_diff_exp = all_de_results

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
	
  #keep only those with logFC > 1 or logFC < 1
  z1 = which(dat_all$logFC >=1)
  z2 = which(dat_all$logFC <=-1)

  dat_all = dat_all[c(z1,z2),]

  dat_all_matrix = acast(dat_all, ID~lnc, value.var="P.Value")
	dat_all_matrix[is.na(dat_all_matrix)] = 1

	file = paste("Aug22_DE_genes_fActivePathways/", canc, "all_up_down_genes_fActivepathways_March25_updFC.rds", sep="_")
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

saveRDS(all_lnc_pathways_df, file="pathways_for_each_lncRNA_Oct30.rds")
all_lnc_pathways_df = readRDS("pathways_for_each_lncRNA_Oct30.rds")

#number of PCGs/lncRNA vs number of Pathways/lncRNAs 

#1. DE PCGs/lncRNA
t = as.data.table(filter(all_de_results, adj.P.Val <= 0.05))
#keep only those with fc >1/<-1
t = as.data.table(filter(t, (logFC >= 1) | (logFC < -1)))
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

sig_des = as.data.table(filter(sig_des, combo %in% sig_paths_sum$combo))
all_lnc_pathways_df = as.data.table(filter(all_lnc_pathways_df, combo %in% sig_paths_sum$combo))

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
census_genes = as.data.table(filter(sig_des, ID %in% census$ensg))


pdf("summary_#DE_pcgs_vs_pathways_pathways_figure_jdan2519.pdf", width=9, height=5)
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

write.csv(sig_paths_sum, file="summary_num_diff_exp_genes_pathways_jan2519.csv", quote=F, row.names=F)

#network part A 
pdf("summary_barplot_#DE_pcgs_vs_pathways_pathways_figure_jan2519.pdf", width=16, height=5)
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
a = ifelse(all_res$font_col == 1, "red", "black")
a <- ifelse(sig_paths_sum$font_col == 1, "red", "black")

#all_res = as.data.table(filter(all_res, num_sig_des >=5))

#117 unique lncRNA-cancer pairs 
#114 unique lncRNAs 

z = which(all_res$gene %in% census_genes$lnc)
census_cands = all_res[z,]
census_cands$type = "census_genes"
census_cands$num_sig_des = ""
for(i in 1:length(census_cands$combo)){
  genes = subset(census_genes, combo %in% census_cands$combo[i])
  num = dim(genes)[1]
  census_cands$num_sig_des[i] = num
}

all_res = rbind(census_cands, all_res)
all_res$combo3 = paste(all_res$combo, all_res$type)
z = which(duplicated(all_res$combo3))
if(!(length(z)==0)){
all_res = all_res[-z,]}

#all_res$num_sig_des = log2(all_res$num_sig_des)
all_res$type = factor(all_res$type, levels = c("genes", "census_genes", "pathways"))

order_l = as.data.table(all_res %>% group_by(combo2) %>% summarise(frequency =  sum(num_sig_des)))
order_l = order_l[order(-frequency)]
all_res$combo2 = factor(all_res$combo2, levels = order_l$combo2)

pdf("summary_barplot_genenames_#DE_pcgs_vs_pathways_pathways_figure_jan2519.pdf", width=8, height=4)
g = ggplot(all_res, aes(combo2, num_sig_des, group=type)) + theme_classic() + 
   geom_col(position = position_stack(), aes(fill=type), color="grey29")+xlab("lncRNA-cancer pair") + ylab("Number of sig genes/pathways")+
   theme(legend.title=element_blank(), legend.position="top", axis.title.x=element_blank(), 
   	axis.text.x = element_text(angle = 90, hjust = 1, size=5)) + scale_fill_manual(values=c("#E69F00", "red", "#56B4E9"))+
   scale_y_continuous(breaks=seq(0,2000,250))
#ggpar(g,
# font.tickslab = c(6,"plain", "black"),
# xtickslab.rt = 45)
print(g)
dev.off()

#stopped here Jan 25

#make barplot just for LGG candidates found in the network
keep = c("ENSG00000253187", "ENSG00000239552", "ENSG00000224950", "ENSG00000254635", "ENSG00000250360", "ENSG00000256482", "ENSG00000257261")
lgg_res = as.data.table(filter(full_diff_exp, cancer == "LGG", adj.P.Val <= 0.05, (logFC >= 1) | (logFC < -1), lnc %in% keep)) 
lgg_res$type[lgg_res$logFC >= 1] = "UpregulatedRisk"
lgg_res$type[lgg_res$logFC < -1] = "DownregulatedRisk"

#get order
lgg_res$gene_name = sapply(lgg_res$lnc, get_name)
ttt = as.data.table(table(lgg_res$gene_name))
ttt = ttt[order(N)]

#how many cancer gene census genes
z = which(lgg_res$ID %in% census$ensg)
lgg_res$census[z] = "census"
lgg_res$census[-z] = "not_census"

barplot1 = as.data.table(table(lgg_res$type, lgg_res$gene_name))
barplot2 = as.data.table(table(lgg_res$census, lgg_res$gene_name))
barplot2 = as.data.table(filter(barplot2, V1 == "census"))
barplot = rbind(barplot1, barplot2)
colnames(barplot) = c("type", "lncRNA", "num")
barplot = as.data.table(filter(barplot, num >0))
barplot$lncRNA = factor(barplot$lncRNA, levels = ttt$V1)

pdf("figure6A_lgg_only.pdf", width=5, height=6)
g = ggplot(barplot, aes(x=lncRNA, y=num, fill=type)) + theme_classic() + 
   geom_bar(aes(fill=type), stat="identity", position=position_dodge())+xlab("LGG lncRNAs") + ylab("Number of genes")+
   theme(legend.title=element_blank(), legend.position="top", axis.title.x=element_blank(), 
    axis.text.x = element_text(angle = 90, hjust = 1, size=8)) + scale_fill_manual(values=c("#56B4E9","red", "#E69F00"))
ggpar(g,
 font.tickslab = c(9,"plain", "black"),
 xtickslab.rt = 45)
dev.off()

###---------------Make heatmap-------------------------------###

#liver gluco related pathways 

load("_LGG_all_up_down_genes_.2018-12-13.rdata")
colnames(res)
#canc = subset(res, res$term.name %in% canc_paths_paths)
canc_paths_paths = res[which(str_detect(res$term.name, "brain")),]$term.name
canc = subset(res, term.name %in% canc_paths_paths)

#1. get all PCGs that are in these pathways 
pcgs_brain = unique(unlist(canc$overlap)) #50 unique PCGs
lncs_brain = unique(unlist(canc$evidence))[1:2] #2 unique lncRNAs to be used as covariates high vs low 

evi = c()
for(i in 1:nrow(res)){
  ress = unlist(res$evidence[i])
  if(ress == "combined"){
    ress = c(1,2)
  }
  evi = c(evi, length(ress))
}


##2-----------------label patients by risk------------------------------

dat = subset(all, type=="LGG")
#lgg_nearby_genes = readRDS("lgg_nearby_genes.rds")

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

write.csv(res_all_brain, file="all_brain_pathways_lgg_march26.csv", quote=F, row.names=F)

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
#heat = heat[,unique(c(z1,z2,z3))]

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
#z = which(all$type == "GBM")
#gbm = all[z,]
gbm = subset(all, type == "GBM")
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

gb_idh = readRDS("gbm_clin_subtypes_glioblastoma.rds")
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


z = which(colnames(all_canc) %in% c("ENSG00000239552", "HOXA10-AS"))
old_canc = all_canc
all_canc = all_canc[,-z]

mat = all_canc[,c(which(!(colnames(all_canc) %in% c("type", "IDH.status"))))]
mat = scale(mat)
mat=t(mat)

df = as.data.frame(old_canc[,c(which(colnames(old_canc) %in% c("type", "IDH.status", "ENSG00000239552", "HOXA10-AS")))])
colnames(df)[2] = get_name("ENSG00000239552")
df$patient = rownames(df)

lnc_med = median(df[,1])
z1=which(df[,1] > lnc_med)
z2=which(df[,1] == lnc_med)
df[z1,1] = 1
df[z2,1] = 0

lnc_med = median(df[,2])
z1=which(df[,2] > lnc_med)
z2=which(df[,2] == lnc_med)
df[z1,2] = 1
df[z2,2] = 0

#which are favourable lncs? (for those ones flip 0s and 1s)
cands_lncs = as.data.table(filter(allCands, Cancer== "Brain Lower Grade Glioma"))
#cands_lncs = as.data.table(filter(cands_lncs, HR <1))
#cands_lncs = cands_lncs$gene

#z = which(colnames(df) %in% cands_lncs)
#for(i in 1:length(z)){
#  k = which(df[,z[i]] == 0)
#  df[k,z[i]] = 1
#  df[-k,z[i]] = 0
#}

df$risk = apply(df[,1:2], 1, sum)

df$patientt = sapply(df$patient, function(x){unlist(strsplit(x, " "))[1]})
surv = all[,c("patient", "OS", "OS.time")]
colnames(surv)[1] = "patientt"
df = merge(df, surv, by="patientt")
get_risk = as.data.table(filter(df, type=="lgg"))
gbm = as.data.table(filter(df, type=="gbm"))
colnames(get_risk)[2:3] = c("ENSG00000253187", "ENSG00000239552")

fitCPH <- coxph(Surv(OS.time, OS) ~ ENSG00000253187 + ENSG00000239552, data=get_risk)    # Cox-PH model
(coefCPH <- coef(fitCPH))    

meanENSG00000253187  <- mean(as.numeric(get_risk$ENSG00000253187) - 1)   # average of financial aid dummy
meanENSG00000239552  <- mean(as.numeric(get_risk$ENSG00000239552)-1)                  
rMean <- exp(coefCPH["ENSG00000253187"] * meanENSG00000253187 +       # e^Xb
           coefCPH["ENSG00000239552"]   *meanENSG00000239552)

all_lgg_pats <- exp(coefCPH["ENSG00000253187"]* (get_risk[1:4, "ENSG00000253187"]-1)
           + coefCPH["ENSG00000239552"] * (get_risk[1:4, "ENSG00000239552"]-1))


relRisk <- predict(fitCPH, get_risk, type="risk")   # relative risk
get_risk$risk = relRisk
colnames(get_risk)[2:3] = c("HOXA10-AS", "HOXB-AS2")

df = rbind(get_risk, gbm)

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
df = df[,2:7]
rownames(df) = df$patient
#df$patient = NULL

identical(rownames(df), colnames(mat))

df$patient <- factor(df$patient, levels=colnames(mat))
df = df[order(df$patient),]
df = as.data.frame(df)
rownames(df) = df$patient
df$patient = NULL
identical(rownames(df), colnames(mat))

order = colnames(mat)
library(ComplexHeatmap)
library(circlize)
values = df$risk
df$risk = NULL
typevals = df$type
idh = as.character(df$IDH.status)
hoxa10as = df[,1]
hoxbas2 = df[,2]

ha = HeatmapAnnotation(relative_risk = anno_points(values, gp = gpar(size=0.3), axis = TRUE),
  type = typevals, m = idh, hoxa10as = hoxa10as, hoxbas2=hoxbas2, 
  col = list(type = c("lgg" = "orange", "gbm" = "purple"), m = c("Mutant" = "black", "WT" = "gray88"),
  hoxa10as = c("1" = "limegreen", "0" = "mistyrose2"), hoxbas2=c("1" = "limegreen", "0" = "mistyrose2")), 
    show_annotation_name = TRUE,
    annotation_name_offset = unit(3, "mm"))

#pdf("developmental_genes_lgg_gbm_all_small_version_genes_new_march26.pdf", width=12, height=4)
pdf("developmental_genes_lgg_gbm_all_small_version_genes_new_april11.pdf", width=9, height=7)
Heatmap(mat, column_names_gp = gpar(fontsize = 1), top_annotation = ha, show_column_names = FALSE,
  heatmap_legend_param = list(legend_height = unit(2, "cm"), legend_width = unit(2, "cm")),
  top_annotation_height = unit(3, "cm"), row_names_gp = gpar(fontsize = 5), clustering_distance_rows = "pearson", clustering_distance_columns = "pearson")
dev.off()


#tiff("developmental_genes_lgg_gbm_all_small_version_genes_new_march26.tiff", height = 10, width = 25, units= "cm",  
#     compression = "lzw", res = 300)
#Heatmap(mat, column_names_gp = gpar(fontsize = 1), top_annotation = ha, show_column_names = FALSE,
#  heatmap_legend_param = list(legend_height = unit(3, "cm"), legend_width = unit(3, "cm")),
#  top_annotation_height = unit(2, "cm"), clustering_distance_rows = "pearson", row_names_gp = gpar(fontsize = 5))
#dev.off()


#get km plots for all the PCGs
genes = rownames(mat)[3:nrow(mat)]
genes = sapply(genes, get_ensg_pcg)

#pdf("LGG_heatmap_KM_plots.pdf")
#list = sapply(genes,get_km_plot, cancer = "LGG")
#dev.off()

#genes = c("ENSG00000253187", "ENSG00000239552")
#pdf("LGG_heatmaps_KM_plots_lncrnas.pdf")
#list = sapply(genes,get_km_plot, cancer = "LGG")
#dev.off()

library(gridExtra)
#n <- length(list)
#nCol <- floor(sqrt(n))

#pdf("hoxlncs.pdf", width=8, height=11)
#do.call("grid.arrange", c(list, ncol=nCol))
#dev.off()

dim(lgg_idh)
head(lgg_idh)
lgg = as.data.table(subset(rna, type == "LGG"))
lgg = lgg[,c("ENSG00000253187", "ENSG00000239552", "patient", "type", "OS", "OS.time")]
lgg = merge(lgg, lgg_idh, by="patient")

lgg$ENSG00000253187[lgg$ENSG00000253187 > 0] = "High"
lgg$ENSG00000253187[lgg$ENSG00000253187 == 0] = "Low"

lgg$ENSG00000239552[lgg$ENSG00000239552 > 0] = "High"
lgg$ENSG00000239552[lgg$ENSG00000239552 == 0] = "Low"

colnames(lgg)[3] = "HOXB-AS2"
colnames(lgg)[2] = "HOXA10-AS"

#get km plots

lgg$OS.time = lgg$OS.time/365

gene_name = colnames(lgg)[2]
colnames(lgg)[2] = "gene"

lgg$IDH.status = factor(lgg$IDH.status, levels=c("WT", "Mutant"))

fit <- survfit(Surv(OS.time, OS) ~ gene + IDH.status, data = lgg)

pdf("lgg_two_lncRNAs_cands_figure6d.pdf", width=8, height=7)

s <- ggsurvplot(
          title = gene_name, 
          fit, 
          xlab = "Time (Years)", 
          #surv.median.line = "hv",
          font.main = c(14, "bold", "black"),
          font.x = c(12, "plain", "black"),
          font.y = c(12, "plain", "black"),
          font.tickslab = c(12, "plain", "black"),
          font.legend = 10,
          risk.table.fontsize = 5, 
          #legend.labs = c("High Expression", "Low Expression"),             # survfit object with calculated statistics.
          data = lgg,      # data used to fit survival curves. 
          risk.table = TRUE,       # show risk table.
          legend = "right", 
          pval = TRUE,             # show p-value of log-rank test.
          conf.int = FALSE,        # show confidence intervals for 
                            # point estimaes of survival curves.
          xlim = c(0,5),        # present narrower X axis, but not affect
                            # survival estimates.
          break.time.by = 1,     # break X axis in time intervals by 500.
          #palette = colorRampPalette(mypal)(14), 
          #palette = mypal[c(4,1)],
          palette = "npg", 
           legend.labs = c("High exp & IDH WT", "High exp & IDH Mutant", "Low exp & IDH WT", "Low exp & IDH Mutant"), 
          #ggtheme = theme_bw(), # customize plot and risk table with a theme.
          risk.table.y.text.col = T, # colour risk table text annotations.
          risk.table.y.text = FALSE # show bars instead of names in text annotations
                            # in legend of risk table
          )

print(s)

gene_name = colnames(lgg)[3]
colnames(lgg)[2] = "lnc"

colnames(lgg)[3] = "gene"

lgg$IDH.status = factor(lgg$IDH.status, levels=c("WT", "Mutant"))

fit <- survfit(Surv(OS.time, OS) ~ gene + IDH.status, data = lgg)

s <- ggsurvplot(
          title = gene_name, 
          fit, 
          xlab = "Time (Years)", 
          #surv.median.line = "hv",
          font.main = c(14, "bold", "black"),
          font.x = c(12, "plain", "black"),
          font.y = c(12, "plain", "black"),
          font.tickslab = c(12, "plain", "black"),
          font.legend = 10,
          risk.table.fontsize = 5, 
          #legend.labs = c("High Expression", "Low Expression"),             # survfit object with calculated statistics.
          data = lgg,      # data used to fit survival curves. 
          risk.table = TRUE,       # show risk table.
          legend = "right", 
          pval = TRUE,             # show p-value of log-rank test.
          conf.int = FALSE,        # show confidence intervals for 
                            # point estimaes of survival curves.
          xlim = c(0,5),        # present narrower X axis, but not affect
                            # survival estimates.
          break.time.by = 1,     # break X axis in time intervals by 500.
          #palette = colorRampPalette(mypal)(14), 
          #palette = mypal[c(4,1)],
          palette = "npg", 
           legend.labs = c("High exp & IDH WT", "High exp & IDH Mutant", "Low exp & IDH WT", "Low exp & IDH Mutant"), 
          #ggtheme = theme_bw(), # customize plot and risk table with a theme.
          risk.table.y.text.col = T, # colour risk table text annotations.
          risk.table.y.text = FALSE # show bars instead of names in text annotations
                            # in legend of risk table
          )
print(s)

dev.off()



lgg = as.data.table(subset(rna, type == "LGG"))
lgg = lgg[,c("ENSG00000253187", "ENSG00000239552", "patient", "type", "OS", "OS.time")]
lgg = merge(lgg, lgg_idh, by="patient")

z = which(is.na(lgg$IDH.status))
lgg = lgg[-z,]

lgg$ENSG00000253187 = log1p(lgg$ENSG00000253187)
lgg$ENSG00000253187_tag[lgg$ENSG00000253187 > 0] = "High"
lgg$ENSG00000253187_tag[lgg$ENSG00000253187 == 0] = "Low"
lgg$ENSG00000253187_tag = factor(lgg$ENSG00000253187_tag, levels=c("Low", "High"))

lgg$ENSG00000239552 = log1p(lgg$ENSG00000239552)
lgg$ENSG00000239552_tag[lgg$ENSG00000239552 > 0] = "High"
lgg$ENSG00000239552_tag[lgg$ENSG00000239552 == 0] = "Low"
lgg$ENSG00000239552_tag = factor(lgg$ENSG00000239552_tag, levels=c("Low", "High"))

lgg$IDH.status = factor(lgg$IDH.status, levels=c("WT", "Mutant"))
lgg$ENSG00000253187_tag = paste(lgg$ENSG00000253187_tag, lgg$IDH.status)
lgg$ENSG00000253187_tag = factor(lgg$ENSG00000253187_tag, levels=c("High WT", "High Mutant", "Low WT", "Low Mutant"))

lgg$ENSG00000239552_tag = paste(lgg$ENSG00000239552_tag, lgg$IDH.status)
lgg$ENSG00000239552_tag = factor(lgg$ENSG00000239552_tag, levels=c("High WT", "High Mutant", "Low WT", "Low Mutant"))

pdf("lgg_two_lncRNAs_cands_figure6d_boxplots.pdf", width=6, height=6)

#make boxplots
   p <- ggboxplot(lgg, x = "ENSG00000253187_tag", y = "ENSG00000253187",
          color = "ENSG00000253187_tag",
         palette = "npg", title = "HOXA10-AS expression", 
          add = "jitter", ylab = "log1p(FPKM-UQ)",  ggtheme = theme_classic())
        # Change method
  p = p + stat_n_text() + scale_color_npg() + theme(legend.position="none")+xlab("HOXA10-AS expression & IDH mutation status")
  ggpar(p, font.tickslab = c(12, "plain", "black")) 

#make boxplots
   p <- ggboxplot(lgg, x = "ENSG00000239552_tag", y = "ENSG00000239552",
          color = "ENSG00000239552_tag",
         palette = "npg", title = "HOXB-AS2 expression", 
          add = "jitter", ylab = "log1p(FPKM-UQ)",  ggtheme = theme_classic())
        # Change method
  p = p + stat_n_text() + scale_color_npg() + theme(legend.position="none") +xlab("HOXB-AS2 expression & IDH mutation status")
  ggpar(p, font.tickslab = c(12, "plain", "black")) 
dev.off()



















