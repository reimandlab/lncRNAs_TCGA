library(survAUC)
#source("source_code_Cox_MonteCarlo_CV_Mar13.R")
require(caTools)
#check if this person is in my analysis: TCGA-61-2095
library(glmnet)
library(survcomp)
library(caret)
library(stringr)

source("universal_LASSO_survival_script.R")

library(ggpubr)
library(ggrepel)
library(viridis)
library(patchwork)
library(caret)  
library(Rtsne)
library("FactoMineR")
library(ComplexHeatmap)
library(circlize)
library(reshape2)

#------FEATURES-----------------------------------------------------

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = filter(allCands, data=="TCGA") #168 unique lncRNA-cancer combos, #166 unique lncRNAs 

colnames(fantom)[1] = "gene"

canc_conv = rna[,which(colnames(rna) %in% c("Cancer", "type"))]
canc_conv = canc_conv[!duplicated(canc_conv), ]

mypal = c("#E5DFD9","#EAD286" ,"#D1EB7B", "#96897F" ,"#E5C0A6" ,
  "#72A93B", "#74DAE3" ,"#49B98D" ,"#D97B8F" ,"#70A2A4", "#64709B" ,"#DFBF38" ,"#61EA4F" ,
  "#C7CBE7", "#786DDA",
"#CFA0E0" ,"#67E9D0" ,"#7C9BE1", "#D94753" ,
"#AAE6B0", "#D13BDF" ,"#DEAEC7" ,"#BBE6DF" ,"#B2B47A" ,"#E6ECBA", "#C86ED7",
 "#7BEE95" ,"#6F46E6" ,"#65B9E0", "#C0EC3E",
"#DE8D54" ,"#DF4FA6")

#remove cancer types with less than 50 patients 
pats_num = as.data.table(table(rna$Cancer))
pats_num = filter(pats_num, N <50)
canc_rm = pats_num$V1

#remove those ones
cancers = unlist(cancers[which(!(cancers %in% canc_rm))])
rna = subset(rna, Cancer %in% cancers)

#-------------------------------------------------------------------
#-----------------PCA using just all expressed lncRNA --------------
#-------------------------------------------------------------------

table(rna$type)
rna= as.data.frame(rna)
#seperate into labels and data used for clustering 
z = which(str_detect(colnames(rna), "ENSG"))
rna_data = rna[,z]
rna_labels = rna[,-z]
rna_labels = rna_labels[,"type"]

rna_data = log1p(rna_data)

library(umap)
umap = umap(rna_data)



function(x, labels,
       main="A UMAP visualization of the Iris dataset",
          pad=0.1, cex=0.65, pch=19, add=FALSE, legend.suffix="",
          cex.main=1, cex.legend=1) { 
   layout = x
  if (class(x)=="umap") {
    layout = x$layout
  } 
  
  xylim = range(layout)
  xylim = xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)
  if (!add) {
     par(mar=c(0.2,0.7,1.2,0.7), ps=10)
     plot(xylim, xylim, type="n", axes=F, frame=F)
     rect(xylim[1], xylim[1], xylim[2], xylim[2], border="#aaaaaa", lwd=0.25)  
   }
   points(layout[,1], layout[,2], col=iris.colors[as.integer(labels)],
          cex=cex, pch=pch)
   mtext(side=3, main, cex=cex.main)
 
   labels.u = unique(labels)
   legend.pos = "topright"
   legend.text = as.character(labels.u)
   if (add) {
     legend.pos = "bottomright"
     legend.text = paste(as.character(labels.u), legend.suffix)
   }
   legend(legend.pos, legend=legend.text,
          col=iris.colors[as.integer(labels.u)],
          bty="n", pch=pch, cex=cex.legend) 
}



















#1. get 1000 most variable lncRNAs across cancer types 

z = which(str_detect(colnames(rna), "ENSG"))
vars = as.data.frame(apply(rna[,z], 2, var))
vars$gene = rownames(vars)
vars = as.data.table(vars)
colnames(vars)[1] = "variance" 
vars = vars[order(-variance)]

vars_lncs = vars$gene[1:500]

z = which(colnames(rna) %in% c(vars_lncs, "type"))
vars_rna= rna[,z]

z = which(str_detect(colnames(vars_rna), "ENSG"))
vars_rna[,z] = log1p(vars_rna[,z])

rownames(vars_rna) = paste(vars_rna$type, 1:length(vars_rna$type), sep="")
df = as.data.frame(vars_rna$type)
df$id = rownames(vars_rna)

vars_rna = vars_rna[,which(str_detect(colnames(vars_rna), "ENSG"))]
xmat = scale(vars_rna)
xmat=t(xmat)

colnames(df)[1] = "cancer"
identical(df$id, colnames(xmat))

order = colnames(xmat)
library(ComplexHeatmap)
library(circlize)

cancers = as.character(df$cancer)

ha = HeatmapAnnotation(
  m = cancers)

pdf("testing_clustering_patients_lncRNAs_top500.pdf", width=20)
Heatmap(xmat, column_names_gp = gpar(fontsize = 1), top_annotation = ha, show_column_names = FALSE, show_row_names = FALSE,
  heatmap_legend_param = list(legend_height = unit(3, "cm"), legend_width = unit(3, "cm")),
  top_annotation_height = unit(1, "cm"), clustering_distance_rows = "pearson")
dev.off()


#get canc data
cancers = unique(vars_rna$type)
get_canc = function(canc){
	canc_dat = subset(vars_rna, type == canc)
	return(canc_dat)
}

canc_datas = llply(cancers, get_canc)

#get median lncRNA exp for each cancer type 

get_medians = function(dtt){
	
	z1 = which(str_detect(colnames(dtt), "ENSG"))	
	#meds 
	lncs = colnames(dtt)[z1]
	meds = apply(dtt[,z1], 2, median)
	meds = as.data.frame(meds)
	meds$canc = dtt$type[1]
	meds$gene = rownames(meds)
	meds = merge(meds, fantom, by="gene")
	return(meds)
}

meds_cancers = llply(canc_datas, get_medians, .progress="text")
meds_cancers1 = ldply(meds_cancers, data.frame)
meds_dat = meds_cancers1[,c("gene", "canc", "meds")]
meds_dat = acast(meds_dat, canc ~ gene)
meds_dat_scale = apply(meds_dat, 2, scale)
#meds_dat = t(meds_dat)
#rownames(meds_dat_scale) = rownames(meds_dat)

Heatmap(meds_dat, km = 10)
dev.off()





