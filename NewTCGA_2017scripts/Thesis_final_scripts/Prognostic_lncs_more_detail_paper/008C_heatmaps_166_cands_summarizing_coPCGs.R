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

val_cands = read.csv("175_lncRNA_cancers_combos_23_cancer_types_july5.csv")
val_cands = as.data.table(val_cands)
val_cands = subset(val_cands, data == "PCAWG") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 
val_cands$combo = unique(paste(val_cands$gene, val_cands$cancer, sep="_"))
val_cands = subset(val_cands, top_pcawg_val == "YES") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 

#Combined into one dataframe because need to get ranks 
all <- merge(rna, pcg, by = c("patient", "Cancer"))
all = all[,1:25170]

#--------This script ------------------------------------------------

#summarize results from co-expression analysis of PCGs
#how many per lcnRNA
#how many pathways per lncRNA
#how many cancer genes, RBPs, TFs...

#--------------------------------------------------------------------

#write function that adds tag to whole data group 
#and does survival analysis on whole group

cands_dups = unique(allCands$gene[which(duplicated(allCands$gene))])

#--------------------------------------------------------------------
#RESULTS-------------------------------------------------------------
#--------------------------------------------------------------------

coexp = readRDS("coexpression_results_processed_july24.rds")

#all these have at least 10 sig DE PCGs
all_de_results = readRDS("coexpression_results_processed_july24.rds")
all_de_results = as.data.table(all_de_results)

#-------------------ANALYSIS--------------------------------------------
#Generate heatmap using those PCGs sig up/down regulated in lncRNA 
#risk or non-risk groups

#For each cancer type get all required data 
#PCG and lncRNA expression
combos = unique(all_de_results$combo)
cancs = sapply(combos, function(x){unlist(strsplit(x, "_"))[2]})
cancs = unique(cancs)

##1-----------------all expression--------------------------------------

get_tissue_specific <- function(combo){
  canc = unlist(strsplit(combo, "_"))[2]
  lnc = unlist(strsplit(combo, "_"))[1]
  tis = all[all$Cancer==canc,]
  tis$combo = combo
  print(combo)
  return(tis)
}
tissues_data <- llply(combos, get_tissue_specific, .progress="text")

##2-----------------label patients by risk------------------------------

get_lnc_canc = function(dat){
  cancer = dat$Cancer[1]
  combo = dat$combo[1]
  lnc = unlist(strsplit(combo, "_"))[1]

  pcgs = colnames(pcg)[2:19351]
  #keep only pcgs that are selected to be in lncRNA signature 
  z = which(all_de_results$combo == combo)
  lnc_pcgs = unique(all_de_results$ID[z])

  dat_keep = dat[,which(colnames(dat) %in% c("patient", lnc, lnc_pcgs))]
  rownames(dat_keep) = dat_keep$patient
  dat_keep$patient = NULL
  #figure out which patients are high risk and which patients low risk
  dat_keep$median <- ""
  median2 <- quantile(as.numeric(dat_keep[,1]), 0.5)

       if(median2 ==0){
        #if median = 0 then anyone greater than zero is 1 
        l1 = which(dat_keep[,1] > 0)
        l2 = which(dat_keep[,1] ==0)
        dat_keep$median[l1] = 1
        dat_keep$median[l2] = 0
        }

      if(!(median2 ==0)){
        l1 = which(dat_keep[,1] >= median2)
        l2 = which(dat_keep[,1] < median2)
        dat_keep$median[l1] = 1
        dat_keep$median[l2] = 0
    }

    #which one is high risk --> need surivval data
    dat_keep$patient = rownames(dat_keep)
    
      dat_keep$median[dat_keep$median ==0] = "Low"
      dat_keep$median[dat_keep$median==1] = "High"

      #cox ph
      z = which((allCands$gene == lnc) & (allCands$cancer == cancer))

      HR = as.numeric(allCands$HR[z])
      
      if(HR <1){
        risk = "Low"
        dat_keep$risk = ""
        dat_keep$risk[dat_keep$median=="High"] ="noRISK"
        dat_keep$risk[dat_keep$median=="Low"] ="RISK"
      }
      if(HR >1){
        risk = "High"
        dat_keep$risk = ""
        dat_keep$risk[dat_keep$median=="High"] ="RISK"
        dat_keep$risk[dat_keep$median=="Low"] ="noRISK"
      }

      dat_keep$lnc = colnames(dat_keep)[1]
      dat_keep$canc = cancer
      colnames(dat_keep)[1] = "lncRNA"

  return(dat_keep)  

}

#all lncRNAs with status 
all_canc_lnc_data = llply(tissues_data, get_lnc_canc, .progress="text")

##3-----------------generate heatmaps-----------------------------------

library(ComplexHeatmap)
library(circlize)

gen_heatmap = function(dat){
  dat$patient = rownames(dat)
  #lnc 
  lnc = dat$lnc[1]
  #canc
  canc = dat$canc[1]
  lnc_canc_combo = paste(lnc, canc, sep="_")
  #get which pcgs enriched in group
  z = which(all_de_results$combo %in% lnc_canc_combo)
 
  #cancer type
  lnc_pcgs = all_de_results[z,]
  print(length(unique(lnc_pcgs$ID))) 
  print(length(unique(lnc_pcgs$ID[lnc_pcgs$pcg_risk == "upregulated_in_risk"])))
  print(length(unique(lnc_pcgs$ID[lnc_pcgs$pcg_risk == "downregulated_in_risk"])))

  #subset gene expression to those pcgs
  #label patients by either high/low lncRNA expression 
  z = which(colnames(dat) %in% c(lnc_pcgs$ID, "patient"))
  #heatmap 
  heat = dat[,z]
  rownames(heat) = heat$patient
  heat$patient = NULL

  heat = log1p(heat)

  if(dim(heat)[2] > 75){
    #get most variable genes and only include them in heatmap 
    #instead of variance order them by absolute fold change!? see if that looks better?
    #vars = data.table(genes = colnames(heat), var = apply(heat, 2, var))
    #vars = vars[order(-var)]
    #vars = vars[1:50,]
    #z = which(colnames(heat) %in% vars$genes)
    #heat = heat[,z]
    z = which(all_de_results$combo %in% lnc_canc_combo)
    pcg_keep = as.data.table(all_de_results[z,])
    pcg_keep = as.data.table(filter(pcg_keep, adj.P.Val <= 0.05))
    pcg_keep = pcg_keep[order(-abs(logFC))]
    pcgs = pcg_keep[1:75,1]
    z = which(colnames(heat) %in% pcgs$ID)
    heat = heat[,z]
  }

  tags <- dat$risk
  color.map <- function(tags) { if (tags=="RISK") "#FF0000" else "#0000FF" }
  patientcolors <- unlist(lapply(tags, color.map))

  #label whether pcg is fav or unfav 
  z = which(all_de_results$combo %in% lnc_canc_combo)
  lnc_pcgs = all_de_results[z,]

  # cluster on correlation
  heat = scale(heat)
  heat = t(heat)

  pcg_order = as.data.frame(matrix(ncol =2))
  for(i in 1:nrow(heat)){
    pcg = rownames(heat)[i]
    type = lnc_pcgs$pcg_risk[which(lnc_pcgs$ID == pcg)]
    row = c(pcg, type)
    pcg_order = rbind(pcg_order, row)
  }
  pcg_order = pcg_order[-1,]
  colnames(pcg_order) = c("pcg", "type")

  #change pcg names
  for(i in 1:nrow(heat)){
    pcg = rownames(heat)[i]
    newname = ucsc$hg19.ensemblToGeneName.value[which(ucsc$hg19.ensGene.name2 == pcg)][1]
    rownames(heat)[i] = newname
  }

  pcgs <- pcg_order$type
  color.map <- function(pcgs) { if (pcgs=="Risk") "Purple" else "Yellow" }
  pcg_cols <- unlist(lapply(pcgs, color.map))

  hc <- hclust(as.dist(1 - cor(t(heat))), method="ward.D2")
  # draw a heatmap
  my_palette <- colorRampPalette(c("blue", "white", "orange"))(n = 100)
  
  lnc_clean = rna$type[which(rna$Cancer == canc)][1]
  canc_clean = fantom$CAT_geneName[which(fantom$CAT_geneID == lnc)]
  risk_lnc = dat$median[which(dat$risk == "RISK")][1]

  title = paste(lnc_clean, canc_clean, "\nRISK = ", risk_lnc, "Expression")    
  ha_column = HeatmapAnnotation(df = data.frame(risk = tags),
    col = list(risk = c("RISK" =  "#FF0000", "noRISK" = "#0000FF")))  

  ha_row = rowAnnotation(df = data.frame(PCGtype = pcgs),
    col = list(PCGtype = c("upregulated_in_risk" = "Yellow", "downregulated_in_risk" = "Purple")))

  h1 = Heatmap(heat, col = colorRamp2(c(-3, 0, 3), c("steelblue1", "white", "orange")), clustering_distance_columns = "pearson", column_title=title, 
    clustering_distance_rows = "spearman", cluster_rows = TRUE, cluster_columns = TRUE, top_annotation = ha_column, show_column_names = FALSE, 
    clustering_method_rows = "complete", clustering_method_columns = "complete", row_names_gp = gpar(fontsize = 2))
  h1 + ha_row

  #heatmap.2(as.matrix(heat), col=my_palette, ColSideColors= patientcolors, cexRow=0.5, cexCol=0.6, Rowv=as.dendrogram(hc), 
  #  RowSideColors= pcg_cols, trace="none", scale="row", dendrogram="row", labCol="", main = title, key=FALSE)

}

pdf("lncs_wSIG_PCGs_heatmaps_sep19_top75_genes.pdf", width=10, height=10)
llply(all_canc_lnc_data, gen_heatmap, .progress="text")
dev.off()

##----DONE SIR----











