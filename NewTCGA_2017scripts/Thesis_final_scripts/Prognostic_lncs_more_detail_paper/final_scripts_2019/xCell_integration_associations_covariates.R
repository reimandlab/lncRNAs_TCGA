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

source("check_lnc_exp_cancers.R")

#final script used to generate clinical analyiss data now figure 4

#------FEATURES-----------------------------------------------------
setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNAs_2019_manuscript")

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = subset(allCands, data == "TCGA") #173 unique lncRNA-cancer combos, #166 unique lncRNAs 
allCands$combo = unique(paste(allCands$gene, allCands$cancer, sep="_"))
cands_dups = unique(allCands$gene[which(duplicated(allCands$gene))])

library(TCGAbiolinks)

#xcell = fread("xCell_TCGA_RSEM.txt")
#clean up patient IDs
#pats=colnames(xcell)
#z1 = which(str_detect(colnames(xcell), "TCGA"))
#xcell=as.data.frame(xcell)
#xcell=xcell[,c(1, z1)]

check=function(pat){
  pat_new=paste(unlist(strsplit(pat, "\\."))[1:3], collapse="-")
  test=unlist(strsplit(pat, "\\."))[4]
  print(test)
  return(pat_new)
}

#colnames(xcell)[2:ncol(xcell)] =  llply(colnames(xcell)[2:ncol(xcell)], check)
#saveRDS(xcell, file="xcell_new_colnames.rds")

xcell=readRDS("xcell_new_colnames.rds")

#--------This script ------------------------------------------------

#--------------------------------------------------------------------
#Clinical files - use TCGAbiolinks
#--------------------------------------------------------------------

#write function that adds tag to whole data group 
#and does survival analysis on whole group

cancers = unique(allCands$canc)
get_canc_dat = function(canc){
  canc_d = subset(rna, Cancer == canc)
  return(canc_d)
}
cancer_data = llply(cancers, get_canc_dat)

get_canc_data_for_plot = function(dtt){
  #get cancer specific candidates 
  z = which(colnames(dtt) %in% c(as.character(allCands$gene[allCands$cancer == dtt$Cancer[1]]), "age_at_initial_pathologic_diagnosis", 
    "OS.time", "OS", "gender", "race", "patient", "clinical_stage", "histological_grade", "treatment_outcome_first_course", 
    "new_tumor_event_type", "Cancer"))
  dtt = dtt[,z]
  return(dtt)
}

filtered_data = llply(cancer_data, get_canc_data_for_plot)

get_clin_lnc_cors = function(dtt){
  canc = dtt$Cancer[1]
  print(canc)
  print(dim(dtt))
  #get lncs
  z = which(str_detect(colnames(dtt), "ENSG")) 
  lncs = colnames(dtt)[z]
  
  #look at individual lncRNAs 
  get_cor = function(lnc){
    z = which((str_detect(colnames(dtt), "ENSG") & !(colnames(dtt) %in% lnc)))
    new_dat = dtt
    if(length(z) > 0){
    new_dat = dtt[,-z]}
    #add 0/1 labels 
    new_dat$lncRNA_tag = ""
    med = median(new_dat[,which(colnames(new_dat) %in% lnc)])
    k = which(colnames(new_dat) %in% lnc)
    if(med ==0){
        #if median = 0 then anyone greater than zero is 1 
        l1 = which(new_dat[,k] > 0)
        l2 = which(new_dat[,k] ==0)
        new_dat$lncRNA_tag[l1] = 1
        new_dat$lncRNA_tag[l2] = 0
        }

        if(!(med ==0)){
        l1 = which(new_dat[,k] >= med)
        l2 = which(new_dat[,k] < med)
        new_dat$lncRNA_tag[l1] = 1
         new_dat$lncRNA_tag[l2] = 0
        }
    #get risk type 
    z = as.numeric(which((allCands$cancer %in% canc) & (allCands$gene %in% lnc) & (allCands$data == "TCGA")))
    hr = as.numeric(allCands$HR[z])
    new_dat$risk = ""
    if(hr >1){new_dat$risk = "HighExp"}
    if(hr <1){new_dat$risk = "LowExp"}
    
    #get xcell data 
    z = which(colnames(xcell) %in% new_dat$patient)
    canc_xcell = t(xcell[,c(1, z)])
    colnames(canc_xcell)=canc_xcell[1,]
    canc_xcell=canc_xcell[-1,]

    #now for each cell type get comparison 
    canc_xcell = as.data.frame(canc_xcell)
    canc_xcell$patient = rownames(canc_xcell)
    canc_xcell = as.data.table(canc_xcell)
    canc_xcell = merge(canc_xcell, new_dat, by="patient")

    #check cell type enrichment between low/high groups
    get_enrich = function(cell_type){
      cell_dat  = as.data.frame(canc_xcell)
      cell_dat = cell_dat[,c(cell_type, 77)]
      cell_t=colnames(cell_dat)[1]
      colnames(cell_dat)[1]="cell_type"
      cell_dat$cell_type = as.numeric(cell_dat$cell_type)
      cell_dat$lncRNA_tag = as.numeric(cell_dat$lncRNA_tag)
      wilcox.test(cell_dat$lncRNA_tag, cell_dat$cell_type, alternative = "g") 
      g = ggboxplot(cell_dat, x="lncRNA_tag", y="cell_type") + stat_compare_means() + ggtitle(paste(lnc, cell_t, canc))
      #print(g)
      cell_dat$lncRNA_tag = as.factor(cell_dat$lncRNA_tag)
      pval_w = as.numeric(tidy(wilcox.test(cell_type ~ lncRNA_tag, data = cell_dat))[2])
      pval_t = as.numeric(tidy(t.test(cell_type ~ lncRNA_tag, data = cell_dat))[2])
      res = c(lnc, cell_t, pval_w, pval_t, mean(cell_dat$cell_type[cell_dat$lncRNA_tag==1]),mean(cell_dat$cell_type[cell_dat$lncRNA_tag==0]),
       median(cell_dat$cell_type[cell_dat$lncRNA_tag==1]),median(cell_dat$cell_type[cell_dat$lncRNA_tag==0]), canc)
      names(res)=c("lnc", "cell_type", "pval_w", "pval_t", "mean_high_lnc", "mean_low_lnc", "median_high_lnc", "median_low", "canc")
      return(res)
    }

    lnc_all_immune = as.data.table(ldply(llply(2:65, get_enrich)))
    lnc_all_immune$lnc_name= sapply(lnc_all_immune$lnc, get_name)
    lnc_all_immune$median_high_lnc = as.numeric(lnc_all_immune$median_high_lnc)
    lnc_all_immune$mean_high_lnc = as.numeric(lnc_all_immune$mean_high_lnc)
    lnc_all_immune$median_low = as.numeric(lnc_all_immune$median_low)
    lnc_all_immune$mean_low_lnc = as.numeric(lnc_all_immune$mean_low_lnc)
    lnc_all_immune$pval_w = as.numeric(lnc_all_immune$pval_w)
    lnc_all_immune$pval_t = as.numeric(lnc_all_immune$pval_t)

    hr = as.numeric(allCands$HR[allCands$gene == lnc])
    if(hr > 1){
      hr="high_exp_bad"
    }
    if(hr < 1){
      hr="low_exp_bad"
    }
    lnc_all_immune$hr = hr

    lnc_all_immune = lnc_all_immune[order(abs(median_high_lnc), abs(mean_high_lnc))]
    canc_xcell_lnc = canc_xcell[,c(1:65,77)]
    heatmap = melt(canc_xcell_lnc, id.vars = c("lncRNA_tag", "patient"))
    heatmap$value=as.numeric(heatmap$value)

    canc_xcell_lnc = canc_xcell[,c(1:65)]
    tags=canc_xcell[,c(1,77)]
    canc_xcell_lnc = as.data.frame(canc_xcell_lnc)
    rownames(canc_xcell_lnc) = canc_xcell_lnc$patient
    canc_xcell_lnc$patient = NULL 

    canc_xcell_lnc[] <- lapply(canc_xcell_lnc, function(x) {
    if(is.character(x)) as.numeric(x) else x
    })
    canc_xcell_lnc = as.matrix(canc_xcell_lnc)

    # Example: grouping from the first letter:
    my_group <- as.character(tags$lncRNA_tag)
    
    library(ComplexHeatmap)

    ha = HeatmapAnnotation(lncrna = my_group, 
    col = list(type = c("0" = "limegreen", "1" = "mistyrose2")), 
    show_annotation_name = TRUE,
    annotation_name_offset = unit(3, "mm"))
    
    heat=Heatmap(t(canc_xcell_lnc), top_annotation = ha, column_names_gp = gpar(fontsize = 1), name=paste(canc, get_name(lnc)))
    print(heat)
    return(lnc_all_immune)
    } #end get_cor

    canc_results = as.data.table(ldply(llply(lncs, get_cor)))
    return(canc_results)
}

pdf("all_lncs_cancers_cell_types_boxplots.pdf", width=15, height=12)
all_cancers_cell_types = as.data.table(ldply(llply(filtered_data, get_clin_lnc_cors)))
saveRDS(all_cancers_cell_types, file="all_cancers_cell_types_xcell_results.rds")
dev.off()








