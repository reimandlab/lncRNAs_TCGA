
#------------------------------------------------------------------------------
#This script takes all the PCGs that are significantly differentially expressed
#between lncRNA risk groups and calculates their correlation 
#plots heatmap of this data as well as lncRNA and pcg prognostic
#relationship 
#------------------------------------------------------------------------------

#source code
source("check_lnc_exp_cancers.R")

#COSMIC cancer gene census
census = read.csv("Census_allFri_Jul_13_16_55_59_2018.csv")

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

#------FUNCTIONS-----------------------------------------------------

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
    )
}

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
coexp$combo1 = paste(coexp$lnc, coexp$pcg, sep="_")
coexp$combo2 = paste(coexp$pcg, coexp$lnc, sep="_")

#PCG lncRNA results
pcg_lnc = readRDS("summary_pcg_analysis_wHRs_jul2y24.rds") #all these have at least 1, 50-pcg signature 
pcg_lnc = pcg_lnc[order(-NumPCGs)]
pcg_lnc$HR = as.numeric(pcg_lnc$HR)
pcg_lnc$lnc_stat = ""
pcg_lnc$lnc_stat[which(pcg_lnc$HR < 0)] = "Favourable"
pcg_lnc$lnc_stat[which(pcg_lnc$HR > 0)] = "Unfavourable"

#-------------------ANALYSIS--------------------------------------------
#Generate heatmap using those PCGs sig up/down regulated in lncRNA 
#risk or non-risk groups

#For each cancer type get all required data 
#PCG and lncRNA expression
combos = unique(pcg_lnc$combo)
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

#PART2 start 

get_lnc_canc = function(dat){
  cancer = dat$Cancer[1]
  combo = dat$combo[1]
  lnc = unlist(strsplit(combo, "_"))[1]

  pcgs = colnames(pcg)[2:19351]
  #keep only pcgs that are selected to be in lncRNA signature 
  z = which(coexp$combo == combo)
  lnc_pcgs = unique(coexp$pcg[z])

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

##3-----------------get correlation pairs-----------------------------------

#PART3 start 
cancer = cancs[4]
library(tidyverse)

prog_pcgs = readRDS("mRNAs_Survival_Results_prognostic_pcgs_July19.rds")
cis_pcgs = readRDS("lncRNA_cands_wPCGs_that_are_in_cis_aug8.rds")
#convert to ensgs
get_ensg = function(name){
  z = which(ucsc$hg19.ensemblToGeneName.value == name)
  return(ucsc$hg19.ensGene.name2[z])
}
cis_pcgs$lnc = sapply(cis_pcgs$lnc, get_ensg)
cis_pcgs$pcg = sapply(cis_pcgs$pcg, get_ensg)
cis_pcgs$combo = paste(cis_pcgs$lnc, cis_pcgs$pcg, sep="_")

get_summary = function(cancer){
  
  #collect data from all lncRNAs in cancer type 
  keep = c()
  for(i in 1:length(all_canc_lnc_data)){
    z = all_canc_lnc_data[[i]]$canc[1] == cancer
    if(z){
      keep = c(keep, i)
    }
  }

  #cancer data
  canc_dats = all_canc_lnc_data[keep]
  canc_dats = reshape::merge_all(canc_dats)

    print(cancer)
    pcgs = colnames(canc_dats)[which(str_detect(colnames(canc_dats), "ENSG"))]
    lncs = unique(canc_dats$lnc)
    genes = c(pcgs, lncs)
    canc_exp = subset(all, Cancer == cancer)
    rownames(canc_exp) = canc_exp$patient
    #gene expression data for all pcgs enriched in at least 1 lncRNA of that
    #cancer type and lncRNAs themselves
    canc_exp = canc_exp[,which(colnames(canc_exp) %in% c(genes))]
    
    #get their correlations
    res2 = rcorr(as.matrix(canc_exp), type="spearman")
    res2 = flattenCorrMatrix(res2$r, res2$P)
    res2$fdr = p.adjust(res2$p, method="fdr")
    res2 = as.data.table(res2)
    res2 = res2[order(fdr)]

    z = which(res2$row %in% lncs)
    res2 = as.data.frame(res2)
    res2$rowgene = ""
    res2$rowgene[z] = "lncRNA"
    res2$rowgene[-z] = "mRNA"
    z = which(res2$column %in% lncs)
    res2$columngene = ""
    res2$columngene[z] = "lncRNA"
    res2$columngene[-z] = "mRNA"
    if(length(z)==0){
      res2$columngene = "mRNA"
    }

    res2 = as.data.table(res2)
    res2 = as.data.table(dplyr::filter(res2, fdr <= 0.05, abs(cor) >= 0.4))

    #all unique genes that appear
    allgenes = unique(c(res2$row, res2$column))
    #this will be the edge table 
    res2$canc = cancer
    res2$combo1 = paste(res2$row, res2$column, sep="_")
    z1 = which(res2$combo1 %in% coexp$combo1)
    z2 = which(res2$combo1 %in% coexp$combo2)
    z = unique(c(z1,z2))
   
    test = merge(res2, coexp, by=c("combo1", "canc"))
    if(!(dim(test)[1] ==0)){

    #for now keep all correlations wtih at least 1 lncRNA
    #lnc_cors = dplyr::filter(res2, rowgene == "lncRNA" | columngene == "lncRNA")
    #lnc_cors = dplyr::filter(res2, row == "ENSG00000230432" | column == "ENSG00000230432")
    #lnc_cors$canc = cancer
    lnc_cors = as.data.table(test)
    lnc_cors = lnc_cors[order(-(abs(cor)))]
    lnc_cors$combo2 = NULL
    lnc_cors$combo1 = NULL
    lnc_cors$combo = NULL
    lnc_cors$combo = paste(lnc_cors$lnc, lnc_cors$pcg, sep="_")
    z = which(lnc_cors$combo %in% cis_pcgs$combo)
    if(!(length(z)==0)){
    lnc_cors$distance[z] = "cis"
    lnc_cors$distance[-z] = "trans"
    }
    if(length(z)==0){
    lnc_cors$distance = "trans"
    }

    #mRNA - 
    z = which(str_detect(colnames(pcg), "ENSG"))
    pcgs = colnames(pcg)[z]

    nodes = data.frame(gene = unique(c(lnc_cors$row, lnc_cors$column)))
    nodes$type = ""
    z = which(nodes$gene %in% lncs)
    #these nodes also in coexp results
    nodes$type[z] = "lncRNA"
    nodes$type[-z] = "mRNA"

    z = which(canc_dats$lnc %in% unique(nodes$gene[nodes$type == "lncRNA"]))
    lnc_dat = canc_dats[z,(ncol(canc_dats)-4):ncol(canc_dats)]

    hrs = allCands[which(allCands$gene %in% unique(nodes$gene[nodes$type == "lncRNA"]))]
    hrs = hrs[,c(1,3, 12)]
    colnames(hrs)[3] = "fdr"

    #get pcg HRs
    progs = prog_pcgs[which(prog_pcgs$canc == cancer),]
    lnc_pcsg = unique(nodes$gene[nodes$type == "mRNA"])
    z = which(progs$gene %in% lnc_pcsg)
    progs = progs[z,]
    hrs_pcgs = progs[,c(1,3,9)]

    all_hrs = rbind(hrs, hrs_pcgs)
    nodes = merge(nodes, all_hrs, by = "gene")
    nodes$fdr = as.numeric(nodes$fdr)
    nodes$HR = as.numeric(nodes$HR)

    nodes$risk_type[nodes$HR > 1] = "Unfavourable"
    nodes$risk_type[nodes$HR < 1] = "Favourable"
    nodes$cancer = cancer
    colnames(nodes)[1] = "id"
    #write node file
    write.table(nodes, file = "node_file_lncRNA_newtwork.txt", sep=";")

    nodes = fread("node_file_lncRNA_newtwork.txt", sep=";")
    nodes$V1 = NULL
    colnames(nodes) = c("id", "type", "HR", "fdr", "risk", "canc")

    #edge file
    lnc_cors$cancer = lnc_cors$canc
    lnc_cors$canc = NULL
    colnames(lnc_cors)[1:2] = c("from", "to")
    write.table(lnc_cors, file = "edges_file_lncRNA_newtwork.txt", sep=";")

    lnc_cors = fread("edges_file_lncRNA_newtwork.txt", sep=";")
    lnc_cors$V1 = NULL
    colnames(lnc_cors) = c("from", "to", "cor", "pval", "fdr_cor", "from_type", "to_type", "lnc", 
      "pcg", "mean_diff", "pvalue", "fdr_lm", "pcg_risk_group", "combo", "distance", "cancer")
    
    lnc_cors$from = sapply(lnc_cors$from, get_name)
    lnc_cors$to = unlist(sapply(lnc_cors$to, get_name))
    nodes$ensg = nodes$id
    nodes$id = sapply(nodes$id, get_name)
    z = which(is.na(nodes$id))
    if(!(length(z)==0)){
      nodes = nodes[-z,]
    }
    z = which(is.na(lnc_cors$to))
    if(!(length(z)==0)){
      lnc_cors = lnc_cors[-z,]
    }

    #make geom_tile
    nodes$HR[nodes$fdr > 0.05] = 1
    nodes$HR = log2(nodes$HR)

    #add which are census genes 
    z = which(lnc_cors$to %in% census$Gene.Symbol)
    lnc_cors$to = as.character(lnc_cors$to)
    if(!(length(z)==0)){lnc_cors$to[z] = paste("*", lnc_cors$to[z], "*")}

    #order by highest correlation
    lnc_cors = lnc_cors[order(from, cor, distance)]
   
    if(dim(lnc_cors)[1] > 100){
      lnc_cors = lnc_cors[1:100,]
    }

    lnc_cors$from = factor(lnc_cors$from, levels=unique(lnc_cors$from))
    lnc_cors$to = factor(lnc_cors$to, levels=unique(lnc_cors$to))

    #START plot
    main_map = ggplot(lnc_cors, aes(from, to)) + geom_tile(aes(fill = cor, color=distance), size=0.5) + 
    scale_colour_manual(values=c("gray25", "lightsteelblue3")) + 
      scale_fill_gradient2(low = "blue", mid = "white", high = "orange") +
      labs(x = "",y = "") + theme_minimal() +
      theme(axis.text.x = element_text(angle=30,hjust=1,vjust=1.0, size=8),
        axis.text.y = element_text(size = 7),
        axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), legend.position="top", legend.box = "horizontal") + 
      ggtitle(cancer) 

    #add covaraite for the lncRNAs --> whether they are fav/unfavour
    lncs_cov = filter(nodes, type =="lncRNA")  
    z = which(lncs_cov$id %in% lnc_cors$from)
    lncs_cov = lncs_cov[z,]
    lncs_cov$id = factor(lncs_cov$id, levels = unique(lnc_cors$from))
    lnc_conv_pl = ggplot(lncs_cov, aes(id, 1)) + geom_tile(aes(fill = HR)) + 
      scale_fill_gradient2(low = "forestgreen", mid = "white", high = "brown3") +
      labs(x = "",y = "") +
      theme(legend.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_blank(),axis.ticks.y=element_blank()) +
      guides(fill=FALSE) +theme_void()

    #add covaraite for the mRNAs --> whether they are fav/unfavour
    mrna_conv = filter(nodes, type =="mRNA")
    z = which(mrna_conv$id %in% lnc_cors$to)
    mrna_conv = mrna_conv[z,] 
    mrna_conv$id = factor(mrna_conv$id, levels = unique(lnc_cors$to)) 
    mrna_conv_pl = ggplot(mrna_conv, aes(id, 1)) + geom_tile(aes(fill = HR), size=0.5) + 
      coord_flip() +
      scale_fill_gradient2(low = "forestgreen", mid = "white", high = "brown3") +
      labs(x = "",y = "") +
         theme(legend.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_blank(),axis.ticks.y=element_blank())+theme_void()
    
    p1 = main_map + mrna_conv_pl + plot_layout(ncol = 2, widths = c(10, 3))
    p2 = p1 + lnc_conv_pl + plot_layout(ncol = 2, heights = c(10, 0.5))
    print(p2)

    return(lnc_cors)
}#if !(dim(test)[1]==0)
}

pdf("top100_regulatory_pcgs_heatmaps_aug10.pdf", height=9, width=9)
canc_results = llply(cancs, get_summary, .progress = "text")
dev.off()

#remove null
canc_results = Filter(Negate(is.null), canc_results)









