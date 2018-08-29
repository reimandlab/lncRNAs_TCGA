
#------------------------------------------------------------------------------
#Check if the cis pcg of antisense lncRNA is also prognostic 
#in the same way or not 
#------------------------------------------------------------------------------

#source code
source("check_lnc_exp_cancers.R")
library(corrplot)

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

#HiC data 
load("hic_data.rsav")
#remove rownames
rownames(hic_data) = c(1:nrow(hic_data))

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
all <- merge(rna, pcg, by = c("patient", "Cancer", "OS", "OS.time", "type"))
#all = all[,1:25170]

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
#this is all lncRNA and pcg expression data for each lncRNA-cancer combo
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
prog_pcgs = as.data.table(prog_pcgs)
cis_pcgs = readRDS("lncRNA_cands_wPCGs_that_are_in_cis_aug8.rds")
#convert to ensgs
get_ensg = function(name){
  z = which(ucsc$hg19.ensemblToGeneName.value == name)
  return(ucsc$hg19.ensGene.name2[z])
}
cis_pcgs$lnc = sapply(cis_pcgs$lnc, get_ensg)
cis_pcgs$pcg = sapply(cis_pcgs$pcg, get_ensg)
cis_pcgs$combo = paste(cis_pcgs$lnc, cis_pcgs$pcg, sep="_")

#check all lncRNA-cis PCG pairs 
combos = unique(cis_pcgs$combo)
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

check_cis_pcg = function(combo){
  lnc = unlist(str_split(combo, "_"))[1]
  cancer = allCands$cancer[which(allCands$gene == lnc)]
  pcgg = unlist(str_split(combo, "_"))[2]
  pcg_surv = filter(prog_pcgs, gene == pcgg, canc == cancer)
  print(pcg_surv)

  #pcg pvalue 
  pcg_pvalue = as.numeric(pcg_surv$fdr)
  if(length(pcg_pvalue)==0){
    pcg_pvalue = 1
  }

  #get surv and exp data for these genes 
  z = which(colnames(all) %in% c(lnc, pcgg, "OS", "OS.time", "type", "patient", "Cancer"))
  exp_dat = all[,z]
  exp_dat = subset(exp_dat, Cancer == cancer)

  #label patients by binary variable 
  z = which(colnames(exp_dat) ==pcgg)
  med = median(exp_dat[,z])
  exp_dat$pcg_median = ""
       if(med ==0){
        #if median = 0 then anyone greater than zero is 1 
        l1 = which(exp_dat[,z] > 0)
        l2 = which(exp_dat[,z] ==0)
        exp_dat$pcg_median[l1] = 1
        exp_dat$pcg_median[l2] = 0
        }

      if(!(med ==0)){
        l1 = which(exp_dat[,z] >= med)
        l2 = which(exp_dat[,z] < med)
        exp_dat$pcg_median[l1] = 1
        exp_dat$pcg_median[l2] = 0
        }
  
  #lncRNA 
  z = which(colnames(exp_dat) ==lnc)
  med = median(exp_dat[,z])
  exp_dat$lnc_median = ""
       if(med ==0){
        #if median = 0 then anyone greater than zero is 1 
        l1 = which(exp_dat[,z] > 0)
        l2 = which(exp_dat[,z] ==0)
        exp_dat$lnc_median[l1] = 1
        exp_dat$lnc_median[l2] = 0
        }

      if(!(med ==0)){
        l1 = which(exp_dat[,z] >= med)
        l2 = which(exp_dat[,z] < med)
        exp_dat$lnc_median[l1] = 1
        exp_dat$lnc_median[l2] = 0
      }

  #generate KM plot for one at a time then both combined
  #km plot just lncRNA
  gene = lnc
  cancer = exp_dat$type[1]
  #get_km_plot(gene, cancer)

  #km plot just pcgg
  gene = pcgg
  cancer = exp_dat$type[1]
  #get_km_plot(gene, cancer)

  #km plot both   
  exp_dat$OS = as.numeric(exp_dat$OS)
  exp_dat$OS.time = as.numeric(exp_dat$OS.time)

  #cox model using both
  cox_lnc = coxph(Surv(OS.time, OS) ~ lnc_median, data = exp_dat)
  cox_pcg = coxph(Surv(OS.time, OS) ~ pcg_median, data = exp_dat)
  cox_both = coxph(Surv(OS.time, OS) ~ lnc_median + pcg_median, data = exp_dat)
  #does pcgg improve lnc?
  res = as.numeric(tidy(anova(cox_lnc, cox_both))[2,4]) #if sig, then yes 

  #compare coefficeints
  lnc_plot = print(ggforest(cox_lnc, data=exp_dat))
  pcg_plot = print(ggforest(cox_pcg, data=exp_dat))
  both_plot = print(ggforest(cox_both, data=exp_dat))

  #arrange
  library(patchwork)
  print(lnc_plot + pcg_plot + both_plot)

  #train model using both of these genes and compare performance to just using lncRNA 
  res = as.data.frame(res)
  res$lnc = lnc
  res$pcg = pcgg
  res$canc = exp_dat$type[1]
  res$improve = ""
  res$improve[res$res <= 0.05] = "yes"
  res$improve[res$res > 0.05] = "no"
  res$pcg_fdr_pval = pcg_pvalue
  print(res)
  exp_dat$OS.time = exp_dat$OS.time/365
  fit <- survfit(Surv(OS.time, OS) ~ lnc_median + pcg_median, data = exp_dat)
          s <- ggsurvplot(
          title = paste(get_name(lnc), get_name(pcgg), exp_dat$type[1]),
          fit, 
          xlab = "Time (Years)", 
          surv.median.line = "hv",
          font.main = c(16, "bold", "black"),
          font.x = c(14, "plain", "black"),
          font.y = c(14, "plain", "black"),
          font.tickslab = c(14, "plain", "black"),
          font.legend = 12,
          risk.table.fontsize = 5, 
          #legend.labs = c("Low Expression", "High Expression"),             # survfit object with calculated statistics.
          data = exp_dat,      # data used to fit survival curves. 
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
          ggtheme = theme_bw(), # customize plot and risk table with a theme.
          risk.table.y.text.col = T, # colour risk table text annotations.
          risk.table.y.text = FALSE # show bars instead of names in text annotations
                            # in legend of risk table
          )
          #print(s)

  return(res)        
          
}#end function


#pdf("all_cis_antisense_lnc_pairs_survival_results_aug20.pdf", width=10)
results = llply(combos, check_cis_pcg, .progress="text")
#dev.off()

results2 = ldply(results)
results2$fdr = p.adjust(results2$res, method="fdr")
results2 = as.data.table(results2)
results2 = results2[order(fdr)]

length(which(results2$res <= 0.05))
length(which(results2$fdr <= 0.05))
length(which(results2$pcg_fdr_pval <= 0.05))

saveRDS(results2, file="110_cis_antisense_pairs_survival_results_aug28.rds")



