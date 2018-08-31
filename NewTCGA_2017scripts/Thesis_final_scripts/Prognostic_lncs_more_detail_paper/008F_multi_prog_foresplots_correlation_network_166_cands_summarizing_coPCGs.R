
#------------------------------------------------------------------------------
#This script takes all the PCGs that are significantly differentially expressed
#between lncRNA risk groups and calculates their correlation 
#plots heatmap of this data as well as lncRNA and pcg prognostic
#relationship 
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
val_cands = subset(val_cands, as.numeric(pval) < 0.05)

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

#--------This script --------------------------------------------------

#take pairs of genes that are correlated where the PCG is also prognostic 
#and evaluate survival models with both used as predictors compared to just one
#at a time 

#----------------------------------------------------------------------

cors = readRDS("12_cancers_PCG_signatures_aug14.rds")

#keep only the sig correlations that have at least 1 sig prog pcg
cors = as.data.table(cors)
cors = filter(cors, pcg_hr == "<1" | pcg_hr==">1") 
cors = as.data.table(cors)

#get list of unique lnc-pcg combos (200 unique combos across 6 cancers?)
combos = unique(cors$combo)

#test...
comb = combos[[1]]

compare_hrs = function(comb){

  dat = filter(cors, combo == comb)
  #get surv and exp data for these genes 
  z = which(colnames(all) %in% c(dat$lnc, dat$pcg, "OS", "OS.time", "type", "patient", "Cancer"))
  exp_dat = all[,z]
  exp_dat = subset(exp_dat, Cancer == dat$cancer)

  #label patients by binary variable 
  #pcg 
  pcg = dat$pcg
  z = which(colnames(exp_dat) ==pcg)
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
  lnc = dat$lnc
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

  #km plot just pcg
  gene = pcg
  cancer = exp_dat$type[1]
  #get_km_plot(gene, cancer)

  #km plot both   
  exp_dat$OS = as.numeric(exp_dat$OS)
  exp_dat$OS.time = as.numeric(exp_dat$OS.time)

  #cox model using both
  cox_lnc = coxph(Surv(OS.time, OS) ~ lnc_median, data = exp_dat)
  cox_pcg = coxph(Surv(OS.time, OS) ~ pcg_median, data = exp_dat)
  cox_both = coxph(Surv(OS.time, OS) ~ lnc_median + pcg_median, data = exp_dat)
  #does pcg improve lnc?
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
  res$pcg = pcg
  res$canc = exp_dat$type[1]
  res$improve = ""
  res$improve[res$res <= 0.05] = "yes"
  res$improve[res$res > 0.05] = "no"

  exp_dat$OS.time = exp_dat$OS.time/365
  fit <- survfit(Surv(OS.time, OS) ~ lnc_median + pcg_median, data = exp_dat)
          s <- ggsurvplot(
          title = paste(get_name(lnc), get_name(pcg), exp_dat$type[1]),
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
          print(s)

  return(res)        
          
}#end function

pdf("all_lnc_pcg_cors_wprog_values_wider.pdf", width=12)
res_all_cors = llply(combos, compare_hrs, .progress="text")
dev.off()


res_all_cors2 = ldply(res_all_cors)
#filter for the ones where there was an improvmenet 
res_all_cors2 = as.data.table(res_all_cors2)
res_all_cors2$fdr = p.adjust(res_all_cors2$res)
res_all_cors2 = as.data.table(filter(res_all_cors2, res <= 0.05))

saveRDS(res_all_cors2, file="surv_results_lncRNA_plus_mRNAs_together_aug15.rds")

















