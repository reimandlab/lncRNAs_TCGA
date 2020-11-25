#------------------------------------------------------------------------------
#Check if the cis pcg of antisense lncRNA is also prognostic 
#in the same way or not 
#------------------------------------------------------------------------------

library(corrplot)

set.seed(911)

source("/u/kisaev/lncRNAs_TCGA/NewTCGA_2017scripts/Thesis_final_scripts/Prognostic_lncs_more_detail_paper/final_scripts_2019/revisions_2020/load_data.R")
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

#setWD
setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNAs_2019_manuscript")

#-------------------ANALYSIS--------------------------------------------

#------FEATURES-----------------------------------------------------

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = subset(allCands, data == "TCGA") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 
allCands$combo = unique(paste(allCands$gene, allCands$cancer, sep="_"))

lncs = allCands$gene

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

check_med = function(lnc){
 
  cancer = allCands$cancer[which(allCands$gene == lnc)]
 
  #get surv and exp data for these genes 
  z = which(colnames(all) %in% c(lnc, "OS", "OS.time", "type", "patient", "Cancer"))
  exp_dat = all[,..z]
  exp_dat = subset(exp_dat, Cancer == cancer)

  #label patients by binary variable 
  #lncRNA 
  z = which(colnames(exp_dat) ==lnc)
  med = median(unlist(exp_dat[,..z]))
  exp_dat$lnc_median = ""
       if(med ==0){
        #if median = 0 then anyone greater than zero is 1 
        l1 = which(exp_dat[,..z] > 0)
        l2 = which(exp_dat[,..z] ==0)
        exp_dat$lnc_median[l1] = 1
        exp_dat$lnc_median[l2] = 0
        }

      if(!(med ==0)){
        l1 = which(exp_dat[,..z] >= med)
        l2 = which(exp_dat[,..z] < med)
        exp_dat$lnc_median[l1] = 1
        exp_dat$lnc_median[l2] = 0
      }
  
  lnc_median = med

  gene = lnc
  cancer = exp_dat$type[1]

  z = which(colnames(exp_dat) == lnc)
  exp_dat[,z] = log1p(exp_dat[,..z])
  colnames(exp_dat)[z] = "LNC"

  #cox model using both
  cox_lnc = coxph(Surv(OS.time, OS) ~ lnc_median, data = exp_dat)

  res = summary(cox_lnc)$coefficients[2]

  #train model using both of these genes and compare performance to just using lncRNA 
  res = as.data.frame(res)
  res$lnc = get_name(lnc)
  res$canc = exp_dat$type[1]
  
  print(res)

  res$lnc_median = lnc_median

  #% of people with zero values for lncRNA 
  exp_dat$lnc_exp = ""
  exp_dat$lnc_exp[exp_dat$LNC == 0] = "no_exp"
  exp_dat$lnc_exp[exp_dat$LNC >0] = "exp"

  res$perc_lnc_off = ""
  res$perc_lnc_off = length(which(exp_dat$lnc_exp == "no_exp"))/dim(exp_dat)[1]
  
  return(res)        
}          

results = llply(lncs, check_med, .progress="text")
results2 = ldply(results)
results2 = as.data.table(results2)
results2 = results2[order(-perc_lnc_off)]
results2$lnc = factor(results2$lnc, levels=results2$lnc)
results2$canc = factor(results2$canc, levels=unique(results2$canc))

#visualize for each pair 
pdf("/u/kisaev/all_lncRNAs_zero_expression_summary.pdf", width=8, height=6)

g = ggbarplot(results2, "lnc", "perc_lnc_off",
  fill = "canc", color = "black",
  palette = mypal)
g = ggpar(g, legend="none", 
 font.tickslab = c(4,"plain", "black"),
 xtickslab.rt = 90)+ylab("% of patients with no detected expression")
g + facet_grid(~canc, scales = "free", space = "free") + theme(text = element_text(size=4))

dev.off()

#make histogram 
pdf("/u/kisaev/all_lncRNAs_zero_expression_summary_histogram.pdf", width=8, height=6)
gghistogram(results2, x="perc_lnc_off", color="black", fill="#00AFBB")
dev.off()