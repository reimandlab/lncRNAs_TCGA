setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/real_elastic_net_runs")

library(data.table)
library(dplyr)
library(plyr)
library(ggpubr)
library(tidyr)
library(broom)
library(ggplot2)
library(stringr)

########################################################################
#1. evaluate c-indicies 
########################################################################

#read in all CINDICES files 
results = list.files(pattern="cindices_.rds")
print(length(results))

#break into cancer types 
get_canc = function(file){
  dat = readRDS(file)
  return(dat)
}

all_res = llply(results, get_canc)
all_res = ldply(all_res)
all_res = as.data.table(all_res)

#MAKE PLOT
#X-axis: cancer type
#Y-axis: cindices 

new_res = as.data.table(all_res %>% gather(all_res, cindex, combined:clinical))

#summary boxplots all cancers 
pdf("cindices_real_march2019.pdf", width=20, height=10)
g = ggboxplot(new_res, "canc", "cindex", fill="all_res", color="black", notch = TRUE)
g =  g + stat_compare_means(aes(group = all_res), label = "p.signif") + theme_minimal()
g = ggpar(g, x.text.angle = 90)
print(g + geom_hline(yintercept=0.5, linetype="dashed", color = "red"))
dev.off()

#save and compare to random shuffled data

saveRDS(all_res, file="all_res_REAL_EN_runs_2203.rds")

#for each cancer type get boxplot 
check_perform = function(cancer){
dat = subset(new_res, canc==cancer)
w = wilcox.test(dat$cindex[dat$all_res=="lncRNAs"], dat$cindex[dat$all_res=="clinical"], alternative="greater")
med_lnc = median(dat$cindex[dat$all_res=="lncRNAs"])
med_clin =  median(dat$cindex[dat$all_res=="clinical"])
t = tidy(w)
t$cancer = cancer
t$med_lnc = med_lnc
t$med_clin = med_clin
imp = med_lnc - med_clin
t$imp = imp
return(t)
}

cancers = unique(new_res$canc)

#RUN
wil = llply(cancers, check_perform)
wil = ldply(wil)
wil = as.data.table(wil)
colnames(wil)[2] = "pval"
wil = wil[order(pval, -imp)]

wil_sig = wil

wil_sig$imp = round(wil_sig$imp, digits=2)
wil_sig = as.data.table(wil_sig)
wil_sig = wil_sig[order(med_lnc)]

########################################################################
#2. get genes 
########################################################################

genes = list.files(pattern = "genes_.rds")

#break into cancer types 
get_canc = function(file){
  dat = readRDS(file)
  return(dat)
}

all_res = llply(genes, get_canc)
all_res = ldply(all_res)
all_res = as.data.table(all_res)

rounds = unique(all_res$round)

get_fdr = function(r){
  canc = as.data.table(filter(all_res, round == r))
  #get only those with sig inference post selective p-values 
  canc$c = p.adjust(canc$wald_p, method="fdr")
  canc = as.data.table(filter(canc, c < 0.05))
  return(canc)
}

rounds = as.data.table(ldply(llply(rounds, get_fdr)))






