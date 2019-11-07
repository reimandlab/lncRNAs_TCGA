setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNAs_2019_manuscript/surv_shuffle_0315")

library(data.table)
library(dplyr)
library(plyr)
library(ggpubr)
library(tidyr)
library(broom)
library(ggplot2)
library(stringr)

#read in real files 
real_dat = readRDS("all_res_REAL_EN_10000_runs_0606.rds")

########################################################################
#1. evaluate c-indicies 
########################################################################

#read in all CINDICES files 
results = list.files(pattern="cindices_.rds")
print(length(results))

#break into cancer types 
get_canc = function(file){
  dat = readRDS(file)
  round = paste(unlist(strsplit(file, "_"))[6:8], collapse="-")
  print(round)
  dat$round = round
  return(dat)
}

all_res = llply(results, get_canc)
all_res = ldply(all_res)
all_res = as.data.table(all_res)

fake_res = all_res
all_res = as.data.table(subset(fake_res, canc %in% real_dat$canc))

#MAKE PLOT
#X-axis: cancer type
#Y-axis: cindices 

new_res = as.data.table(all_res %>% gather(all_res, cindex, combined:clinical))

#summary boxplots all cancers 
pdf("cindices_permutations_march2019.pdf", width=20, height=10)
g = ggboxplot(new_res, "canc", "cindex", fill="all_res", color="black", notch = TRUE)
g =  g + stat_compare_means(aes(group = all_res), label = "p.signif") + theme_minimal()
g = ggpar(g, x.text.angle = 90)
print(g + geom_hline(yintercept=0.5, linetype="dashed", color = "red"))
dev.off()

#save and compare to random shuffled data

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
  canc$c = p.adjust(canc$inference_pvals, method="fdr")
  canc = as.data.table(filter(canc, c < 0.05))
  return(canc)
}

rounds = as.data.table(ldply(llply(rounds, get_fdr)))


########################################################################
#3. COMPARE REAL TO SHUFFLED DATA  
########################################################################

dim(real_dat)
dim(all_res)

cancers = unique(real_dat$canc)

get_comp = function(cancer){
  real = as.data.table(filter(real_dat, canc == cancer))
  fake = as.data.table(filter(fake_res, canc == cancer))
  real$round = "real"
  #need to get median cindex for real dat 
  lncrna_med_real = median(real$lncRNAs)
  lncrna_med_fake = median(fake$lncRNAs)
 
  #get median c-index for each round of permutations 
  # median of groups
  get_med = function(r){
    med = median(as.data.table(filter(fake, round == r))$lncRNAs)
    return(med)
  }
  rounds = unique(fake$round)
  fake_meds = unlist(llply(rounds, get_med))
  t = as.data.table(fake_meds)
  print(paste(cancer, dim(t)[1]))

  #pval = wilcox.test(fake$lncRNAs, real$lncRNAs)$p.value
  pval = wilcox.test(t$fake_meds, lncrna_med_real)$p.value

  p = ggplot(t, aes(x=fake_meds)) + 
  geom_histogram() + xlim(0,1)+
  geom_vline(xintercept=lncrna_med_real, linetype="dashed", color = "red")+
  annotate("text", x = 0.5, y = 5, label = paste("P-val=",round(pval, digits=4)), size=4)+
  ggtitle(cancer)
  print(pval)
  print(cancer)
  #print(p)
  return(p)

}

#
list_plots = llply(cancers, get_comp)

library(gridExtra)
n <- length(list_plots)
nCol <- floor(sqrt(n))

pdf("summary_perms_real_vs_fake.pdf", width=10, height=12)
do.call("grid.arrange", c(list_plots, ncol=nCol))
dev.off()

#dev.off()







