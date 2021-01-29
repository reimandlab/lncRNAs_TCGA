#setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/real_elastic_net_runs_1000")
setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNAs_2021_manuscript/real_elastic_net_runs_1000")

library(data.table)
library(dplyr)
library(plyr)
library(ggpubr)
library(tidyr)
library(broom)
library(ggplot2)
library(stringr)
library(wesanderson)

########################################################################
#1. evaluate c-indicies
########################################################################

#read in all CINDICES files
results = list.files()
print(length(results))

results = results[which(str_detect(results, "cindices_.rds"))]

#break into cancer types
get_canc = function(file){
  dat = readRDS(file)
  print(summary(dat))
  return(dat)
}

all_res = llply(results, get_canc)
all_res = ldply(all_res)
all_res = as.data.table(all_res)

#MAKE PLOT
#X-axis: cancer type
#Y-axis: cindices

new_res = as.data.table(all_res %>% gather(all_res, cindex, combined:clinical))

#save and compare to random shuffled data

#get cancer types
canc_conv = readRDS("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/canc_conv.rds")
colnames(canc_conv)[2] = "canc"
canc_conv = as.data.table(merge(canc_conv, new_res, by="canc"))
canc_conv$cindex= as.numeric(canc_conv$cindex)
##order by medians
res = as.data.table(canc_conv %>%
	group_by(type, all_res) %>% summarise_each(funs(max, min, mean, median, sd), cindex))

res = as.data.table(filter(res, all_res == "lncRNAs"))
res = res[order(median)]
canc_conv$type = factor(canc_conv$type, levels = res$type)
#summary boxplots all cancers

#mypal = wes_palette("FantasticFox")

pdf("/u/kisaev/Jan2021/cindices_real_march2019_1000.pdf", width=9, height=3)
g = ggboxplot(canc_conv, "type", "cindex", bxp.errorbar	=TRUE,
fill="all_res", color="black", error.plot = "errorbar", notch = TRUE, width = 0.5, palette=c("grey", "dodgerblue4", "orange"))
g =  g + theme_minimal()
g = ggpar(g, x.text.angle = 65, legend.title="Predictors")
print(g + geom_hline(yintercept=0.5, linetype="dashed", color = "red"))
dev.off()

saveRDS(all_res, file="all_res_REAL_EN_1000_runs_0606.rds")
#scp all_res_REAL_EN_1000_runs_0606.rds /u/kisaev/
#scp kisaev@chickenwire.oicr.on.ca:/u/kisaev/all_res_REAL_EN_1000_runs_0606.rds /Users/kisaev/Documents/lncRNAs

#for each cancer type get boxplot
check_perform = function(cancer){
dat = subset(new_res, canc==cancer)
w = wilcox.test(dat$cindex[dat$all_res=="lncRNAs"], dat$cindex[dat$all_res=="clinical"], alternative="greater")
med_lnc = median(dat$cindex[dat$all_res=="lncRNAs"])
med_clin =  median(dat$cindex[dat$all_res=="clinical"])
med_combo = median(dat$cindex[dat$all_res=="combined"])
w_combo_clin = wilcox.test(dat$cindex[dat$all_res=="combined"], dat$cindex[dat$all_res=="clinical"], alternative="greater")

 t = tidy(w)
  t$cancer = cancer
  t$med_lnc = med_lnc
  t$med_clin = med_clin
  t$med_combo = med_combo
  imp = med_lnc - med_clin
  imp_combo_clin = med_combo - med_clin
  t$imp = imp
  t$imp_combo_clin = imp_combo_clin
  t$w_combo_clin = w_combo_clin$p.value

  return(t)
}

cancers = unique(new_res$canc)

#RUN
wil = llply(cancers, check_perform)
wil = ldply(wil)
wil = as.data.table(wil)
colnames(wil)[2] = "pval"
wil = wil[order(pval, -imp)]
wil$fdr = p.adjust(wil$pval, method="fdr")
wil$fdr_w_combo_clin = p.adjust(wil$w_combo_clin, method="fdr")

wil_sig = wil

wil_sig$imp = round(wil_sig$imp, digits=2)
wil_sig = as.data.table(wil_sig)
wil_sig = wil_sig[order(med_lnc)]
sig = filter(wil_sig, imp >0, fdr <0.05, med_lnc >=0.5)
z = which(canc_conv$canc %in% sig$cancer)
canc_conv$sig = ""
canc_conv$sig[z] = "V"
canc_conv = canc_conv[z,]

canc_conv$all_res = factor(canc_conv$all_res, levels = c("clinical", "lncRNAs", "combined"))

pdf("/u/kisaev/Jan2021/cindices_real_march2019_1000_pvalues.pdf", width=10, height=3)
g = ggplot(canc_conv, aes(type, cindex, fill=all_res)) +
  geom_boxplot(outlier.alpha = 0.1, color="black") + theme_classic() #+
  #stat_compare_means(aes(group = all_res), label = "p.signif")
g = ggpar(g, x.text.angle = 45, legend.title="Predictors")
print(g + geom_hline(yintercept=0.5, linetype="dashed", color = "red") + scale_fill_brewer(c("grey", "green", "orange")) +
 xlab("Cancer") + ylab("c-index"))
dev.off()


########################################################################
#2. get genes
########################################################################

results = list.files()
print(length(results))
genes = results[which(str_detect(results, "genes_.rds"))]

#break into cancer types
get_canc = function(file){
  dat = readRDS(file)
  dat = as.data.table(filter(dat, num_rounds_selected >=500))
  #dat = as.data.table(filter(dat, num_rounds_selected >=100))
  return(dat)
}

all_res = llply(genes, get_canc)
all_res = ldply(all_res)
all_res = as.data.table(all_res)

rounds = unique(all_res$round)

#remove candidates that are duplicated
#dups = all_res$gene[which(duplicated(all_res$gene))]
#all_res = as.data.table(filter(all_res, !(gene %in% dups)))
#all_res = subset(all_res, canc == "Glioblastoma multiforme")
all_res = all_res[,c("gene", "type", "cancer", "num_rounds_selected")]
#write.csv(all_res, file="gbm_min10perc_rounds_candidates.csv", quote=F, row.names=F)
saveRDS(all_res, file="/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNAs_2019_manuscript/lncRNAs_selected_by_EN_april14.rds")

#qsub -cwd -b y -N CVs -l h_vmem=75g "module load R/3.4.0;Rscript Figure2c_analysis.R"
