library(data.table)
library(dplyr)
library(plyr)
library(ggpubr)
library(tidyr)
library(broom)
library(ggplot2)
library(stringr)
library(wesanderson)

#data 
cindices = readRDS("all_res_REAL_EN_1000_runs_0606.rds")
new_res = as.data.table(cindices %>% gather(cindices, cindex, combined:clinical))

#genes 
genes = readRDS("lncRNAs_selected_by_EN_april14.rds")

#get cancer types 
canc_conv = readRDS("canc_conv.rds")
colnames(canc_conv)[2] = "cancer"

#----------------------------------------------------------------
#ANALYSIS 
#----------------------------------------------------------------

#for each cancer type get boxplot 
check_perform = function(cancer){
  dat = subset(new_res, canc==cancer)
  w = wilcox.test(dat$cindex[dat$cindices=="lncRNAs"], dat$cindex[dat$cindices=="clinical"], alternative="greater")
  med_lnc = median(dat$cindex[dat$cindices=="lncRNAs"])
  med_clin =  median(dat$cindex[dat$cindices=="clinical"])
  med_combo = median(dat$cindex[dat$cindices=="combined"])
  t = tidy(w)
  t$cancer = cancer
  t$med_lnc = med_lnc
  t$med_clin = med_clin
  t$med_combo = med_combo
  imp = med_lnc - med_clin
  imp_combo_clin = med_combo - med_clin
  t$imp = imp
  t$imp_combo_clin = imp_combo_clin
  return(t)
}

cancers = unique(new_res$canc)

#RUN
wil = llply(cancers, check_perform)
wil = ldply(wil)
wil = as.data.table(wil)
colnames(wil)[2] = "pval"
wil$fdr = p.adjust(wil$pval, method="fdr")
wil = wil[order(pval, -imp)]

wil_sig = wil

wil_sig$imp = round(wil_sig$imp, digits=2)
wil_sig = as.data.table(wil_sig)
wil_sig = wil_sig[order(med_lnc)]
sig = filter(wil_sig, imp >0, fdr <0.05, med_lnc >=0.5)
z = which(canc_conv$canc %in% sig$cancer)
canc_conv$sig = ""
canc_conv$sig[z] = "V"
canc_conv = canc_conv[z,]
colnames(canc_conv)[2] = "cancer"
#canc_conv = merge(canc_conv, sig, "cancer")
colnames(new_res)[1] = "cancer"
new_res = as.data.table(filter(new_res, cancer %in% canc_conv$cancer))
new_res = merge(new_res, canc_conv, by="cancer")

#pdf("cindices_real_march2019_1000.pdf", width=9, height=3)
g = ggboxplot(new_res, "type", "cindex", bxp.errorbar	=TRUE, 
              fill="cindices", color="black", error.plot = "errorbar", notch = TRUE, width = 0.5, palette=c("grey", "dodgerblue4", "orange"))
g =  g + theme_minimal()
g = ggpar(g, x.text.angle = 65, legend.title="Predictors")
print(g + geom_hline(yintercept=0.5, linetype="dashed", color = "red"))

#dev.off()