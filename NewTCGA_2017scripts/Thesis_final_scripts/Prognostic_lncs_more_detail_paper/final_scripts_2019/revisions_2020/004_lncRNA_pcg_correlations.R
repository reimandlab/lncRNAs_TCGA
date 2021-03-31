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


#-------------------ANALYSIS--------------------------------------------

#------FEATURES-----------------------------------------------------

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = filter(allCands, data=="TCGA") #179 unique lncRNA-cancer combos, #166 unique lncRNAs
cands_dups = unique(allCands$gene[which(duplicated(allCands$gene))])
allCands$combo=paste(allCands$gene, allCands$cancer, sep="_")
allCands$gene_name = allCands$gene_symbol

cis_pcgs  =readRDS("lncRNA_cands_wPCGs_that_are_in_cis_10kb_nov16.rds")

#convert to ensgs
get_ensg = function(name){
  z = which(allCands$gene_symbol == name)
  return(allCands$gene[z][1])
}
cis_pcgs$lnc = sapply(cis_pcgs$lnc, get_ensg)

get_ensg = function(name){
  z = which(hg38$symbol == name)
  return(hg38$ensgene[z][1])
}
cis_pcgs$pcg = sapply(cis_pcgs$pcg, get_ensg)
cis_pcgs$combo = paste(cis_pcgs$lnc, cis_pcgs$pcg, sep="_")

#check all lncRNA-cis PCG pairs
combos = unique(cis_pcgs$combo) #126 pairs

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

all = as.data.table(all)

check_cis_pcg = function(combo){
  lnc = unlist(str_split(combo, "_"))[1]
  cancer = allCands$cancer[which(allCands$gene == lnc)]
  pcgg = unlist(str_split(combo, "_"))[2]
  res_all = as.data.frame(matrix(ncol=20, nrow=1))

 colnames(res_all) = c(
    "res","res2","lnc","pcg",
    "canc","pcg_improvelnc" ,"lnc_improvepcg" ,"rho",
    "rho_p","perc_lnc_off","perc_pcg_off","lncAIC",
   "lncConcordance" ,"pcgAIC","pcgConcordance" ,"combo_concord",
   "lnc_pval", "pcg_pval",
   "lnc_median","pcg_median")

  cancers=cancer
  print(combo)

  for(i in 1:length(cancers)){
  print(cancers[i])

  #get surv and exp data for these genes
  z = which(colnames(all) %in% c(lnc, pcgg, "OS", "OS.time", "type", "patient", "Cancer"))
  exp_dat = all[,..z]
  exp_dat = subset(exp_dat, Cancer == cancers[i])

  #label patients by binary variable
  z = which(colnames(exp_dat) ==pcgg)
  med = median(unlist(exp_dat[,..z]))
  exp_dat$pcg_median = ""
       if(med ==0){
        #if median = 0 then anyone greater than zero is 1
        l1 = which(exp_dat[,..z] > 0)
        l2 = which(exp_dat[,..z] ==0)
        exp_dat$pcg_median[l1] = 1
        exp_dat$pcg_median[l2] = 0
        }

      if(!(med ==0)){
        l1 = which(exp_dat[,..z] >= med)
        l2 = which(exp_dat[,..z] < med)
        exp_dat$pcg_median[l1] = 1
        exp_dat$pcg_median[l2] = 0
        }
  pcg_median = med
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

  z = which(colnames(exp_dat) == pcgg)
  exp_dat[,z] = log1p(exp_dat[,..z])
  colnames(exp_dat)[z] = "PCG"

  z = which(colnames(exp_dat) == lnc)
  exp_dat[,z] = log1p(exp_dat[,..z])
  colnames(exp_dat)[z] = "LNC"

  #get correlation between them (1)
  rho = rcorr(exp_dat$LNC, exp_dat$PCG, type="spearman")$r[2]
  rho_p = rcorr(exp_dat$LNC, exp_dat$PCG, type="spearman")$P[2]

  g =  ggscatter(exp_dat, x = "LNC", y = "PCG",
          color = "black", shape = 21, size = 3, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "spearman", label.sep = "\n"))+
  ggtitle(paste(allCands$gene_symbol[allCands$gene==lnc][1], get_name_pcg(pcgg), cancers[i]))
  print(g)

  check1 = table(exp_dat$pcg_median)[1] >= 10
  check2 = table(exp_dat$pcg_median)[2] >= 10

  if(dim(table(exp_dat$pcg_median))>1){

  #cox model using both
  cox_lnc = coxph(Surv(OS.time, OS) ~ lnc_median, data = exp_dat)
  cox_pcg = coxph(Surv(OS.time, OS) ~ pcg_median, data = exp_dat)
  cox_both = coxph(Surv(OS.time, OS) ~ lnc_median + pcg_median, data = exp_dat)
  #does pcgg improve lnc?
  res = as.numeric(tidy(anova(cox_lnc, cox_both))[2,4]) #if sig, then yes

  #does lnc improve pcg?
  res2 = as.numeric(tidy(anova(cox_pcg, cox_both))[2,4]) #if sig, then yes

  #arrange
  library(patchwork)

  #train model using both of these genes and compare performance to just using lncRNA
  res = as.data.frame(res)
  res$res2 = res2
  res$lnc = allCands$gene_symbol[allCands$gene==lnc][1]
  res$pcg = get_name_pcg(pcgg)
  res$canc = exp_dat$type[1]

  res$pcg_improvelnc = ""
  res$pcg_improvelnc[res$res <= 0.05] = "yes"
  res$pcg_improvelnc[res$res > 0.05] = "no"

  res$lnc_improvepcg = ""
  res$lnc_improvepcg[res$res2 <= 0.05] = "yes"
  res$lnc_improvepcg[res$res2 > 0.05] = "no"
  res$rho = rho
  res$rho_p = rho_p

   #% of people with zero values for lncRNA
  exp_dat$lnc_exp = ""
  exp_dat$lnc_exp[exp_dat$LNC == 0] = "no_exp"
  exp_dat$lnc_exp[exp_dat$LNC >0] = "exp"
  exp_dat$pcg_exp = ""
  exp_dat$pcg_exp[exp_dat$PCG == 0] = "no_exp"
  exp_dat$pcg_exp[exp_dat$PCG >0] = "exp"

  res$perc_lnc_off = ""
  res$perc_lnc_off = length(which(exp_dat$lnc_exp == "no_exp"))/dim(exp_dat)[1]

  #% of people with zero values for pcg
  res$perc_pcg_off = ""
  res$perc_pcg_off = length(which(exp_dat$pcg_exp == "no_exp"))/dim(exp_dat)[1]

  #get concordance and AIC values
  res$lncAIC = as.data.table(glance(cox_lnc))$AIC
  res$lncConcordance = as.data.table(glance(cox_lnc))$concordance
  res$pcgAIC= as.data.table(glance(cox_pcg))$AIC
  res$pcgConcordance = as.data.table(glance(cox_pcg))$concordance
  res$combo_concord = as.data.table(glance(cox_both))$concordance

  res$lnc_pval = as.data.table(glance(cox_lnc))$p.value.wald
  res$pcg_pval = as.data.table(glance(cox_pcg))$p.value.wald

  res$lnc_median = lnc_median
  res$pcg_median = pcg_median
  res_all = rbind(res_all, res)
  }
}

res_all = res_all[-1,]
return(res_all)

}#end function

pdf("/u/kisaev/Jan2021/all_cis_antisense_lnc_pairs_survival_results_10kb_nov16.pdf", width=6, height=6)
results = llply(combos, check_cis_pcg, .progress="text")
dev.off()

#res = does pcg improve lnc
#res2 = does lnc improve pcg

results2 = ldply(results)
results2$fdr_pcg_improve_lnc_anova = p.adjust(results2$res, method="fdr")
results2$fdr_lnc_improve_pcg_anova = p.adjust(results2$res2, method="fdr")
results2$rho_fdr = p.adjust(results2$rho_p, method="fdr")

results2$lnc_pval_fdr = p.adjust(results2$lnc_pval, method="fdr")
results2$pcg_pval_fdr = p.adjust(results2$pcg_pval, method="fdr")

results2 = as.data.table(results2)
results2 = results2[order(fdr_pcg_improve_lnc_anova)]

length(which(results2$res <= 0.05))
length(which(results2$fdr_pcg_improve_lnc_anova <= 0.05)) #7/129 pairs, the lncRNA beneifts from signal from neighboring gene when accounting for multiple testing correction
length(which(results2$fdr_lnc_improve_pcg_anova <= 0.05)) #even after multiple testing correction, 129/129 pairs improve from adding lncRNA expression to PCG expression

results2$pcg_improvelnc = ""
results2$pcg_improvelnc[results2$fdr_pcg_improve_lnc_anova <= 0.05] = "yes"
results2$pcg_improvelnc[results2$fdr_pcg_improve_lnc_anova > 0.05] = "no"

results2$lnc_improvepcg = ""
results2$lnc_improvepcg[results2$fdr_lnc_improve_pcg_anova <= 0.05] = "yes"
results2$lnc_improvepcg[results2$fdr_lnc_improve_pcg_anova > 0.05] = "no"

#saveRDS(results2, file="110_cis_antisense_pairs_survival_results_aug28.rds")
saveRDS(results2, file="127_cis_antisense_pairs_survival_results_10kb_nov16.rds")
write.csv(results2, file="/u/kisaev/151_cis_antisense_pairs_survival_results_10kb.csv", quote=F, row.names=F)
write.csv(results2, file="/u/kisaev/Jan2021/129_cis_antisense_pairs_survival_results_10kb.csv", quote=F, row.names=F)

###START HERE###-----------------------------------------------------------------

#which of these are candidates?
r = readRDS("127_cis_antisense_pairs_survival_results_10kb_nov16.rds")
cancs = rna[,c("type", "Cancer")]
cancs = unique(cancs)
colnames(cancs) = c("canc", "Cancer")

r = merge(r, cancs, by="canc")
r$combo = paste(r$lnc, r$Cancer, sep="_")

#z = which(r$combo %in% allCands$combo)
#cands_pairs = r[z,]
cands_pairs = r
cands_pairs$pcg_combo = paste(cands_pairs$pcg, cands_pairs$Cancer, sep="_")

#make summary plot
#any of these in cancer gene census list?
#COSMIC cancer gene census
#census = read.csv("Census_allFri_Jul_13_16_55_59_2018.csv")
#get ensg
get_census_ensg = function(genes){
  glist = unlist(strsplit(genes, ","))
  z = which(str_detect(glist, "ENSG"))
  ensg = glist[z]
  return(ensg)
}
census$ensg = sapply(census$Synonyms, get_census_ensg)

z = which(cands_pairs$pcg %in% census$Gene.Symbol)
cands_pairs$census = ""
cands_pairs$census[z] = "YES"

#cands_pairs$pcg = unlist(llply(cands_pairs$pcg, get_name_pcg))
#cands_pairs$lnc = unlist(llply(cands_pairs$lnc, get_name))

prog_pcgs = cands_pairs

prog_pcgs$cor[prog_pcgs$rho <0] = "Negative"
prog_pcgs$cor[prog_pcgs$rho >0] = "Positive"
prog_pcgs$cor[prog_pcgs$rho_fdr > 0.05] = "NS"

prog_pcgs$rho_fdr[prog_pcgs$rho_fdr <0.0000000001] = 0.000001
prog_pcgs$lnc_pcg = paste(prog_pcgs$lnc, prog_pcgs$pcg, prog_pcgs$canc, sep="/")


get_ensg = function(name){
  z = which(allCands$gene_symbol == name)
  return(allCands$gene[z][1])
}

#add what kind of lncRNA it is
prog_pcgs$type_lnc = ""
for(i in 1:nrow(prog_pcgs)){
  lnc = prog_pcgs$combo[i]
  lnc = unlist(strsplit(lnc, "_"))[1]
  lnc_type = fantom$CAT_geneClass[which(fantom$gene == get_ensg(lnc))]
  prog_pcgs$type_lnc[i] = lnc_type
}

prog_pcgs$diff_cindices = prog_pcgs$lncConcordance - prog_pcgs$pcgConcordance
prog_pcgs$sig_spearman = ""
prog_pcgs$sig_spearman[prog_pcgs$rho_fdr < 0.05] = "sig"

prog_pcgs$wald_pval_fdr[prog_pcgs$pcg_pval_fdr < 0.05] = "sig"
prog_pcgs$wald_pval_fdr[prog_pcgs$pcg_pval_fdr > 0.05] = "pcg_not_sig"
colnames(prog_pcgs)[1] = "type"
prog_pcgs$anova_sig_combo = prog_pcgs$lnc_improvepcg
prog_pcgs$anova_sig_combo[prog_pcgs$anova_sig_combo == "yes"] = "Sig"

pdf("/u/kisaev/Jan2021/figure2E_summary_lncs_pcgs_antisense_10kb_nov16.pdf", width=5,height=5)
g = ggplot(prog_pcgs, aes(pcgConcordance, lncConcordance, label=lnc_pcg)) +
  geom_point(aes(colour=type,
     shape=anova_sig_combo), size=1.75, show.legend = FALSE) +
     scale_shape_manual(values = c(17, 5)) +
     colScale+
    xlab("Neighbour PCG Concordance") + ylab("lncRNA Concordance") +
    theme(axis.text = element_text(size=13),
    legend.position = "none") +
     xlim(0.45,0.75) + ylim(0.45,0.75) + geom_abline(intercept=0) +
     geom_vline(xintercept=0.5, linetype="dashed", color = "red") +
     geom_hline(yintercept=0.5, linetype="dashed", color = "red")+
     theme_bw()
     #geom_text_repel(data = subset(prog_pcgs, lncConcordance > 0.7 | pcgConcordance > 0.65), size=2, nudge_y = 0.1,
      #direction = "x",segment.color = "grey50",
      #segment.size = 0.05)
g
dev.off()

prog_pcgs$lnc_better = ""
prog_pcgs$lnc_better[(prog_pcgs$fdr_lnc_improve_pcg_anova < 0.05) & (prog_pcgs$lncConcordance > prog_pcgs$pcgConcordance)] = "lnc_better"
prog_pcgs$lnc_better[is.na(prog_pcgs$lnc_better)] = ""
prog_pcgs$improv_perc = round((prog_pcgs$lncConcordance - prog_pcgs$pcgConcordance)/ prog_pcgs$pcgConcordance*100, digits=2)

pdf("/u/kisaev/Jan2021/figu2_lncs_pcgs_histogram_spearman_rho.pdf", width=5, height=4)
gghistogram(prog_pcgs, x = "rho", position="stack",
   color = "sig_spearman", fill = "sig_spearman",
   palette = c("grey", "black"))+theme_bw()+xlim(-1,1)+
geom_vline(xintercept=0, linetype="dashed", color = "black")
dev.off()

pdf("/u/kisaev/Jan2021/figu2_lncs_pcgs_histogram_cindex_improvement.pdf", width=6, height=5)
gghistogram(prog_pcgs, x = "diff_cindices",
   color = "lnc_better", fill = "lnc_better",
   palette = c("grey", "black"))+theme_bw()+xlim(-0.5,0.5)+
geom_vline(xintercept=0, linetype="dashed", color = "black")
dev.off()

pdf("/u/kisaev/Jan2021/figure2E_summary_lncs_pcgs_antisense_10kb_wSpearmanRho.pdf", width=5,height=5)
g = ggplot(prog_pcgs, aes(diff_cindices, rho, label=lnc_pcg)) +
 geom_point(aes(shape= sig_spearman, color=lnc_better), alpha=0.6, size=2)+
 theme_bw()+
    scale_colour_manual(values = c("black", "blue", "red", "purple")) +
    ylab("Spearman rho") + xlab("lncRNA cindex - PCG cindex") +
    theme(legend.box = "horizontal", axis.text = element_text(size=13),
      legend.text=element_text(size=10), legend.title=element_text(size=10))+
  geom_vline(xintercept=0, linetype="dashed", color = "red")+
      geom_hline(yintercept=0, linetype="dashed", color = "red")
ggpar(g, legend="bottom")
dev.off()

#summarize anova
anvova_plot = prog_pcgs[,c("lnc", "pcg", "fdr_lnc_improve_pcg_anova", "fdr_pcg_improve_lnc_anova", "canc")]
anvova_plot = melt(anvova_plot)
anvova_plot$pair = paste(anvova_plot$lnc, anvova_plot$pcg, anvova_plot$canc)
anvova_plot$log10 = -log10(anvova_plot$value)
anvova_plot = anvova_plot[order(log10)]
anvova_plot$pair = factor(anvova_plot$pair, levels=unique(anvova_plot$pair))
#anvova_plot$variable = factor(anvova_plot$variable, levels=unique(anvova_plot$variable))

means=as.data.table(anvova_plot %>% group_by(canc) %>% dplyr::summarize(mean = mean(log10)))
means = means[order(-mean)]
anvova_plot$canc = factor(anvova_plot$canc, levels=unique(means$canc))
anvova_plot = anvova_plot[order(-log10)]
anvova_plot$pair = factor(anvova_plot$pair, levels=unique(anvova_plot$pair))

sums=as.data.table(anvova_plot %>% group_by(canc, pair) %>% dplyr::summarize(tot = sum(log10)))
sums=sums[order(-tot)]
anvova_plot$pair = factor(anvova_plot$pair, levels=unique(sums$pair))

pdf("/u/kisaev/Jan2021/lncs_pcgs_antisense_10kb_anova_tests.pdf", width=10,height=5)

g = ggbarplot(anvova_plot, "pair", "log10",
  fill = "variable", color = "variable",
  palette = c("#00AFBB", "#E7B800"))
g = ggpar(g,
 font.tickslab = c(4,"plain", "black"),
 xtickslab.rt = 90)+ylab("-log10(anova FDR)")
g + facet_grid(~canc, scales = "free", space = "free") + theme(text = element_text(size=4))+
 geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red")

dev.off()

#median on off balance plot how many are binary lncs/pcgs
meds_plot = prog_pcgs[,c("lnc", "pcg", "canc", "perc_lnc_off", "perc_pcg_off")]
meds_plot = melt(meds_plot)
meds_plot$pair = paste(meds_plot$lnc, meds_plot$pcg, meds_plot$canc)
meds_plot$canc = factor(meds_plot$canc, levels=unique(means$canc))
meds_plot$pair = factor(meds_plot$pair, levels=unique(anvova_plot$pair))

pdf("/u/kisaev/Jan2021/lncs_pcgs_antisense_10kb_linearity_fraction_lncs.pdf", width=10,height=3)

g = ggbarplot(filter(meds_plot, variable == "perc_lnc_off"), "pair", "value",
  fill = "variable", color = "variable",
  palette = c("#00AFBB"))
g = ggpar(g,
 font.tickslab = c(4,"plain", "black"),
 xtickslab.rt = 90)+ylab("% of patients with no detected expression")
g + facet_grid(~canc, scales = "free", space = "free") + theme(text = element_text(size=4))
dev.off()

pdf("/u/kisaev/Jan2021/lncs_pcgs_antisense_10kb_linearity_fraction_pcgs.pdf", width=10,height=3)

g = ggbarplot(filter(meds_plot, variable == "perc_pcg_off"), "pair", "value",
  fill = "variable", color = "variable",
  palette = c("#E7B800"))
g = ggpar(g,
 font.tickslab = c(4,"plain", "black"),
 xtickslab.rt = 90)+ylab("% of patients with no detected expression")
g + facet_grid(~canc, scales = "free", space = "free") + theme(text = element_text(size=4))
dev.off()

dot_plot=prog_pcgs[,c("lnc", "pcg", "canc", "perc_lnc_off", "perc_pcg_off")]
dot_plot$max_perc = ""

for(i in 1:nrow(dot_plot)){
  print(i)
  max_perc = max(dot_plot[i,4:5])
  dot_plot$max_perc[i] = max_perc
}
dot_plot$max_perc = as.numeric(dot_plot$max_perc)

#z= which(dot_plot$lnc == "HOXA-AS4")
#dot_plot$lnc[z] = "HOXA10-AS"

mypal = c("#E5DFD9","#EAD286" ,"#D1EB7B", "#96897F" ,"#E5C0A6" ,
  "#72A93B", "#74DAE3" ,"#49B98D" ,"#D97B8F" ,"#70A2A4", "#64709B" ,"#DFBF38" ,"#61EA4F" ,
  "#C7CBE7", "#786DDA",
"#CFA0E0" ,"#67E9D0" ,"#7C9BE1", "#D94753" ,
"#AAE6B0", "#D13BDF" ,"#DEAEC7" ,"#BBE6DF" ,"#B2B47A" ,"#E6ECBA", "#C86ED7",
 "#7BEE95" ,"#6F46E6" ,"#65B9E0", "#C0EC3E",
"#DE8D54" ,"#DF4FA6")

dot_plot = dot_plot[order(-max_perc)]
dot_plot$lnc = factor(dot_plot$lnc, levels=unique(dot_plot$lnc))
dot_plot$pcg = factor(dot_plot$pcg, levels=unique(dot_plot$pcg))
dot_plot$canc = factor(dot_plot$canc, levels=unique(dot_plot$canc))

#visualize for each pair

 g = ggplot(dot_plot, aes(y = pcg,
             x = lnc)) +        ## global aes
  geom_tile(aes(fill = canc)) +         ## to get the rect filled
  geom_point(aes(colour = max_perc), size=1)  +    ## geom_point for circle illusion
  scale_color_gradient(low = "white",
                       high = "black")+       ## color of the corresponding aes
  theme_bw()+scale_fill_manual(values=mypal)

g = ggpar(g, font.tickslab = c(5,"plain", "black"),
 xtickslab.rt = 90)

#g = ggplot(dot_plot, aes(y = pcg,
 #            x = lnc)) +        ## global aes
  #geom_tile(aes(fill = max_perc, colour = canc,, width=0.8, height=0.8), size=1)+ theme_bw()  +
 # scale_fill_gradient(low = "white", high = "black", na.value = 'white') +
 # coord_equal()
 # theme(legend.position="bottom") + theme(text = element_text(size=7))

#g = ggpar(g, font.tickslab = c(4,"plain", "black"),
# xtickslab.rt = 90)

#print(g)
#dev.off()

pdf("/u/kisaev/Jan2021/lncRNA_PCG_correlations_nonlinearity_dotplot.pdf", width=12, height=9)

g = ggplot(data=dot_plot) +
  geom_tile(aes(x=lnc, y=pcg, fill=max_perc), color="purple") +
  theme_bw() +
  scale_color_gradient(low = "white",
                       high = "black")+
  #facet_grid(~canc,
  #           scales="free", space="free") +
  facet_wrap(vars(canc), scales="free") +
  theme(axis.text.x = element_text(angle = 90,vjust=0.5,hjust = 1)) +
  # switch off clipping
  coord_cartesian(clip = "off")+
  scale_fill_gradient(low = "white",
                       high = "black")
g = ggpar(g, font.tickslab = c(5,"plain", "black"),
 xtickslab.rt = 90) + theme(text = element_text(size=6))
print(g)
dev.off()

plot = as.data.table(filter(prog_pcgs, Cancer == "Brain Lower Grade Glioma", rho_fdr < 0.05))
plot$sig_rho = ""
plot$sig_rho[plot$rho_fdr < 0.05] = "*"
plot$diff = plot$lncConcordance-plot$pcgConcordance
z = which((plot$diff > 0) & (plot$fdr_lnc_improve_pcg_anova< 0.05))
plot$lncbetter = ""
plot$lncbetter[z] = "yes"
plot$lncbetter[-z] = "no"

t=as.data.table(table(plot$lnc))
t=t[order(-N)]
plot$lnc = factor(plot$lnc, levels=t$V1)
plot$lncbetter = factor(plot$lncbetter, levels=c("yes", "no"))
plot$rho = round(plot$rho, digits=2)

pdf("/u/kisaev/Jan2021/figure2X_LGG_lncs_pcgs_correlations.pdf", width=6, height=6)
gene_exp = ggplot(plot, aes(lnc, pcg)) +
  geom_tile(aes(fill = rho, width=0.7, height=0.7), size=1) + theme_bw() #colour =lncbetter
gene_exp = ggpar(gene_exp, x.text.angle = 90) +
#scale_fill_gradientn(values=c(1, .6, .5, .4, 0), colours=c("red", "indianred", "white", "lightblue", "royalblue"))+
#scale_fill_gradientn(colours=c("royalblue", "lightblue", "white", "indianred", "red"),
#                       values = scales::rescale(c(-0.5, -0.05, 0, 0.05, 0.5)), na.value = 'black') +
scale_fill_gradient2(midpoint=0, low = "blue", mid = "white", high = "red")+
ylab("Protein-coding gene")+
xlab("lncRNA") + coord_equal()+ scale_color_manual(values=c("black", "white"))+
theme(legend.position="bottom") + geom_text(aes(label = rho), size=2)+
theme(text = element_text(size=10))
gene_exp
dev.off()

#using all cancers

plot = as.data.table(filter(prog_pcgs, rho_fdr < 0.05))
plot$sig_rho = ""
plot$sig_rho[plot$rho_fdr < 0.05] = "*"
plot$diff = plot$lncConcordance-plot$pcgConcordance
z = which((plot$diff > 0) & (plot$fdr_lnc_improve_pcg_anova< 0.05))
plot$lncbetter = ""
plot$lncbetter[z] = "yes"
plot$lncbetter[-z] = "no"
plot$pcg_combo = paste(plot$pcg, plot$canc)
plot$lnc_combo = paste(plot$lnc, plot$canc)

t=as.data.table(table(plot$lnc_combo))
t=t[order(-N)]
plot$lnc_combo = factor(plot$lnc_combo, levels=t$V1)
plot$rho = round(plot$rho, digits=2)

pdf("/u/kisaev/Jan2021/figure2X_all_cancers_lncs_pcgs_correlations.pdf", width=6, height=6)
gene_exp = ggplot(plot, aes(lnc_combo, pcg_combo)) +
  geom_tile(aes(fill = rho, width=0.7, height=0.7), size=1) + theme_bw() #colour =lncbetter
gene_exp = ggpar(gene_exp, x.text.angle = 90) +
scale_fill_gradientn(values=c(1, .6, .5, .4, 0), colours=c("red", "indianred", "white", "lightblue", "royalblue"))+
#scale_fill_gradientn(colours = c("cyan", "black", "red"),
#                       values = scales::rescale(c(-0.5, -0.05, 0, 0.05, 0.5)), na.value = 'white') +
ylab("Protein-coding gene")+
xlab("lncRNA") + coord_equal()+ scale_color_manual(values=c("black", "white"))+
theme(legend.position="bottom") + #geom_text(aes(label = rho), size=2)+
theme(text = element_text(size=3))
gene_exp
dev.off()


#saveRDS(prog_pcgs, file="final_set_126_lncRNAPCG_pairs_nov16.rds")
#write.csv(prog_pcgs, file="126_cis_antisense_pairs_survival_results_10kb_nov16.csv", quote=F, row.names=F)

#make sure it's only updated list of candidates
r = readRDS("127_cis_antisense_pairs_survival_results_10kb_nov16.rds")

#BOXPLOT
#SUMMARY PCG C-INDICES VS LNCRNA C-INDICDES

boxplot=prog_pcgs[,c("pcgConcordance", "lncConcordance", "type_lnc")]
pcg = boxplot
pcg$concordance = pcg$pcgConcordance
pcg$concordance_type = "PCG"
pcg = pcg[,c("concordance", "concordance_type", "type_lnc")]

lnc = boxplot
lnc$concordance = lnc$lncConcordance
lnc$concordance_type = "lncRNA"
lnc = lnc[,c("concordance", "concordance_type", "type_lnc")]

boxplot = rbind(lnc, pcg)

pdf("/u/kisaev/Jan2021/SUPP_FIGURE_BOXPLOT_lncRNAs_vs_NEARBY_PCGs.pdf")
g = ggplot(boxplot, aes(concordance_type, concordance)) +
geom_boxplot(aes(fill=concordance_type), outlier.colour="black", outlier.shape=16,
             outlier.size=2, notch=FALSE)+
    scale_colour_manual(values = c("blue", "dimgrey", "red", "purple")) +
    xlab("Predictor of Survival") + ylab("Concordance") +
    theme(legend.box = "horizontal", axis.text = element_text(size=13),
      legend.text=element_text(size=10), legend.title=element_text(size=10)) + ylim(0,1)+
    geom_hline(yintercept=0.5, linetype="dashed", color = "red")+
    stat_compare_means()
g
dev.off()
