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
allCands = subset(allCands, data == "TCGA") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 
allCands$combo = unique(paste(allCands$gene, allCands$cancer, sep="_"))

cis_pcgs  =readRDS("lncRNA_cands_wPCGs_that_are_in_cis_10kb_nov16.rds")
colnames(cis_pcgs)[c(4,5)] = c("lnc_strand", "lnc")

#convert to ensgs
get_ensg = function(name){
  z = which(ucsc$hg19.ensemblToGeneName.value == name)
  return(ucsc$hg19.ensGene.name2[z][1])
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
  #pcg_surv = filter(prog_pcgs, gene == pcgg, canc == cancer)
  #print(pcg_surv)

  #pcg pvalue 
  #pcg_pvalue = as.numeric(pcg_surv$fdr)
  #if(length(pcg_pvalue)==0){
  #  pcg_pvalue = 1
  #}

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

  z = which(colnames(exp_dat) == pcgg)
  exp_dat[,z] = log1p(exp_dat[,z])
  colnames(exp_dat)[z] = "PCG"

  z = which(colnames(exp_dat) == lnc)
  exp_dat[,z] = log1p(exp_dat[,z])
  colnames(exp_dat)[z] = "LNC"

  #get correlation between them (1)
  rho = rcorr(exp_dat$LNC, exp_dat$PCG, type="spearman")$r[2]
  rho_p = rcorr(exp_dat$LNC, exp_dat$PCG, type="spearman")$P[2]

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

  #compare coefficeints
  #lnc_plot = print(ggforest(cox_lnc, data=exp_dat))
  #pcg_plot = print(ggforest(cox_pcg, data=exp_dat))
  #both_plot = print(ggforest(cox_both, data=exp_dat))

  #arrange
  library(patchwork)

  #train model using both of these genes and compare performance to just using lncRNA 
  res = as.data.frame(res)
  res$res2 = res2
  res$lnc = lnc
  res$pcg = pcgg
  res$canc = exp_dat$type[1]
  
  res$pcg_improvelnc = ""
  res$pcg_improvelnc[res$res <= 0.05] = "yes"
  res$pcg_improvelnc[res$res > 0.05] = "no"
  
  res$lnc_improvepcg = ""
  res$lnc_improvepcg[res$res2 <= 0.05] = "yes"
  res$lnc_improvepcg[res$res2 > 0.05] = "no"

  #res$pcg_fdr_pval = pcg_pvalue #Fdr from cox model
  print(res)

  res$rho = rho
  res$rho_p = rho_p
  
  #get concordance and AIC values
  res$lncAIC = as.data.table(glance(cox_lnc))$AIC
  res$lncConcordance = as.data.table(glance(cox_lnc))$concordance
  res$pcgAIC= as.data.table(glance(cox_pcg))$AIC
  res$pcgConcordance = as.data.table(glance(cox_pcg))$concordance

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
}          
}#end function


#pdf("all_cis_antisense_lnc_pairs_survival_results_10kb_nov16.pdf", width=10)
results = llply(combos, check_cis_pcg, .progress="text")
#dev.off()

results2 = ldply(results)
results2$fdr_res = p.adjust(results2$res, method="fdr")
results2$fdr_res2 = p.adjust(results2$res2, method="fdr")
results2$rho_fdr = p.adjust(results2$rho_p, method="fdr")

results2 = as.data.table(results2)
results2 = results2[order(fdr_res)]

length(which(results2$res <= 0.05)) #in 16/127 pairs , the lncRNA benefits from the signal from its neigboring gene 
length(which(results2$fdr_res <= 0.05)) #2/127 pairs, the lncRNA beneifts from signal from neighboring gene when accounting for multiple testing correction
length(which(results2$fdr_res2 <= 0.05)) #even after multiple testing correction, 124/127 pairs improve from adding lncRNA expression to PCG expression 

#saveRDS(results2, file="110_cis_antisense_pairs_survival_results_aug28.rds")
saveRDS(results2, file="127_cis_antisense_pairs_survival_results_10kb_nov16.rds")
write.csv(results2, file="127_cis_antisense_pairs_survival_results_10kb_nov16.csv", quote=F, row.names=F)

###START HERE###-----------------------------------------------------------------

#which of these are candidates? 
r = readRDS("127_cis_antisense_pairs_survival_results_10kb_nov16.rds")
cancs = rna[,c("type", "Cancer")]
cancs = unique(cancs)
colnames(cancs) = c("canc", "Cancer")

r = merge(r, cancs, by="canc")
r$combo = paste(r$lnc, r$Cancer, sep="_")

z = which(r$combo %in% allCands$combo)
cands_pairs = r[z,]

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

z = which(cands_pairs$pcg %in% census$ensg)
cands_pairs$census = ""
cands_pairs$census[z] = "YES"

cands_pairs$pcg = unlist(llply(cands_pairs$pcg, get_name_pcg))
cands_pairs$lnc = unlist(llply(cands_pairs$lnc, get_name))

prog_pcgs = cands_pairs

prog_pcgs$cor[prog_pcgs$rho <0] = "Negative"
prog_pcgs$cor[prog_pcgs$rho >0] = "Positive"
prog_pcgs$cor[prog_pcgs$rho_fdr > 0.05] = "NS"

prog_pcgs$rho_fdr[prog_pcgs$rho_fdr <0.0000000001] = 0.000001
prog_pcgs$lnc_pcg = paste(prog_pcgs$lnc, prog_pcgs$pcg, prog_pcgs$type, sep="/")

#add what kind of lncRNA it is
prog_pcgs$type_lnc = ""
for(i in 1:nrow(prog_pcgs)){
  lnc = prog_pcgs$combo[i]
  lnc = unlist(strsplit(lnc, "_"))[1]
  lnc_type = fantom$CAT_geneClass[which(fantom$CAT_geneID == lnc)]
  prog_pcgs$type_lnc[i] = lnc_type
}

prog_pcgs$diff_cindices = prog_pcgs$lncConcordance - prog_pcgs$pcgConcordance

pdf("figure2E_summary_lncs_pcgs_antisense_10kb_nov16.pdf", width=5,height=5)
g = ggplot(prog_pcgs, aes(pcgConcordance, lncConcordance, label=lnc_pcg)) +
 geom_point(color = "black")+
    scale_colour_manual(values = c("blue", "dimgrey", "red", "purple")) + 
    xlab("Neighbour PCG Concordance") + ylab("lncRNA Concordance") + 
    theme(legend.box = "horizontal", axis.text = element_text(size=13), 
      legend.text=element_text(size=10), legend.title=element_text(size=10))+
     xlim(0.45,0.75) + ylim(0.45,0.75) + geom_abline(intercept=0) +
     geom_vline(xintercept=0.5, linetype="dashed", color = "red") + 
     geom_hline(yintercept=0.5, linetype="dashed", color = "red")
     #geom_text_repel(data = subset(prog_pcgs, lncConcordance > 0.7 | pcgConcordance > 0.65), size=2, nudge_y = 0.1,
      #direction = "x",segment.color = "grey50",
      #segment.size = 0.05)
g
dev.off()

prog_pcgs$sig_spearman = ""
prog_pcgs$sig_spearman[prog_pcgs$rho_fdr < 0.05] = "sig"
prog_pcgs$lnc_better = ""
prog_pcgs$lnc_better[(prog_pcgs$fdr_res2 < 0.05) & (prog_pcgs$lncConcordance > prog_pcgs$pcgConcordance)] = "lnc_better"
prog_pcgs$lnc_better[is.na(prog_pcgs$lnc_better)] = ""
prog_pcgs$improv_perc = round((prog_pcgs$lncConcordance - prog_pcgs$pcgConcordance)/ prog_pcgs$pcgConcordance*100, digits=2)

pdf("/u/kisaev/figu2_lncs_pcgs_histogram_spearman_rho.pdf", width=6, height=5)
gghistogram(prog_pcgs, x = "rho",
   add = "mean", 
   color = "sig_spearman", fill = "sig_spearman",
   palette = c("grey", "black"))+theme_bw()+xlim(-0.5,1)
dev.off()

pdf("/u/kisaev/figu2_lncs_pcgs_histogram_cindex_improvement.pdf", width=6, height=5)
gghistogram(prog_pcgs, x = "improv_perc",
   add = "mean", 
   color = "sig_spearman", fill = "sig_spearman",
   palette = c("grey", "black"))+theme_bw()+xlim(-15,60)
dev.off()

pdf("/u/kisaev/figure2E_summary_lncs_pcgs_antisense_10kb_wSpearmanRho.pdf", width=5,height=5)
g = ggplot(prog_pcgs, aes(diff_cindices, rho, label=lnc_pcg)) +
 geom_point(aes(shape= sig_spearman, color=lnc_better), alpha=0.6, size=2)+
    scale_colour_manual(values = c("black", "blue", "red", "purple")) + 
    ylab("Spearman rho") + xlab("lncRNA cindex - PCG cindex") + 
    theme(legend.box = "horizontal", axis.text = element_text(size=13), 
      legend.text=element_text(size=10), legend.title=element_text(size=10))+
  geom_vline(xintercept=0, linetype="dashed", color = "red")+
      geom_hline(yintercept=0, linetype="dashed", color = "red")
ggpar(g, legend="bottom")
dev.off()

plot = as.data.table(filter(prog_pcgs, Cancer == "Brain Lower Grade Glioma"))
plot$sig_rho = ""
plot$sig_rho[plot$rho_fdr < 0.05] = "*"
plot$diff = plot$lncConcordance-plot$pcgConcordance
z = which((plot$diff > 0) & (plot$fdr_res2< 0.05))
plot$lncbetter = ""
plot$lncbetter[z] = "lnc"

t=as.data.table(table(plot$lnc))
t=t[order(-N)]
plot$lnc = factor(plot$lnc, levels=t$V1)

pdf("/u/kisaev/figure2X_LGG_lncs_pcgs_correlations.pdf", width=6, height=6)
gene_exp = ggplot(plot, aes(lnc, pcg)) +
  geom_tile(aes(fill = rho, colour =lncbetter, width=0.7, height=0.7), size=1) + theme_bw()
gene_exp = ggpar(gene_exp, x.text.angle = 45) + 
scale_fill_gradient(low = "grey", high = "red", na.value = 'white') + ylab("Protein-coding gene")+
xlab("lncRNA") + coord_equal()+ scale_color_manual(values=c("white", "black"))+
theme(legend.position="bottom") + geom_text(aes(label = sig_rho), size=4)+
theme(text = element_text(size=10)) 
gene_exp
dev.off()

#saveRDS(prog_pcgs, file="final_set_126_lncRNAPCG_pairs_nov16.rds")
#write.csv(prog_pcgs, file="126_cis_antisense_pairs_survival_results_10kb_nov16.csv", quote=F, row.names=F)

#make sure it's only updated list of candidates 
r = readRDS("final_set_126_lncRNAPCG_pairs_nov16.rds")
filter(r, combo %in% allCands$combo)

#127_cis_antisense_pairs_survival_results_10kb_nov16.rds
#lncRNA_cands_wPCGs_that_are_in_cis_10kb_nov16.rds

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

pdf("SUPP_FIGURE_BOXPLOT_lncRNAs_vs_NEARBY_PCGs.pdf")
g = ggplot(boxplot, aes(concordance_type, concordance)) +
geom_boxplot(aes(fill=concordance_type), outlier.colour="black", outlier.shape=16,
             outlier.size=2, notch=FALSE)+
    scale_colour_manual(values = c("blue", "dimgrey", "red", "purple")) + 
    xlab("Predictor of Survival") + ylab("Concordance") + 
    theme(legend.box = "horizontal", axis.text = element_text(size=13), 
      legend.text=element_text(size=10), legend.title=element_text(size=10)) + ylim(0,1)+ 
    geom_hline(yintercept=0.5, linetype="dashed", color = "red") 
g
dev.off()

