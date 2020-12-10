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

#convert to ensgs
get_ensg = function(name){
  z = which(hg38$symbol == name)
  return(hg38$ensgene[z][1])
}
cis_pcgs$lnc = sapply(cis_pcgs$lnc, get_ensg)
cis_pcgs$pcg = sapply(cis_pcgs$pcg, get_ensg)
cis_pcgs$combo = paste(cis_pcgs$lnc, cis_pcgs$pcg, sep="_")

#check all lncRNA-cis PCG pairs
combos = unique(cis_pcgs$combo)

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
all = as.data.table(all)

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
  exp_dat = all[,..z]
  exp_dat = subset(exp_dat, Cancer == cancer)

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

  #train model using both of these genes and compare performance to just using lncRNA
  res = as.data.frame(res)
  res$res2 = res2
  res$lnc = get_name(lnc)
  res$pcg = get_name_pcg(pcgg)
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

  res$lnc_median = lnc_median
  res$pcg_median = pcg_median

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

  return(res)
}
}#end function


results = llply(combos, check_cis_pcg, .progress="text")
results2 = ldply(results)
results2 = as.data.table(results2)

#visualize for each pair
#pdf("/u/kisaev/Dec2020/lncRNA_PCG_correlations_nonlinearity.pdf", width=8)

#ggplot(results2, aes(y = factor(lnc),
#             x = factor(pcg))) +        ## global aes
#  geom_tile(aes(fill = rectheat)) +         ## to get the rect filled
#  geom_point(aes(colour = circlefill,
#                   size =circlesize))  +    ## geom_point for circle illusion
#  scale_color_gradient(low = "yellow",
#                       high = "red")+       ## color of the corresponding aes
#  scale_size(range = c(1, 20))+             ## to tune the size of circles
#  theme_bw()

#dev.off()

#visualize for each pair
mypal = c("#E5DFD9","#EAD286" ,"#D1EB7B", "#96897F" ,"#E5C0A6" ,
  "#72A93B", "#74DAE3" ,"#49B98D" ,"#D97B8F" ,"#70A2A4", "#64709B" ,"#DFBF38" ,"#61EA4F" ,
  "#C7CBE7", "#786DDA",
"#CFA0E0" ,"#67E9D0" ,"#7C9BE1", "#D94753" ,
"#AAE6B0", "#D13BDF" ,"#DEAEC7" ,"#BBE6DF" ,"#B2B47A" ,"#E6ECBA", "#C86ED7",
 "#7BEE95" ,"#6F46E6" ,"#65B9E0", "#C0EC3E",
"#DE8D54" ,"#DF4FA6")

results2$perc_lnc_off = as.numeric(results2$perc_lnc_off)
lncs_zero = unique(results2[,c("lnc", "perc_lnc_off", "canc")]) #79 unique genes

pdf("/u/kisaev/Dec2020/lncRNAs_wnearby_PCG_zero_expression_summary.pdf", width=7, height=5)

g = ggbarplot(lncs_zero, "lnc", "perc_lnc_off",
  fill = "canc", color = "black",
  palette = mypal)
g = ggpar(g, legend="none",
 font.tickslab = c(4,"plain", "black"),
 xtickslab.rt = 90)+ylab("% of patients with no detected expression")
g + facet_grid(~canc, scales = "free", space = "free") + theme(text = element_text(size=4))

dev.off()

#make histogram
#pdf("/u/kisaev/Dec2020/all_lncRNAs_zero_expression_summary_histogram.pdf", width=8, height=6)
#gghistogram(results2, x="perc_lnc_off", color="black", fill="#00AFBB")
#dev.off()
