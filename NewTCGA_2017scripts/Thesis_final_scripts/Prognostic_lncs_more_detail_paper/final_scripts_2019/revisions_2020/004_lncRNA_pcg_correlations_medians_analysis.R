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
  res_all = as.data.frame(matrix(ncol=13, nrow=1))
  colnames(res_all) = c("res", "res2", "lnc", "pcg",
 "canc", "pcg_improvelnc" , "lnc_improvepcg" ,"rho",
  "rho_p" , "lnc_median" , "pcg_median" , "perc_lnc_off",
 "perc_pcg_off")
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


    #km plot just lncRNA
    gene = lnc
    cancer = exp_dat$type[1]

    #km plot just pcgg
    gene = pcgg
    cancer = exp_dat$type[1]

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
    res$lnc = allCands$gene_symbol[allCands$gene==lnc][1]
    res$pcg = get_name_pcg(pcgg)
    res$canc = exp_dat$type[1]

    res$pcg_improvelnc = ""
    res$pcg_improvelnc[res$res <= 0.05] = "yes"
    res$pcg_improvelnc[res$res > 0.05] = "no"

    res$lnc_improvepcg = ""
    res$lnc_improvepcg[res$res2 <= 0.05] = "yes"
    res$lnc_improvepcg[res$res2 > 0.05] = "no"

    #res$pcg_fdr_pval = pcg_pvalue #Fdr from cox model
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
    res_all = rbind(res_all, res)
    #print(res_all)

    }

  }#for(i in 1:length(cancer))
  res_all = res_all[-1,]
  return(res_all)

}#end function

results = llply(combos, check_cis_pcg, .progress="text")
results2 = ldply(results)
results2 = as.data.table(results2)

#visualize for each pair
#add colours for each cancer
colnames(colours_palette)=c("canc", "color")
results2$perc_lnc_off = as.numeric(results2$perc_lnc_off)
lncs_zero = unique(results2[,c("lnc", "perc_lnc_off", "canc")]) #92 instances
lncs_zero = merge(lncs_zero, colours_palette, by="canc")
t=as.data.table(table(lncs_zero$canc))
t=t[order(-N)]
lncs_zero$canc = factor(lncs_zero$canc, levels=t$V1)
colnames(t)[1] = "canc"
t=merge(t, colours_palette, by="canc")
t=t[order(-N)]
lncs_zero = lncs_zero[order(-perc_lnc_off)]
lncs_zero$lncRNA = paste(lncs_zero$lnc, lncs_zero$canc)

pdf("/u/kisaev/Jan2021/lncRNAs_wnearby_PCG_zero_expression_summary.pdf", width=8, height=4)

g = ggbarplot(lncs_zero, "lncRNA", "perc_lnc_off",
  fill = "canc", color = "black", sort.val = "desc",          # Sort the value in descending order
  palette = unique(t$color))
g = ggpar(g, legend="none",
 font.tickslab = c(4,"plain", "black"),
 xtickslab.rt = 90)+ylab("% of patients with no detected expression")+
 geom_hline(yintercept=0.9, linetype="dashed", color = "grey")+ylim(c(0,1))
g + facet_grid(~canc, scales = "free", space = "free") + theme(text = element_text(size=4))

dev.off()

#make histogram
#pdf("/u/kisaev/Dec2020/all_lncRNAs_zero_expression_summary_histogram.pdf", width=8, height=6)
#gghistogram(results2, x="perc_lnc_off", color="black", fill="#00AFBB")
#dev.off()
