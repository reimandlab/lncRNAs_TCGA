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

#pcgs = c("IDH1", "IDH2", "MGMT", "TERT", "ERBB2", "ESR1", "ATRX", "PGR",
 # "CDKN2A", "SETD2", "BAP1", "PBRM1", "PIK3CA", "ARID1A")

lncs = unique(allCands$gene_name)

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

check_cis_pcg = function(lnc){
  cancer = allCands$cancer[which(allCands$gene_name == lnc)]
  canc_type=canc_conv$type[canc_conv$Cancer == cancer]

  print(lnc)
  print(canc_type)

  if(canc_type == "LGG"){
    pcgs = c("IDH1", "IDH2", "MGMT", "TERT", "ATRX")
  }

  if(canc_type == "BRCA"){
    pcgs = c("ERBB2", "ESR1", "PGR")
  }

  if(canc_type == "KIRP"){
    pcgs = c("BAP1", "SETD2", "PBRM1")
  }

  if(canc_type == "STAD"){
    pcgs = c("CDKN2A", "PIK3CA", "ARID1A")
  }

  if(canc_type %in% c("LGG", "BRCA", "KIRP", "STAD")){

  for(i in 1:length(pcgs)){
  #get surv and exp data for these genes 
  pcgg  = pcgs[i]
  print(get_ensg_pcg(pcgg))
  pcgg_en = get_ensg_pcg(pcgg)
  lnc_en = get_ensg(lnc)
  z = which(colnames(all) %in% c(lnc_en, pcgg_en, "OS", "OS.time", "type", "patient", "Cancer"))
  exp_dat = all[,..z]
  exp_dat = subset(exp_dat, Cancer == cancer)
  
  #lncRNA 
  z = which(colnames(exp_dat) ==lnc_en)
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

  z = which(colnames(exp_dat) == pcgg_en)
  exp_dat[,z] = log1p(exp_dat[,..z])
  colnames(exp_dat)[z] = "PCG"

  z = which(colnames(exp_dat) == lnc_en)
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
  ggtitle(paste(lnc, pcgg, canc_conv$type[canc_conv$Cancer == cancer]))

  exp_dat$lnc_median = factor(exp_dat$lnc_median, levels=c(0,1))
  box = ggboxplot(exp_dat, x = "lnc_median", y = "PCG",
        color = "lnc_median", palette =c("#00AFBB", "#FC4E07"),
        add = "jitter") +  ggtitle(paste(lnc, pcgg, canc_conv$type[canc_conv$Cancer == cancer]))

  box = box + stat_compare_means() + stat_n_text()
  p = plot_grid(g, box, labels = c('A', 'B'), label_size = 12)
  print(p)

}}}

pdf("/u/kisaev/lncRNA_vs_biomarker_PCGs.pdf", width=10, height=9)
llply(lncs, check_cis_pcg, .progress="text")
dev.off()



