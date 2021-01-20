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

get_name=function(g){
  z=which(allCands$gene == g)
  name=allCands$gene_symbol[z]
  name=name[1]
  return(name)
}

#setWD
setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNAs_2019_manuscript")

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = subset(allCands, data == "TCGA") #173 unique lncRNA-cancer combos, #166 unique lncRNAs
allCands$combo = unique(paste(allCands$gene, allCands$cancer, sep="_"))

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

#pcgs = c("IDH1", "IDH2", "MGMT", "TERT", "ERBB2", "ESR1", "ATRX", "PGR",
 # "CDKN2A", "SETD2", "BAP1", "PBRM1", "PIK3CA", "ARID1A")

lncs = unique(allCands$gene_symbol)

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

check_cis_pcg = function(lnc){
  cancer = allCands$cancer[which(allCands$gene_symbol == lnc)]
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

  pcg_res = as.data.frame(matrix(ncol=5))
  colnames(pcg_res)=c("lnc", "pcg", "canc", "spear_rho", "spear_p")

  for(i in 1:length(pcgs)){
  #get surv and exp data for these genes
  pcgg  = pcgs[i]
  print(pcgg)
  print(get_ensg_pcg(pcgg))
  pcgg_en = get_ensg_pcg(pcgg)
  lnc_en = allCands$gene[allCands$gene_symbol == lnc]
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
  exp_dat=as.data.frame(exp_dat)
  exp_dat[,z] = log1p(exp_dat[,z])
  colnames(exp_dat)[z] = "PCG"

  z = which(colnames(exp_dat) == lnc_en)
  exp_dat[,z] = log1p(exp_dat[,z])
  colnames(exp_dat)[z] = "LNC"

  #get correlation between them (1)
  rho = rcorr(exp_dat$LNC, exp_dat$PCG, type="spearman")$r[2]
  rho_p = rcorr(exp_dat$LNC, exp_dat$PCG, type="spearman")$P[2]

  lnc_name=lnc
  if(lnc_name=="HOXA-AS4"){
    lnc_name="HOXA10-AS"
  }
  g =  ggscatter(exp_dat, x = "LNC", y = "PCG",
          color = "black", shape = 21, size = 3, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "spearman", label.sep = "\n"))+
  ggtitle(paste(lnc_name, pcgg, canc_conv$type[canc_conv$Cancer == cancer]))

  exp_dat$lnc_median = factor(exp_dat$lnc_median, levels=c(0,1))
  box = ggboxplot(exp_dat, x = "lnc_median", y = "PCG",
        color = "lnc_median", palette =c("#00AFBB", "#FC4E07"),
        add = "jitter") +  ggtitle(paste(lnc_name, pcgg, canc_conv$type[canc_conv$Cancer == cancer]))

  box = box + stat_compare_means() + stat_n_text()
  p = plot_grid(g, box, labels = c('A', 'B'), label_size = 12)
  print(p)
  canc=canc_conv$type[canc_conv$Cancer==cancer]
  res=c(lnc, pcgg, canc, rho, rho_p)
  names(res)=c("lnc", "pcg", "canc", "spear_rho", "spear_p")
  pcg_res = rbind(pcg_res, res)
}
pcg_res=pcg_res[-1,]
return(pcg_res)

}
}

pdf("/u/kisaev/Jan2021/lncRNA_vs_biomarker_PCGs.pdf", width=10, height=9)
all_res = as.data.table(ldply(llply(lncs, check_cis_pcg, .progress="text")))
dev.off()

#adjust p-values
all_res$fdr = p.adjust(all_res$spear_p, method="fdr")
all_res$sig=""
all_res$sig[all_res$fdr < 0.05] = "*"
all_res$lnc[all_res$lnc=="HOXA-AS4"] = "HOXA10-AS"
all_res$spear_rho = as.numeric(all_res$spear_rho)
all_res$canc=factor(all_res$canc, levels=c("BRCA", "LGG", "STAD", "KIRP"))

#make summary plot
pdf("/u/kisaev/Jan2021/lncRNA_vs_biomarker_PCGs_geom_tile_plot_summary.pdf", height=6)
ggplot(all_res, aes(lnc, pcg)) +
  geom_tile(aes(fill = spear_rho, width=0.7, height=0.7), size=0.55, color="grey") +
  theme_bw() + geom_text(aes(label = sig), size=3) +
  theme(legend.title=element_blank(), legend.position="bottom", axis.title.x=element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size=6)) +
    scale_fill_gradientn(colours = c("blue", "white", "red"),
                       values = scales::rescale(c(-0.75, -0.25, 0, 0.25, 0.75)))+
   #scale_fill_gradient(low = "blue", high = "red")+
    facet_grid(cols = vars(canc), scales = "free", space = "free")+
     theme(strip.background = element_rect(colour="black", fill="white",
                                       size=1.5, linetype="solid"))
dev.off()
