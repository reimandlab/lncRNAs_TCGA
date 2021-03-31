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

cancers= unique(rna$type)

get_lnc_percentiles = function(canc){

  print(canc)
  #get median of each pcg
  rna_exp = filter(rna, type == canc)
  lncrnas = colnames(rna_exp)[which(str_detect(colnames(rna_exp), "ENSG"))]

  get_median_lnc = function(gene){
    #print(gene)
    z = which(colnames(rna_exp) == gene)
    median=median(unlist(rna_exp[,..z]))
    if(median == 0){
       median = median(unlist(rna_exp[,..z])[-(which(unlist(rna_exp[,..z]) ==0))])
     }
    if(!(is.na(median))){
    return(c(gene, median))
  }}

  all_lncs_medians = as.data.table(ldply(llply(lncrnas, get_median_lnc, .progress="text")))
  colnames(all_lncs_medians) = c("gene", "median")
  all_lncs_medians$type_gene = "lncRNA"
  all_lncs_medians$canc = canc
  print(head(all_lncs_medians))
  return(all_lncs_medians)
}

all_lncs_medians_order = as.data.table(ldply(llply(cancers, get_lnc_percentiles, .progress="text")))

#make violin plots of all lncNRA median expression across each cancer
all_lncs_medians_order$median = as.numeric(all_lncs_medians_order$median)
meds = as.data.table(all_lncs_medians_order %>% group_by(canc) %>% dplyr::summarize(med=median(median)))
meds = meds[order(med)]
all_lncs_medians_order$canc = factor(all_lncs_medians_order$canc, levels=meds$canc)
all_lncs_medians_order$median = log1p(all_lncs_medians_order$median)

#ordered barplot
pdf("/u/kisaev/Jan2021/lncRNA_expression_violin_plot_all_cancers.pdf", width=7, height=6)
ggplot(all_lncs_medians_order, aes(canc, median)) +
geom_violin()+
				ylab("Median lncRNA log1p(FPKM-UQ)")+
        stat_summary(fun.y=median, geom="point", size=1, color="red")+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=7))+xlab("Cancer")
dev.off()
