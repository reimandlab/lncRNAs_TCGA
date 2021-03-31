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

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = filter(allCands, data=="TCGA") #179 unique lncRNA-cancer combos, #166 unique lncRNAs
cands_dups = unique(allCands$gene[which(duplicated(allCands$gene))])
allCands$combo=paste(allCands$gene, allCands$cancer, sep="_")
allCands$gene_name = allCands$gene_symbol
combos = unique(allCands$combo)
cancers = unique(allCands$cancer)

get_pcgs_percentiles = function(canc){
  #get median of each pcg
  pcg_exp = filter(pcg, Cancer == canc)
  pcgs = colnames(pcg_exp)[which(str_detect(colnames(pcg_exp), "ENSG"))]

  get_median_pcg = function(gene){
    #print(gene)
    z = which(colnames(pcg_exp) == gene)
    median=median(unlist(pcg_exp[,..z]))
    if(median == 0){
       median = median(unlist(pcg_exp[,..z])[-(which(unlist(pcg_exp[,..z]) ==0))])
     }
    if(!(is.na(median))){
    return(c(gene, median))
  }}

  all_pcgs_medians = as.data.table(ldply(llply(pcgs, get_median_pcg, .progress="text")))
  colnames(all_pcgs_medians) = c("gene", "median")
  all_pcgs_medians$type_gene = "pcg"
  all_pcgs_medians$canc = canc
  return(all_pcgs_medians)
}

all_pcgs_medians_order = as.data.table(ldply(llply(cancers, get_pcgs_percentiles, .progress="text")))

get_exp_vs_pcgs = function(comb){

  print(comb)
  gene = unlist(strsplit(comb, "_"))[1]
  cancer = unlist(strsplit(comb, "_"))[2]
  print(cancer)

  #get lnc expression and all pcgs
  rna_exp = filter(rna, Cancer == cancer)
  z = which(colnames(rna_exp) == gene)
  lnc_median=median(unlist(rna_exp[,..z]))
  if(lnc_median == 0){
     lnc_median = median(unlist(rna_exp[,..z])[-(which(unlist(rna_exp[,..z]) ==0))])
  }

  lnc_row = as.data.frame(matrix(c(gene, lnc_median, "lnc", cancer), ncol=4))
  all_pcgs_medians = filter(all_pcgs_medians_order, canc == cancer)
  names(lnc_row) = colnames(all_pcgs_medians)
  all_pcgs_medians = rbind(all_pcgs_medians, lnc_row)
  all_pcgs_medians$median = as.numeric(all_pcgs_medians$median)
  all_pcgs_medians = all_pcgs_medians[order(median)]
  all_pcgs_medians$rank = 1:nrow(all_pcgs_medians)
  all_pcgs_medians$perc = all_pcgs_medians$rank/nrow(all_pcgs_medians)
  lnc_perc = filter(all_pcgs_medians, type_gene=="lnc")

  return(lnc_perc)
}

all_lncs_percentiles = as.data.table(ldply(llply(combos, get_exp_vs_pcgs, .progress="text")))
saveRDS(all_lncs_percentiles, file="all_lncs_percentiles.rds")

pdf("/u/kisaev/Jan2021/lncRNA_expression_vs_pcgs_in_same_cancer.pdf", width=6, height=5)
gghistogram(all_lncs_percentiles, x = "perc",
   color = "grey", fill = "grey",
   palette = c("grey", "black"))+theme_bw()+xlim(0,1)
dev.off()

all_lncs_percentiles$gene_name = ""
for(i in 1:nrow(all_lncs_percentiles)){
  print(i)
  gene_e=all_lncs_percentiles$gene[i]
  gene_name=allCands$gene_name[allCands$gene == gene_e][1]
  all_lncs_percentiles$gene_name[i] = gene_name
}

colnames(canc_conv)[2] = "canc"
all_lncs_percentiles = merge(all_lncs_percentiles, canc_conv, by="canc")
all_lncs_percentiles$combo = paste(all_lncs_percentiles$gene_name, all_lncs_percentiles$type)
all_lncs_percentiles = all_lncs_percentiles[order(perc)]
all_lncs_percentiles$combo = factor(all_lncs_percentiles$combo, levels=all_lncs_percentiles$combo)

#ordered barplot
pdf("/u/kisaev/Jan2021/lncRNA_expression_vs_pcgs_in_same_cancer_barplot.pdf", width=7, height=6)
ggplot(all_lncs_percentiles, aes(combo, perc)) +  geom_col(fill = "grey")+
				ylab("Percentile relative to median protein-coding gene expression")+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4))+
        geom_hline(yintercept=median(all_lncs_percentiles$perc), linetype="dashed", color = "black")
dev.off()

all_lncs_percentiles$test = "allcancers"
all_lncs_percentiles$gene_plot[all_lncs_percentiles$gene_name == "HOXA10-AS"] = "HOXA10-AS"

#add labelled point for HOXA10-AS
pdf("/u/kisaev/Jan2021/lncRNA_expression_vs_pcgs_in_same_cancer_violinplot.pdf", width=3, height=6)
ggplot(all_lncs_percentiles, aes(test, perc, label=gene_plot)) +
theme_bw()+
geom_violin()+
				ylab("Percentile relative to median protein-coding gene expression")+
        geom_jitter(aes(colour = type), shape=16) +
        colScale +
        geom_text_repel(data = subset(all_lncs_percentiles,
         gene_plot == c("HOXA10-AS")), min.segment.length = unit(0, 'lines'),
                      nudge_x = 0.1,   size=2)+
        theme(legend.position = "none")+
        #stat_summary(fun.y=median, geom="point", size=1, color="red")+
        #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=7))+
        xlab("All lncRNAs")
dev.off()
