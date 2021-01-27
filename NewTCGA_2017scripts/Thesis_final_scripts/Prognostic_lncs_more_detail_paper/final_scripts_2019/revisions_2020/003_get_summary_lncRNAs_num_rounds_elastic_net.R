set.seed(911)

source("/u/kisaev/lncRNAs_TCGA/NewTCGA_2017scripts/Thesis_final_scripts/Prognostic_lncs_more_detail_paper/final_scripts_2019/revisions_2020/load_data.R")
#setWD
setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNAs_2019_manuscript")

#get candidates files

#Data--------------------------------------------
allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = filter(allCands, data=="TCGA") #142 unique lncRNA-cancer combos

#results from elastic net
cands = readRDS("lncRNAs_selected_by_EN_april14.rds") #1000 runs of cross-validations using new updated dataset (GBM=150, OV and LUAD)
cands$combo = paste(cands$gene, cands$cancer, sep="_")
cands = as.data.table(filter(cands, combo %in% allCands$combo))

get_name = function(gene){
    z = which(allCands$gene == gene)
    return(allCands$gene_symbol[z][1])
}

#get gene names
cands$gene_symbol = sapply(cands$gene, get_name)

#get order of cancers
t= as.data.table(table(cands$type))
t=t[order(-N)]
cands$type = factor(cands$type, levels=unique(t$V1))

cands=cands[order(-num_rounds_selected)]
cands$combo = paste(cands$gene_symbol, cands$type)
cands$combo = factor(cands$combo, levels=unique(cands$combo))

#make barplot summary
pdf("/u/kisaev/Jan2021/suppFigure2_num_rounds_each_lncRNA_selected.pdf", width=14, height=6)

g = ggplot(data=cands, aes(x=combo, y=num_rounds_selected, order = -num_rounds_selected)) +
  geom_bar(stat="identity", fill="grey",colour="black") +
  facet_grid(~ type, scale="free", space = "free")+
 theme_bw()
ggpar(g, xtickslab.rt=90, font.tickslab=c(6, "plain", "black")) +
geom_hline(yintercept=500, color = "red") + ylim(0,1000)+
  theme(strip.text.x = element_text(size = 9, colour = "Black", angle=90), axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1)) +
  xlab("lncRNA") + ylab("Number of rounds selected in Elastic Net cross-validations")

dev.off()
