#------------------------------------------------------------------------------
#This script takes all the PCGs that are significantly differentially expressed
#between lncRNA risk groups and calculates their correlation 
#plots heatmap of this data as well as lncRNA and pcg prognostic
#relationship 
#------------------------------------------------------------------------------

#source code
source("check_lnc_exp_cancers.R")

#COSMIC cancer gene census
census = read.csv("Census_allFri_Jul_13_16_55_59_2018.csv")

#------FEATURES-----------------------------------------------------

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = subset(allCands, data == "TCGA") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 
allCands$combo = unique(paste(allCands$gene, allCands$cancer, sep="_"))

val_cands = read.csv("175_lncRNA_cancers_combos_23_cancer_types_july5.csv")
val_cands = as.data.table(val_cands)
val_cands = subset(val_cands, data == "PCAWG") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 
val_cands$combo = unique(paste(val_cands$gene, val_cands$cancer, sep="_"))
val_cands = subset(val_cands, top_pcawg_val == "YES") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 

#Combined into one dataframe because need to get ranks 
all <- merge(rna, pcg, by = c("patient", "Cancer"))
all = all[,1:25170]

#canc conversion
canc_conv = rna[,c(which(colnames(rna) %in% c("Cancer", "type")))]
canc_conv = canc_conv[!duplicated(canc_conv),]
colnames(canc_conv) = c("type", "cancer")

#------DATA-----------------------------------------------------

#lncRNAs with 5% improvemenet over clinical 
lncs_perc = readRDS("148_combos_robust_5perc_increase_internal_validation_survival_lncRNAs_aug9.rds")

#-----FIGURE 2E-------------------------------------------------
#summary of how many individual lncRNAs were selected in each cancer type
#these are the ones that were studied throughout the paper
#show how many of them had a 5% improvement over clinical 
allCands = filter(allCands, data == "TCGA")
head(lncs_perc)
allCands$better_than_clin = ""
z = which(allCands$combo %in% lncs_perc$combo)
allCands$better_than_clin[z] = "yes"
allCands$better_than_clin[-z] = "no"
allCands$multivar_sig = ""
z = which(allCands$fdr <= 0.05)
allCands$multivar_sig[z] = "yes"
allCands$multivar_sig[-z] = "no"
allCands = merge(allCands, canc_conv, by = c("cancer"))

#get summary 
summ = as.data.table(table(allCands$type, allCands$better_than_clin))
colnames(summ) = c("Cancer", "BetterThanClin", "NumCandidates")
summ = summ[order(-BetterThanClin, -NumCandidates)]

#order of cancer types
od = as.data.table(table(allCands$type))
od = od[order(-N)]

summ$Cancer = factor(summ$Cancer, levels = unique(od$V1))
summ$BetterThanClin = factor(summ$BetterThanClin, levels = c("yes", "no"))

#barplot summary
# Create the barplot
pdf("figure2E_aug13.pdf")
ggplot(data=summ, aes(x=Cancer, y=NumCandidates, fill=BetterThanClin)) +
  theme_bw()+
  geom_bar(stat="identity", color="black")+ 
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  scale_fill_manual(values=c('#E69F00', '#999999')) + 
  xlab("Cancer") + ylab("# of lncRNAs selected in more than 50% of cross validations")
  dev.off()

