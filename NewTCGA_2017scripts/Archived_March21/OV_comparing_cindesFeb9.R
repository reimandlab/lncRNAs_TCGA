###---------------------------------------------------------------
###Load libraries and data 
###---------------------------------------------------------------

source("source_code_Cox_MonteCarlo_CV_Jan12.R")
require(caTools)

#start with only lncRNA_intergenic
lincs = subset(fantom, CAT_geneClass == "lncRNA_intergenic")
z = which(colnames(rna) %in% lincs$gene)
rna = as.data.frame(rna)
rna = rna[,c(z, 5786:5790)]


###---------------------------------------------------------------
###Results - using just 202 patients 
###---------------------------------------------------------------

#1. 1000 CV, clinical variables combined + lncRNA candidates 
m1 = readRDS("cinds_OV_clinical_vs_all_lncRNAs_feb9_1000_runs.rds")
df <- ldply(m1, data.frame)
colnames(df) = c("Predictor", "cindex")

#2. just individual lncRNAs 
m2 = readRDS("cinds_OV_individual_lncRNAs_feb9_1000_runs.rds")
df2 <- ldply(m2, data.frame)
colnames(df2) = c("Predictor", "cindex")

cinds_all = rbind(df, df2)
#cinds_all$Predictor = factor(cinds_all$Predictor, levels=c("Clinical",
	#"Combined_lncRNAs","Combined_lncRNAs_and_Clinical", "RP13-188A5.1", "AC018647.3", "AP001057.1", "RP11-555H7.2",
	# "RP11-321E2.4","RP11-443B7.3", "U3"))

###---------------------------------------------------------------
###Plot C-Indeces for the different models 
###---------------------------------------------------------------

#change gene ids to names 
for(i in 1:nrow(cinds_all)){
  z = which(fantom$gene == cinds_all$Predictor[i])
  if(!(length(z)==0)){
    cinds_all$Predictor[i] = fantom$CAT_geneName[z]
  }
}

pall = c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", 
  "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", 
  "#DDAA77", "#771122", "#AA4455", "#DD7788")

pdf("MODEL_comparison_OV_FEB9.pdf", width=16, height=12)
g = ggboxplot(cinds_all, x="Predictor", y= "cindex", fill="Predictor", palette=pall, 
	notch=F, main="305 Ovarian Cancer Patients, 1000 Cross Validations", 
	#order=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"), 
	legend="right", font.main = c(18, "bold", "black"),
          font.x = c(15, "plain", "black"),
          font.y = c(15, "plain", "black"),
          font.tickslab = c(15, "plain", "black"),
          font.legend = 15, ylab="C-index") + 
		stat_boxplot(geom = "errorbar", width = 0.5) +
		stat_compare_means(method = "wilcox.test",
                     ref.group = "just_clinical", label = "p.signif")+
                     geom_hline(yintercept = 0.5, linetype = 2, colour="red")+ 
                     theme(axis.text.x = element_blank())
                     
g = g + rremove("x.text")
print(g)
dev.off()





























