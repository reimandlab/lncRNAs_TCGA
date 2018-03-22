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

#1. OV
m1 = readRDS("cinds_OV_clinical_vs_all_lncRNAs_feb9_1000_runs.rds")
ov <- ldply(m1, data.frame)
colnames(ov) = c("Predictor", "cindex")
#2. just individual lncRNAs 
m2 = readRDS("cinds_OV_individual_lncRNAs_feb9_1000_runs.rds")
ov2 <- ldply(m2, data.frame)
colnames(ov2) = c("Predictor", "cindex")
ov_cinds_all = rbind(ov, ov2)
ov_cinds_all$canc = "OV"
preds = unique(ov_cinds_all$Predictor)
preds = preds[-(which(preds %in% c("just_clinical", "just_all_lncRNAs", "clinical_and_lncRNAs_all")))]
get_med = function(gene){
  med = median(ov_cinds_all$cindex[which(ov_cinds_all$Predictor ==  gene)])
  return(med)
}
preds_meds = llply(preds, get_med)
z = which(preds_meds == max(unlist(preds_meds)))
keep = c(preds[z])

#2. PAAD
m1 = readRDS("cinds_PAAD_clinical_vs_all_lncRNAs_feb9_1000_runs.rds")
paad <- ldply(m1, data.frame)
colnames(paad) = c("Predictor", "cindex")
#2. just individual lncRNAs 
m2 = readRDS("cinds_PAAD_individual_lncRNAs_feb9_1000_runs.rds")
paad2 <- ldply(m2, data.frame)
colnames(paad2) = c("Predictor", "cindex")
paad2_cinds_all = rbind(paad, paad2)
paad2_cinds_all$canc = "PAAD"
preds = unique(paad2_cinds_all$Predictor)
preds = preds[-(which(preds %in% c("just_clinical", "just_all_lncRNAs", "clinical_and_lncRNAs_all")))]
z = which(is.na(paad2_cinds_all$cindex))
paad2_cinds_all = paad2_cinds_all[-z,]

get_med = function(gene){
   med = median(paad2_cinds_all$cindex[which(paad2_cinds_all$Predictor ==  gene)])
  return(med)
}
preds_meds = llply(preds, get_med)
z = which(preds_meds == max(unlist(preds_meds)))
keep = c(keep,preds[z])

#3. LIHC
m1 = readRDS("cinds_LIHC_clinical_vs_all_lncRNAs_feb9_1000_runs.rds")
lihc <- ldply(m1, data.frame)
colnames(lihc) = c("Predictor", "cindex")
#2. just individual lncRNAs 
m2 = readRDS("cinds_LIHC_individual_lncRNAs_feb9_1000_runs.rds")
lihc2 <- ldply(m2, data.frame)
colnames(lihc2) = c("Predictor", "cindex")
lihc_cinds_all = rbind(lihc, lihc2)
lihc_cinds_all$canc = "LIHC"
preds = unique(lihc_cinds_all$Predictor)
preds = preds[-(which(preds %in% c("just_clinical", "just_all_lncRNAs", "clinical_and_lncRNAs_all")))]
z = which(is.na(lihc_cinds_all$cindex))
lihc_cinds_all = lihc_cinds_all[-z,]
get_med = function(gene){
  med = median(lihc_cinds_all$cindex[which(lihc_cinds_all$Predictor ==  gene)])
  return(med)
}
preds_meds = llply(preds, get_med)
z = which(preds_meds == max(unlist(preds_meds)))
keep = c(keep,preds[z])

#4. KIRC
m1 = readRDS("cinds_KIRC_clinical_vs_all_lncRNAs_feb9_1000_runs.rds")
kirc <- ldply(m1, data.frame)
colnames(kirc) = c("Predictor", "cindex")
#2. just individual lncRNAs 
m2 = readRDS("cinds_KIRC_individual_lncRNAs_feb9_1000_runs.rds")
kirc2 <- ldply(m2, data.frame)
colnames(kirc2) = c("Predictor", "cindex")
kirc_cinds_all = rbind(kirc, kirc2)
kirc_cinds_all$canc = "KIRC"
preds = unique(kirc_cinds_all$Predictor)
preds = preds[-(which(preds %in% c("just_clinical", "just_all_lncRNAs", "clinical_and_lncRNAs_all")))]
get_med = function(gene){
  med = median(kirc_cinds_all$cindex[which(kirc_cinds_all$Predictor ==  gene)])
  return(med)
}
preds_meds = llply(preds, get_med)
z = which(preds_meds == max(unlist(preds_meds)))
keep = c(keep,preds[z])

cinds_all = rbind(ov_cinds_all, paad2_cinds_all, lihc_cinds_all, kirc_cinds_all)
cinds_all = subset(cinds_all, Predictor %in% c("just_clinical", "just_all_lncRNAs", "clinical_and_lncRNAs_all", keep))

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

cinds_all$Predictor[cinds_all$Predictor == "just_clinical"] = "Clinical"
cinds_all$Predictor[cinds_all$Predictor == "just_all_lncRNAs"] = "CandidateLncRNAs"
cinds_all$Predictor[cinds_all$Predictor == "clinical_and_lncRNAs_all"] = "Combined"

pdf("CINDS_comparison_4cancers_facet_Feb16.pdf", width=16, height=12)
ov = ggboxplot(cinds_all[cinds_all$canc=="OV",], x="Predictor", y= "cindex", fill="Predictor", palette=mypal, 
	notch=F, main = "OV", 
	#order=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"), 
	legend="top", font.main = c(18, "bold", "black"),
          font.x = c(15, "plain", "black"),
          font.y = c(15, "plain", "black"),
          font.tickslab = c(15, "plain", "black"),
          font.legend = 15, ylab="C-index") + 
		stat_boxplot(geom = "errorbar", width = 0.5) +
		stat_compare_means(method = "wilcox.test",
                     ref.group = "Clinical", label = "p.signif")+
                     geom_hline(yintercept = 0.5, linetype = 2, colour="red")+ 
                     theme(axis.text.x = element_blank())
ov = ov + rremove("x.text")

lihc = ggboxplot(cinds_all[cinds_all$canc=="LIHC",], x="Predictor", y= "cindex", fill="Predictor", palette=mypal, 
  notch=F, 
  #order=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"), 
  legend="top", font.main = c(18, "bold", "black"),
  main = "LIHC", 
          font.x = c(15, "plain", "black"),
          font.y = c(15, "plain", "black"),
          font.tickslab = c(15, "plain", "black"),
          font.legend = 15, ylab="C-index") + 
    stat_boxplot(geom = "errorbar", width = 0.5) +
    stat_compare_means(method = "wilcox.test",
                     ref.group = "Clinical", label = "p.signif")+
                     geom_hline(yintercept = 0.5, linetype = 2, colour="red")+ 
                     theme(axis.text.x = element_blank())
lihc = lihc + rremove("x.text")

paad = ggboxplot(cinds_all[cinds_all$canc=="PAAD",], x="Predictor", y= "cindex", fill="Predictor", palette=mypal, 
  notch=F, main = "PAAD", 
  #order=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"), 
  legend="top", font.main = c(18, "bold", "black"),
          font.x = c(15, "plain", "black"),
          font.y = c(15, "plain", "black"),
          font.tickslab = c(15, "plain", "black"),
          font.legend = 15, ylab="C-index") + 
    stat_boxplot(geom = "errorbar", width = 0.5) +
    stat_compare_means(method = "wilcox.test",
                     ref.group = "Clinical", label = "p.signif")+
                     geom_hline(yintercept = 0.5, linetype = 2, colour="red")+ 
                     theme(axis.text.x = element_blank())
paad = paad + rremove("x.text")

kirc = ggboxplot(cinds_all[cinds_all$canc=="KIRC",], x="Predictor", y= "cindex", fill="Predictor", palette=mypal, 
  notch=F, main = "KIRC", 
  #order=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"), 
  legend="top", font.main = c(18, "bold", "black"),
          font.x = c(15, "plain", "black"),
          font.y = c(15, "plain", "black"),
          font.tickslab = c(15, "plain", "black"),
          font.legend = 15, ylab="C-index") + 
    stat_boxplot(geom = "errorbar", width = 0.5) +
    stat_compare_means(method = "wilcox.test",
                     ref.group = "Clinical", label = "p.signif")+
                     geom_hline(yintercept = 0.5, linetype = 2, colour="red")+ 
                     theme(axis.text.x = element_blank())
kirc = kirc + rremove("x.text")

#combine
g = ggarrange(ov, lihc, paad, kirc , 
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)

print(g)
dev.off()





























