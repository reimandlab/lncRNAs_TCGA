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
##cindeces
m1 = readRDS("MODEL1_305_OV_pats_CAND_genes_AND_STAGE_GRADE_AGE_CV_1000timesJan26.RDS")
m1 = as.data.frame(m1)
colnames(m1) = "cindex"
m1$Predictor = "Combined_lncRNAs_and_Clinical"
m1$Model = 1

#2. 1000 CV, all 7 cands combined 
##cindeces
m2 = readRDS("MODEL2_305_OV_pats_CAND_genes_CV_1000timesJan26.RDS")
m2 = as.data.frame(m2)
colnames(m2) = "cindex"
m2$Predictor = "Combined_lncRNAs"
m2$Model = 2

#3. 1000 CV, just clinical combined 
##cindeces
m3 = readRDS("MODEL3_305_OV_pats_grade_stage_age_CV_1000timesJan26.RDS")
m3 = as.data.frame(m3)
colnames(m3) = "cindex"
m3$Predictor = "Clinical"
m3$Model = 3

#4. 1000 CV, just lncRNA 1 
##cindeces
m4 = readRDS("MODEL4_305_OV_pats_CAND_genes_CV_1000timesJan26.RDS")
m4 = as.data.frame(m4)
colnames(m4) = "cindex"
m4$Predictor = "RP13-188A5.1"
m4$Model = 4

#5. 1000 CV, just lncRNA 2 
##cindeces
m5 = readRDS("MODEL5_305_OV_pats_CAND_genes_CV_1000timesJan26.RDS")
m5 = as.data.frame(m5)
colnames(m5) = "cindex"
m5$Predictor = "AC018647.3"
m5$Model = 5

#6. 1000 CV, just lncRNA 3
##cindeces
m6 = readRDS("MODEL6_305_OV_pats_CAND_genes_CV_1000timesJan26.RDS")
m6 = as.data.frame(m6)
colnames(m6) = "cindex"
m6$Predictor = "AP001057.1"
m6$Model = 6

#7. 1000 CV, just lncRNA 4
##cindeces
m7 = readRDS("MODEL7_305_OV_pats_CAND_genes_CV_1000timesJan26.RDS")
m7 = as.data.frame(m7)
colnames(m7) = "cindex"
m7$Predictor = "RP11-555H7.2"
m7$Model = 7

#8. 1000 CV, just lncRNA 5
##cindeces
m8 = readRDS("MODEL8_305_OV_pats_CAND_genes_CV_1000timesJan26.RDS")
m8 = as.data.frame(m8)
colnames(m8) = "cindex"
m8$Predictor = "RP11-321E2.4"
m8$Model = 8

#9. 1000 CV, just lncRNA 6
##cindeces
m9 = readRDS("MODEL9_305_OV_pats_CAND_genes_CV_1000timesJan26.RDS")
m9 = as.data.frame(m9)
colnames(m9) = "cindex"
m9$Predictor = "RP11-443B7.3"
m9$Model = 9

#10. 1000 CV, just lncRNA 7
##cindeces
m10 = readRDS("MODEL10_305_OV_pats_CAND_genes_CV_1000timesJan26.RDS")
m10 = as.data.frame(m10)
colnames(m10) = "cindex"
m10$Predictor = "U3"
m10$Model = 10

cinds_all = rbind(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10)
#cinds_all$Predictor = factor(cinds_all$Predictor, levels=c("Clinical",
	"Combined_lncRNAs","Combined_lncRNAs_and_Clinical", "RP13-188A5.1", "AC018647.3", "AP001057.1", "RP11-555H7.2",
	 "RP11-321E2.4","RP11-443B7.3", "U3"))

###---------------------------------------------------------------
###Plot C-Indeces for the different models 
###---------------------------------------------------------------

pdf("305_ovarian_cancer_patients_comparing_model_performanceJan26.pdf", width=16, height=12)
g = ggboxplot(cinds_all, x="Model", y= "cindex", fill="Predictor", palette=mypal, 
	notch=F, main="305 Ovarian Cancer Patients, 1000 Cross Validations", 
	#order=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"), 
	legend="right", font.main = c(18, "bold", "black"),
          font.x = c(15, "plain", "black"),
          font.y = c(15, "plain", "black"),
          font.tickslab = c(15, "plain", "black"),
          font.legend = 15)+ 
		stat_boxplot(geom = "errorbar", width = 0.5) +
		stat_compare_means(method = "wilcox.test",
                     ref.group = "3", label = "p.format")+
                     geom_hline(yintercept = 0.5, linetype = 2, colour="red")
ggpar(g, xlab ="Model", ylab="C-index")
print(g)
dev.off()



























