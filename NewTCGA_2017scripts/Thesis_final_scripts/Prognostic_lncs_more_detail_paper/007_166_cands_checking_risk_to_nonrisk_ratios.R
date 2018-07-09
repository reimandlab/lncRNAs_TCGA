library(survAUC)
#source("source_code_Cox_MonteCarlo_CV_Mar13.R")
require(caTools)
#check if this person is in my analysis: TCGA-61-2095
library(glmnet)
library(survcomp)
library(caret)
library(stringr)
library(EnvStats)
library(patchwork)

source("universal_LASSO_survival_script.R")

library(ggpubr)
library(ggrepel)
library(viridis)
library(patchwork)
library(caret)  
library(Rtsne)

#------FEATURES-----------------------------------------------------

cands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_May4.rds")
#cands = filter(cands, data == "PCAWG", pval <=0.05)
cands = filter(cands, AnalysisType == "noFDR")
#colnames(cands)[7] = "canc"
cands$Cancer = NULL
all_cands = cands


#--------This script ------------------------------------------------

#just make KM plots for TCGA 
#whatever data is available for PCAWG
#make them KM plots as well 
#just get list of genes that are significant in both data sets
#also check Cox PH assumptions within each data-set

#--------------------------------------------------------------------

#write function that adds tag to whole data group 
#and does survival analysis on whole group

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = filter(allCands, data=="TCGA") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 

cands_dups = unique(allCands$gene[which(duplicated(allCands$gene))])

z = which(cancers %in% allCands$cancer)
cancer_data = canc_datas[z] #cancers list and canc_datas list should be the same 

get_canc_data_for_plot = function(dtt){
  #get cancer specific candidates 
  z = which(colnames(dtt) %in% c(as.character(allCands$gene[allCands$cancer == dtt$Cancer[1]]), "age_at_initial_pathologic_diagnosis", 
    "OS.time", "OS", "gender", "race", "patient", "clinical_stage", "histological_grade", "treatment_outcome_first_course", 
    "new_tumor_event_type", "Cancer"))
  dtt = dtt[,z]
  return(dtt)
}
filtered_data = llply(cancer_data, get_canc_data_for_plot)

add_tags = function(dtt){
  print(dtt$Cancer[1])
  rownames(dtt) = dtt$patient
  dtt$patient = NULL
  #get ratio of #of risk patients to #non risk patients 
  #if 50/50 (when median split) then ratio should be 1 or log2(1) = 0
  #log1p 
  z = which(str_detect(colnames(dtt), "ENSG"))
  lncs = colnames(dtt)[z]
  get_ratio = function(lnc){
    lnc_data = dtt[,which(colnames(dtt) %in% c(lnc, "Cancer"))] 
      median2 <- quantile(as.numeric(lnc_data[,1]), 0.5)
      lnc_data$median =""
       if(median2 ==0){
        #if median = 0 then anyone greater than zero is 1 
        l1 = which(lnc_data[,1] > 0)
        l2 = which(lnc_data[,1] ==0)
        lnc_data$median[l1] = 1
        lnc_data$median[l2] = 0
        }

      if(!(median2 ==0)){
        l1 = which(lnc_data[,1] >= median2)
        l2 = which(lnc_data[,1] < median2)
        lnc_data$median[l1] = 1
        lnc_data$median[l2] = 0
        }

      #which one is high risk --> need surivval data
      lnc_data$patient = rownames(lnc_data)
    
      lnc_data$median[lnc_data$median ==0] = "Low"
      lnc_data$median[lnc_data$median==1] = "High"

      #cox ph
      z = which((allCands$gene == lnc) & (allCands$cancer == lnc_data$Cancer[1]))

      HR = as.numeric(allCands$HR[z])
      lnc_data$risk =""
      if(HR <1){
        risk = "Low"
        lnc_data$risk = ""
        lnc_data$risk[lnc_data$median=="High"] ="noRISK"
        lnc_data$risk[lnc_data$median=="Low"] ="RISK"
      }
      if(HR >1){
        risk = "High"
        lnc_data$risk = ""
        lnc_data$risk[lnc_data$median=="High"] ="RISK"
        lnc_data$risk[lnc_data$median=="Low"] ="noRISK"
      }

    order = c("noRISK", "RISK")
    lnc_data$risk <- factor(lnc_data$risk, levels = order)
    lnc_data[,1] = log1p(lnc_data[,1])
    colnames(lnc_data)[1] = "lncRNA"
    ratio = round(length(which(lnc_data$risk == "RISK"))/length(which(lnc_data$risk=="noRISK")), digits=3)
    z = which((allCands$gene == lnc) & (allCands$cancer == lnc_data$Cancer[1]))
    pval = round(as.numeric(allCands$fdr_pval[z]), digits=5)
    lowci = round(as.numeric(allCands$low95[z]), digits=5)
    highci = round(as.numeric(allCands$upper95[z]), digits=5)

    g = ggboxplot(lnc_data, x = "risk", y="lncRNA", fill="risk", palette=mypal[c(2,1)]) + stat_n_text(size=5) + theme_bw() + 
    stat_compare_means(method = "wilcox.test") + 
    ggtitle(paste(lnc, lnc_data$Cancer[1], "\nratio is", ratio, "\np-val is", pval)) + ylab("log1p(FPKM)")
    print(g)
    row = c(lnc, lnc_data$Cancer[1], HR, pval, ratio, lowci, highci)
    names(row) = c("lncRNA", "Cancer", "HR", "pval", "ratio", "lowci", "highci")
    return(row)
  }
  lnc_data = llply(lncs, get_ratio, .progress="text")
  lnc_data = do.call(rbind.data.frame, lnc_data)
  colnames(lnc_data) = c("lncRNA", "Cancer", "HR", "pval", "ratio", "lowci", "highci")
  return(lnc_data)
}    

pdf("166_candidates_ratio_risk_tononrisk_TCGA_june28.pdf")
filtered_data_tagged = llply(filtered_data, add_tags, .progress="text")
dev.off()

#----------------- PART 2 ------------------------------------------------

#Summarize the relationship between ratio, p-value and Hazard ratio 
#X-AXIS is lncRNAs..., have covariate for cancer type.... 

filtered_data_tagged = do.call(rbind.data.frame, filtered_data_tagged)
saveRDS(filtered_data_tagged, file="007_166_cands_checking_risk_to_nonrisk_ratios_data_june28.rds")


#START----------------------------------------------------------------------------------------------

filtered_data_tagged = readRDS("007_166_cands_checking_risk_to_nonrisk_ratios_data_june28.rds")

colnames(filtered_data_tagged) = c("lncRNA", "Cancer", "HR", "pval", "ratio", "lowci", "highci")
filtered_data_tagged$HR = as.numeric(filtered_data_tagged$HR)
filtered_data_tagged$pval = as.numeric(filtered_data_tagged$pval)
filtered_data_tagged$lowci = as.numeric(filtered_data_tagged$lowci)
filtered_data_tagged$highci = as.numeric(filtered_data_tagged$highci)
filtered_data_tagged$HR = log2(filtered_data_tagged$HR)
filtered_data_tagged$lowci = log2(filtered_data_tagged$lowci)
filtered_data_tagged$highci = log2(filtered_data_tagged$highci)
filtered_data_tagged$sig = ""
filtered_data_tagged$sig[filtered_data_tagged$pval <= 0.05] = "Sig"
filtered_data_tagged$sig[filtered_data_tagged$pval > 0.05] = "NotSig"

##--------Check which ones significant/were tested in PCAWG---------------------------
all_results = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
all_results = all_results[!duplicated(all_results), ]
pcawg_tested = filter(all_results, data == "PCAWG")
pcawg_tested = as.data.table(pcawg_tested)
filtered_data_tagged$pcawg_pval = ""
filtered_data_tagged$pcawg_match = ""

all_results = filter(all_results, pval <= 0.05)
lncs = as.list(as.character(unique(all_results$gene[all_results$data == "PCAWG"])))
check_match = function(lnc){
  z = which(all_results$gene == lnc)
  res = as.data.table(all_results[z,])
  
  if(dim(res)[1] > 2){
  
  canc = res$cancer[which(duplicated(res$cancer))]
  res = filter(res, cancer == canc)
  }

  test = which(as.numeric(res$HR) >= 1)
  test = length(test)
  if(test ==2){
    match = "match"
  }
  if(test ==0){
    match = "match"
  }
  if(test ==1){
    match = "nomatch"
  }
  canc = unique(res$canc)
  return(c(lnc, match, canc))
}

matches = llply(lncs, check_match)
matches = do.call(rbind.data.frame, matches)
matches = as.data.table(matches)
colnames(matches) = c("lnc", "match", "cancer")
matches = filter(matches, match == "match")
colnames(matches)[1] = "gene"
matches = merge(matches, all_results, by=c("gene", "cancer"))

#remove weird KIRC one 
z = which(matches$gene == "ENSG00000250360")
matches = matches[-z,]
matches = as.data.table(matches)

for(i in 1:nrow(filtered_data_tagged)){
  canc = filtered_data_tagged$Cancer[i]
  lnc = filtered_data_tagged$lncRNA[i]
  new_dat = filter(pcawg_tested, gene == lnc, cancer == canc)
  pcawg_dat = filter(matches, gene == lnc, cancer == canc)
  if(dim(new_dat)[1] >=1){
    pval = new_dat$pval
    filtered_data_tagged$pcawg_pval[i] = pval
    if(dim(pcawg_dat)[1] >=1){
      filtered_data_tagged$pcawg_match[i] = "MATCH"
    }
  }
}


#convert duplicated gene id to unique id 
if(!(length(cands_dups)==0)){
  z = which(filtered_data_tagged$lncRNA %in% cands_dups)
  for(y in 1:length(z)){
    filtered_data_tagged$lncRNA[z][y] = paste(filtered_data_tagged$lncRNA[z][y], y, sep="_")
  } 
}


hr = filtered_data_tagged[,c(1:2, 3, 4, 5, 8,9, 10)]
hr$hr_type = "HR"
CI = filtered_data_tagged[,c(1:2, 6, 4, 5, 8,9, 10)] ; colnames(CI)[3] = "HR"
CI$hr_type = "CI"
CI2 = filtered_data_tagged[,c(1:2, 7, 4, 5, 8,9, 10)] ; colnames(CI2)[3] = "HR"
CI2$hr_type = "CI"

hr = as.data.table(hr)
hr = hr[order(HR)]
orderhr= (hr$lncRNA)

filtered_data_tagged = rbind(hr, CI, CI2)
filtered_data_tagged = as.data.table(filtered_data_tagged)
plotting_data = filtered_data_tagged[order(HR)]
order = c("CI", "HR")
plotting_data$hr_type <- factor(plotting_data$hr_type, levels = order)
plotting_data$lncRNA <- factor(plotting_data$lncRNA, levels = orderhr)


#Plot summary of HR Confidence Intervals for lncRNA 

##---------Main plot with Hazard Ratios Horizontally spread out-----------------------

g  = ggplot(plotting_data, aes(x=lncRNA, y=HR, col=sig)) + 
  scale_shape_manual(values=c(20, 15)) +
  geom_path(size = 0.2, aes(color=sig)) + 
  geom_point(size=0.65, aes(shape=hr_type), color="black") + 
  theme_classic() + 
  geom_hline(yintercept=0, linetype="dashed", color = "red") + coord_flip() +
  theme(axis.title.y=element_blank()) + theme_minimal()

g = ggpar(g, font.ytickslab = c(3, "plain", "black"), legend = "none") 
main_plot = g


##---------Covariate heatmap for cancer type-------------------------------------------

#color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
#mypal = sample(color, 23)
#saveRDS(mypal, file="palette_23_cancer_types.rds")

mypal = readRDS("palette_23_cancer_types.rds")

cancers_conv = rna[,which(colnames(rna) %in% c("type", "Cancer"))]
cancers_conv = cancers_conv[!duplicated(cancers_conv), ]
plotting_data = merge(plotting_data, cancers_conv, by="Cancer")

cancers = ggplot(plotting_data, aes(lncRNA, 0.2)) +
    geom_tile(aes(fill = Cancer)) + geom_text(aes(label = type), size=1.5) +
    theme_classic() + scale_fill_manual(values=mypal) +
    coord_flip()

cancers = ggpar(cancers, legend = "none") + theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())


##---------Plot outlining the size of risk/non-risk ratio----------------------------

#x-axis lncRNA, y-axis is the ratio (remember it's log2(#pats in risk group/#pats in non-risk group))
plotting_data$ratio = as.numeric(plotting_data$ratio)
#plotting_data$lncRNA = as.character(plotting_data$lncRNA)
plotting_data$ratio = log2(plotting_data$ratio)

ratios = ggplot(plotting_data, aes(x=lncRNA, y=ratio, col=pval)) +
  geom_path(size = 0.2) + 
  geom_point(size=0.8) + 
  scale_color_gradient(low="darkcyan", high="red") +
  theme_classic() + xlab("log2(#risk/#non-risk") +
  geom_hline(yintercept=0, linetype="dashed", color = "azure4") + coord_flip() + theme_minimal() +
  geom_hline(yintercept=log2(0.05), linetype="dashed", color = "azure4")+
  geom_hline(yintercept=log2(5), linetype="dashed", color = "azure4")

ratios = ggpar(ratios, legend = "right") + theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

#pdf("test.pdf", height=15)
#ratios
#dev.off()

##---------Combine cancers + HRs ----------------------------------------------------

pdf("Summary_lncRNAs_hazard_ratios_Risk_patient_ratio_June28.pdf", width=8, height=10)
cancers + main_plot + ratios + plot_layout(ncol = 3, widths = c(1, 12, 4)) 
dev.off()

saveRDS(plotting_data, file="166_candidates_many_levels_of_data_summarized_June27.rds")


##---------Correlation between Ratio and Hazard Ratio---------------------------------

#1. Subset to just HR not CI

plotting_data = as.data.table(plotting_data)
just_hrs = filter(plotting_data, hr_type == "HR")

#2. Correlation between hazard ratio from Cox and size of risk group 
#color point by pvalue?
just_hrs$pcawg_pval = as.numeric(just_hrs$pcawg_pval)
just_hrs$pcawg_pval =  -log10(just_hrs$pcawg_pval)

pdf("correlation_HR_ratio_risk_tononrisk.pdf")

ggscatter(just_hrs, x = "HR", y = "ratio", alpha=0.5,
   fill = "pcawg_match", shape = 21, size = 2,  # Points color, shape and size
   palette=c("snow", "magenta3"), 
   label = "pcawg_match", 
   repel = TRUE, font.label = c(5, "plain"), 
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   #conf.int = TRUE, # Add confidence interval
   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
   ) + 
geom_hline(yintercept=1, linetype="dashed", color = "azure4") +
geom_vline(xintercept=0, linetype="dashed", color = "azure4")
dev.off()


##---------Summarize how many lncRNAs are stable split favourable vs unfavourable------




















