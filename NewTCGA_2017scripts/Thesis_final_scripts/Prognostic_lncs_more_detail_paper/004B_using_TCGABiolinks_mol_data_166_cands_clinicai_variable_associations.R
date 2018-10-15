library(survAUC)
#source("source_code_Cox_MonteCarlo_CV_Mar13.R")
require(caTools)
#check if this person is in my analysis: TCGA-61-2095
library(glmnet)
library(survcomp)
library(caret)
library(stringr)

source("universal_LASSO_survival_script.R")

library(ggpubr)
library(ggrepel)
library(viridis)
library(patchwork)
library(caret)  
library(Rtsne)
library(EnvStats)

source("check_lnc_exp_cancers.R")

#------FEATURES-----------------------------------------------------

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = subset(allCands, data == "TCGA") #173 unique lncRNA-cancer combos, #166 unique lncRNAs 
allCands$combo = unique(paste(allCands$gene, allCands$cancer, sep="_"))
cands_dups = unique(allCands$gene[which(duplicated(allCands$gene))])

val_cands = read.csv("175_lncRNA_cancers_combos_23_cancer_types_july5.csv")
val_cands = as.data.table(val_cands)
val_cands = subset(val_cands, data == "PCAWG") #173 unique lncRNA-cancer combos, #166 unique lncRNAs 
val_cands$combo = unique(paste(val_cands$gene, val_cands$cancer, sep="_"))
val_cands = subset(val_cands, top_pcawg_val == "YES") #173 unique lncRNA-cancer combos, #166 unique lncRNAs 

library(TCGAbiolinks)

#--------This script ------------------------------------------------

#include additional clinical variables from more detailed
#files for each cancer type 
#fit multivariate models using those variables for each candidate
#foresplots?
#correlation plots?

#--------------------------------------------------------------------
#Clinical files - use TCGAbiolinks
#--------------------------------------------------------------------

#write function that adds tag to whole data group 
#and does survival analysis on whole group

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

#subtypes available from biolinks
subtypes_data = toupper(c("acc", "brca", "coad", "gbm", "hnsc", "kich", "kirp", 
  "kirc", "lgg", "luad", "lusc", "prad", "pancan", "read", "skcm", "stad", "thca", "ucec"))

#--------ADD CLINICAL VARIABLES----------------------------------------

add_clin_vars = function(dtt){
  canc = dtt$Cancer[1]
  canc = rna$type[rna$Cancer == canc][1]

  #Check if TCGA has 
  z = which(subtypes_data %in% canc)
  if(!(length(z)==0)){
  clin_subtypes <- TCGAquery_subtype(tumor = canc)

    clin = clin_subtypes
    z = which(clin$patient %in% dtt$patient)
    clin = clin[z,]
    print(length(colnames(clin)))
    print(length(unique(clin$patient)))
    for(i in 1:ncol(clin)){
      print(table(clin[,i]))
    }

    #check which columns have enoguh contrasts
    #remove columns where # of NAs is greater than 50% of patietnt cohort
    
    check_nas = function(col){
      check = length(which((col == "[Not Applicable]") | (col == "[Not Available]") | (col == "Unknown")))
        if(check < (dim(clin)[1] *0.5)){
          return("keep")
        }
      }
    
    keep_cols1 = unlist(apply(clin, 2, check_nas))
    clin = clin[,which(colnames(clin) %in% names(keep_cols1))]

    check_contrasts = function(col){
      check = dim(table(col))
        if(check >1){
          return("keep")
        }
      }
    
    keep_cols2 = unlist(apply(clin, 2, check_contrasts))
    clin = clin[,which(colnames(clin) %in% names(keep_cols2))]

    cols = colnames(clin)[which(colnames(clin) %in% colnames(dtt))]

    dtt = merge(dtt, clin, by=cols)
    return(dtt)

    } #end add_clin_vars 

  #if not in molecular profiles subset of biolinks
  #just look at whatever clinical variables are available
  #if(length(z)==0){
  #  clinical <- GDCquery_clinic(project = paste("TCGA-", canc, sep=""), type = "clinical")
  #}

}

clin_data_lncs = llply(filtered_data, add_clin_vars)

#remove Nulls
clin_data_lncs = Filter(Negate(is.null), clin_data_lncs)
#save this file can work on at home 
#saved file --- below
#saveRDS(clin_data_lncs, file="clin_data_lncs_new_variables_July19_tcgabiolinks_data.rds")

#--------LOOK AT ASSOCIATIONS BETWEEN EXPRESSION-------------------------------

#For each clinical variable -> xaxis is the clinical variable
#y-axis it each lncRNAs expression 
#x-axis is continous if variable is continous such as age...

get_clin_lnc_cors = function(dtt){
  canc = dtt$Cancer[1]
  
  #get lncs
  z = which(str_detect(colnames(dtt), "ENSG")) 
  lncs = colnames(dtt)[z]
  
  #look at individual lncRNAs 
  get_cor = function(lnc){
    z = which((str_detect(colnames(dtt), "ENSG") & !(colnames(dtt) %in% lnc)))
    new_dat = dtt
    if(length(z) > 0){
    new_dat = dtt[,-z]}
    #add 0/1 labels 
    new_dat$lncRNA_tag = ""
    med = median(new_dat[,which(colnames(new_dat) %in% lnc)])
    k = which(colnames(new_dat) %in% lnc)
    if(med ==0){
        #if median = 0 then anyone greater than zero is 1 
        l1 = which(new_dat[,k] > 0)
        l2 = which(new_dat[,k] ==0)
        new_dat$lncRNA_tag[l1] = 1
        new_dat$lncRNA_tag[l2] = 0
        }

        if(!(med ==0)){
        l1 = which(new_dat[,k] >= med)
        l2 = which(new_dat[,k] < med)
        new_dat$lncRNA_tag[l1] = 1
         new_dat$lncRNA_tag[l2] = 0
        }
    #get risk type 
    z = as.numeric(which((allCands$cancer %in% canc) & (allCands$gene %in% lnc) & (allCands$data == "TCGA")))
    hr = as.numeric(allCands$HR[z])
    new_dat$risk = ""
    if(hr >1){new_dat$risk = "HighExp"}
    if(hr <1){new_dat$risk = "LowExp"}
    
    #each clinical variable
    canc_col_results = as.data.frame(matrix(ncol=16)) ; colnames(canc_col_results)=c("canc", "lnc", "colname", "cor", "pval", "test", "chisq", "kw_pval",
    "clin_pval", "anova_both_vs_lnc", "lnc_concordance", "clin_concordance", "lnc_HR", "clin_HR", "concordance_combo_model", "clin_vs_combo_anova")

    for(i in 1:(ncol(new_dat)-2)){
      print(i)    
      col = colnames(new_dat)[i]
      if((!(col == lnc)) & (!(str_detect(col, "RPPA")))){

      if(!(is.numeric(new_dat[,i]))){
      new_dat[,i] = as.character(new_dat[,i])}
      
      print(col)
      if(!(col %in% c("patient", "patient_id", "bcr_patient_uuid", "tissue_source_site", 
        "last_contact_days_to", "days_to_initial_pathologic_diagnosis", "tumor_tissue_site", 
        "form_completion_date", "OS.time", "OS", "days_to_death", "Signet.Ring", "MACIS"))){

        new_dat_plot = new_dat[,c("patient", col, lnc, "lncRNA_tag", "risk")]
        test = as.numeric(new_dat_plot[,2])  
        
        if(str_detect(col, "year")){
          test[1] = 5
        }

        if(!(length(which(is.na(test))) == length(test))){
        
          z = test[which(!(is.na(test)))]
          if((length(z) > 10) & !(length(z) == length(which(test==0)))){

          #if(!(is.na(test[1]))){
          new_dat_plot[,2] = as.numeric(new_dat_plot[,2])
          colnames(new_dat_plot)[2] = "Clinical"
          new_dat_plot[,3] = log1p(new_dat_plot[,3])
          colnames(new_dat_plot)[3] = "lncRNA_exp"

          #get correlation results and save results into file
          cor = rcorr(new_dat_plot$lncRNA_exp, new_dat_plot$Clinical, "spearman")$r[2]
          pval_cor = rcorr(new_dat_plot$lncRNA_exp, new_dat_plot$Clinical, "spearman")$P[2]
          chisq_pval = "nochisq"
          kw_pval = "nokw"
          
          #how good of a predictor of survial is the clinical variable itself? 
          surv_dat = rna[,which(colnames(rna) %in% c("patient", "OS", "OS.time"))]
          new_dat_plot = merge(new_dat_plot, surv_dat, by = c("patient"))
          new_dat_plot$OS = as.numeric(new_dat_plot$OS)
          new_dat_plot$OS.time = as.numeric(new_dat_plot$OS.time)

          cox_lnc = coxph(Surv(OS.time, OS) ~ lncRNA_tag, data = new_dat_plot)
          cox_clin = coxph(Surv(OS.time, OS) ~ Clinical, data = new_dat_plot)
          both = coxph(Surv(OS.time, OS) ~ lncRNA_tag + Clinical, data = new_dat_plot)
          
          clin_concordance = glance(cox_clin)$concordance
          lnc_concordance = glance(cox_lnc)$concordance
          combo_concordance = glance(both)$concordance

          hr_clin = summary(cox_clin)$coefficients[2]

          clin_pval = glance(cox_clin)[6]
          z = which(is.na(new_dat_plot[,2]))
          if(!(length(z)==0)){
            anov_pval = "cant_calc"
          }

          if(length(z)==0){
          anov_pval = anova(cox_lnc, both)[2,4]
          clin_vs_combo_anova = anova(cox_clin, both)[2,4]
          }

          #add lnc, clinical Hazard Ratios and conordance of combined model 

          row = c(canc, lnc, col, cor, pval_cor, "Ftest", chisq_pval, kw_pval, 
          clin_pval, anov_pval, lnc_concordance, clin_concordance, hr, hr_clin, combo_concordance, clin_vs_combo_anova)
          names(row) = colnames(canc_col_results)

          print(ggforest(both, main = paste(lnc, col, canc), data=new_dat_plot))

          canc_col_results = rbind(canc_col_results, row)
          #scatter plot 
          sp <- ggscatter(new_dat_plot, x = "Clinical", y = "lncRNA_exp",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE # Add confidence interval
          )
          # Add correlation coefficient
          sp = sp + stat_cor(method = "spearman") + theme_bw() + ggtitle(paste(canc, lnc, col))
          print(sp)

        }
        }

        #if(is.na(test)[1]){
        if(length(which(is.na(test))) == length(test)){
        #boxplot
        
        #remove catgeories with less than 5 patients 
        t = as.data.table(table(new_dat_plot[,2]))
        t = filter(t, N < 10)
        rm = unique(t$V1)
        if(!(length(rm) ==0)){
          new_dat_plot = new_dat_plot[-(which(new_dat_plot[,2] %in% rm)),]
        }

        if(!(length(unique(new_dat_plot$Clinical)) > 10)){

        #palette
        colourCount = length(unique(new_dat_plot$Clinical))
        getPalette = colorRampPalette(brewer.pal(9, "Set1"))

        check = dim(table(new_dat_plot[,which(colnames(new_dat_plot) %in% col)]))
        if(check >1){

        #remove any NAs 
        #remove NAs

        new_dat_plot[,3] = log1p(new_dat_plot[,3])
        med = median(new_dat_plot[,3])
        colnames(new_dat_plot)[3] = "lncRNA_exp"

        z1 = which(is.na(new_dat_plot[,which(colnames(new_dat_plot) %in% col)]))  
        z2 = which(new_dat_plot[,which(colnames(new_dat_plot) %in% col)] %in% c("#N/A", "Unknown", "N/A", "NA", "Not Available", "Not performed", "Performed but Not Available"))  
        z3 = which(new_dat_plot[,which(colnames(new_dat_plot) %in% col)] %in% c("[Unknown]", "[Not Available]", "[Not Evaluated]", "[Discrepancy]"))  

        z = unique(c(z1, z2,z3))
        
        if(!(length(z)==0)){
        new_dat_plot = new_dat_plot[-z,]}

        unq = length(unique(new_dat_plot[,2]))

        if(unq > 1){

        if(dim(new_dat_plot)[1] > 10){

        colnames(new_dat_plot)[2] = "Clinical"
        
        m1 = lm(new_dat_plot$lncRNA_exp ~1)
        m2 = lm(new_dat_plot$lncRNA_exp ~ 1 + new_dat_plot$Clinical)
        anova = anova(m1, m2)
        anova = anova[2,6]

        new_dat_plot$Clinical = as.factor(new_dat_plot$Clinical)
        kw_pval = as.numeric(tidy(kruskal.test(lncRNA_exp ~ Clinical, data = new_dat_plot))[2])

        #do Chisq test of independence 
        tb = table(new_dat_plot$lncRNA_tag, new_dat_plot$Clinical)
        chisq_pval = as.numeric(tidy(chisq.test(tb))[2])

        #how good of a predictor of survial is the clinical variable itself? 
        surv_dat = rna[,which(colnames(rna) %in% c("patient", "OS", "OS.time"))]
        new_dat_plot = merge(new_dat_plot, surv_dat, by = c("patient"))
        new_dat_plot$OS = as.numeric(new_dat_plot$OS)
        new_dat_plot$OS.time = as.numeric(new_dat_plot$OS.time)

        num_high = length(which(new_dat$lncRNA_tag ==1))
        num_low = length(which(new_dat$lncRNA_tag ==0))

        lncheck = ((num_low >=10) & (num_high >=10))

        #make sure medians of two groups aren't the same
        #ie both are 0 then effect isn't really significant 
        med_check = as.data.table(new_dat_plot %>% group_by(Clinical) %>% summarise_each(funs(median),lncRNA_exp))
        med_check = !(med_check$lncRNA_exp[1] == med_check$lncRNA_exp[2])

        if((dim(table(new_dat_plot$lncRNA_tag)) > 1) & lncheck & med_check){

        cox_lnc = coxph(Surv(OS.time, OS) ~ lncRNA_tag, data = new_dat_plot)
        cox_clin = coxph(Surv(OS.time, OS) ~ Clinical, data = new_dat_plot)
        both = coxph(Surv(OS.time, OS) ~ lncRNA_tag + Clinical, data = new_dat_plot)

        clin_pval = glance(cox_clin)[6]
        
        clin_concordance = glance(cox_clin)$concordance
        lnc_concordance = glance(cox_lnc)$concordance
        combo_concordance = glance(both)$concordance
        hr_clin = summary(cox_clin)$coefficients[2]
        anov_pval = anova(cox_lnc, both)[2,4]
        clin_vs_combo_anova = anova(cox_clin, both)[2,4]
        print(ggforest(both, main = paste(lnc, col, canc), data=new_dat_plot))

        }

        if((dim(table(new_dat_plot$lncRNA_tag))) <= 1 & (!(check))){
        clin_pval = "cant_calc"
        anov_pval = "cant_calc"
        clin_concordance = "cant_calc"
        lnc_concordance = "cant_calc"
        combo_concordance = "cant_calc"
        hr_clin = "cant_calc"
        }

        row = c(canc, lnc, col, "nocor", anova, "Ftest", chisq_pval, kw_pval, 
          clin_pval, anov_pval, lnc_concordance, clin_concordance, hr, hr_clin, combo_concordance, clin_vs_combo_anova)
        names(row) = colnames(canc_col_results)
        canc_col_results = rbind(canc_col_results, row)

        meds = as.data.table(aggregate(lncRNA_exp ~ Clinical, new_dat_plot, median))
        meds = meds[order(lncRNA_exp)]
        new_dat_plot$Clinical = factor(new_dat_plot$Clinical, levels = unique(meds$Clinical))

        p <- ggboxplot(new_dat_plot, x = "Clinical", y = "lncRNA_exp",
          color = "Clinical",
          title = paste(canc, lnc, col), 
          add = "jitter", ylab = "lncRNA expression",  ggtheme = theme_bw()) +
          stat_compare_means() + geom_hline(yintercept=med, linetype="dashed", color = "red") + 
          stat_n_text(size=5)

        p = ggpar(p,
          font.xtickslab = c(14,"plain", "black"),font.tickslab=c(14,"plain", "black"), 
          xtickslab.rt = 55, legend="none")
        print(p)
          
          }
          } #check >1 
        }
      }
    }
    }
    }
    } #for i in 1:ncol(new_dat)
    canc_col_results = canc_col_results[-1,]
    return(canc_col_results)

  } #end get_cor

  pdf(paste(canc, "clinical_plots.pdf", sep="_"))
  all_canc_lncs_results = llply(lncs, get_cor)
  all_canc_lncs_results = do.call(rbind.data.frame, all_canc_lncs_results)
  dev.off()
  return(all_canc_lncs_results)

} #end get_clin_lnc_cors

#clin_data_lncs_cors = llply(clin_data_lncs, get_clin_lnc_cors)
d1 = get_clin_lnc_cors(clin_data_lncs[[1]])
d2 = get_clin_lnc_cors(clin_data_lncs[[2]])
d3 = get_clin_lnc_cors(clin_data_lncs[[3]])
d4 = get_clin_lnc_cors(clin_data_lncs[[4]])
d5 = get_clin_lnc_cors(clin_data_lncs[[5]])
d6 = get_clin_lnc_cors(clin_data_lncs[[6]])
d7 = get_clin_lnc_cors(clin_data_lncs[[7]])
d8 = get_clin_lnc_cors(clin_data_lncs[[8]])
d9 = get_clin_lnc_cors(clin_data_lncs[[9]])
d10 = get_clin_lnc_cors(clin_data_lncs[[10]]) #<------ this one didn't work
d11 = get_clin_lnc_cors(clin_data_lncs[[11]])
d12 = get_clin_lnc_cors(clin_data_lncs[[12]])
d13 = get_clin_lnc_cors(clin_data_lncs[[13]])
d14 = get_clin_lnc_cors(clin_data_lncs[[14]])

all_clin = list(d1, d2,d3,d4,d5,d6,d7,d8,d9,d11,d12,d13,d14)
saveRDS(all_clin, file="13_data_sets_biolinks_results.rds")

#--------FDR & Summarize Results-------------------------------------
all_clin = readRDS("13_data_sets_biolinks_results.rds")

fdr_sum = function(dtt){

  #first add fdr for p-values
  dtt$fdr = p.adjust(dtt$pval, method="fdr")
  dtt = as.data.table(dtt)
  #dtt = filter(dtt, fdr < 0.05)

  #remove OS.time, OS, lncRNA tag... 
  z = which(dtt$colname %in% c("OS", "OS.time", "lncRNA_tag", "vital_status", 
    "Tissue.source.site", "Whole.genome", "SNP6", "HM450", "HM27", "Vital.status..1.dead.", 
    "COC", "Status", "Survival..months.", "OS.event", "Days.to.date.of.Death", "OS.event", "BCR", 
    "Tumor", "X2009stagegroup", "time_of_follow.up", "CDE_ID.3045435", "batch", "Survival", "Exome.data", "CDE_ID.3104937","OS.Time",
    "Country", "os_days", "CDE_ID.2006657", "icd_o_3_site", "WGS"))
  if(!(length(z)==0)){
    dtt = dtt[-z,]
  }

 #remove if it's just the same lncRNA correalted with itself 
 z = which(dtt$lnc == dtt$colname)
 if(!(length(z)==0)){
 dtt = dtt[-z,]}
 dtt = as.data.table(dtt)
 dtt = dtt[order(fdr)]
 print(unique(dtt$colname))
 return(dtt)

}#end fdr_sum

clean_up = llply(all_clin, fdr_sum)
clean_up = ldply(clean_up, data.table)
clean_up = as.data.table(clean_up)
clean_up = clean_up[order(fdr)]

#look at just chisq tests
clean_up$chisq = as.numeric(clean_up$chisq)
z = which(clean_up$chisq <= 0.05)
clean_up$sig_chisq[z] = "*"

canc_conv = rna[,c("type", "Cancer")]
canc_conv = unique(canc_conv)
colnames(canc_conv)[2] = "canc"
clean_up = merge(clean_up, canc_conv, by="canc")
clean_up = clean_up[order(fdr)]
clean_up$type = factor(clean_up$type, levels=unique(clean_up$type))

clean_up$colname[which(str_detect(clean_up$colname, "age_at"))] = "Age"
clean_up$colname[clean_up$colname == "age"] = "Age"

fdr = as.data.table(filter(clean_up, fdr < 0.05))
fdr$cor = as.numeric(fdr$cor)
fdr$kw_pval = as.numeric(fdr$kw_pval)
fdr = filter(fdr, (is.na(cor) | (abs(cor) >= 0.3)), kw_pval < 0.05)

write.csv(fdr, file="cleaned_clinical_variables_associations_data_sept19_precleanup.csv", quote=F, row.names=F)

pdf("summary_biolinks_subtypes_lncRNA_exp.pdf", height=10)
#make geom_tile plot
ggplot(clean_up, aes(type, colname)) +
  geom_tile(aes(fill = sig_chisq), colour = "grey50") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size=5))
dev.off()

saveRDS(clean_up, file="correlation_results_clinical_lncRNA_exp_July19_using_biolinks.rds")
write.table(clean_up, file="correlation_results_clinical_lncRNA_exp_July19_using_biolnks.txt", row.names=F, quote=F)

#-------PLOT summary results-------------------------------------------

#clin_results = readRDS("correlation_results_clinical_lncRNA_exp_July19_using_biolinks.rds")
clin_results = read.csv("cleaned_clinical_variables_associations_data_sept19_precleanup.csv")
clin_results$combo = paste(clin_results$lnc, clin_results$type, sep="_")
clean_up = clin_results

clean_up = as.data.table(filter(clean_up, fdr < 0.05))
clean_up$cor = as.numeric(clean_up$cor)
clean_up$kw_pval = as.numeric(clean_up$kw_pval)

clean_up = filter(clean_up, (is.na(cor) | (abs(cor) >= 0.3)), kw_pval < 0.05)
clean_up$better = ""
z = which((clean_up$concordance_combo_model > clean_up$lnc_concordance) & (clean_up$anova_both_vs_lnc < 0.05))
clean_up$better[z] = "V"

get_name = function(ensg){
    z = which(fantom$CAT_geneID == ensg)
    return(fantom$CAT_geneName[z][1])
}

clean_up$name = unlist(llply(clean_up$lnc, get_name))

#split by cancer type
cancer_clins = split(clean_up, clean_up$type)

get_nodes_edges = function(canc){

  canc = as.data.table(filter(canc, pval < 0.05))

  #make node file
  lncs = as.data.frame(unique(canc$name))
  lncs$node_type = "lncRNA"
  colnames(lncs)[1] = "node"

  clinical = as.data.frame(unique(canc$colname))
  clinical$node_type = "clinical"
  colnames(clinical)[1] = "node"
  nodes = rbind(clinical, lncs)
  nodes$canc = canc$type[1]
  write.table(nodes, file=paste(canc$type[1], "clinical_biolinks_analysis_nodes_table_sept14.txt", sep="_"), quote=F, row.names=F, sep="\t")
  
  #make edge file 
  #node1  node2   edge_type
  edges = canc[,c("name", "colname", "pval", "lnc_concordance", "clin_concordance")]
  edges = unique(edges)
  edges$edge_type = "lnc_clin"
  write.table(edges, file=paste(canc$type[1], "clinical_biolinks_analysis_edges_table_sept14.txt", sep="_"), quote=F, row.names=F, sep="\t")

  g = ggplot(edges, aes(lnc_concordance, clin_concordance, label=name)) + ggtitle(canc$type[1], "lncRNA vs Clinical variables")+
  geom_point()+
    xlab("lncRNA Concordance") + ylab("Clinical Concordance") + theme_bw() +
      xlim(0.45,0.95) + ylim(0.45,0.95) + geom_abline(intercept=0)
  print(g)
  print(paste("done", canc$type[1]))
}

pdf("lncRNA_vs_Clinical_variables_CIs_sep14.pdf")
llply(cancer_clins, get_nodes_edges)
dev.off()

clean_up$combo = paste(clean_up$name, clean_up$type, sep = " ")
clean_up$sig_tag = ""
clean_up$clin_pval = p.adjust(clean_up$clin_pval, method="fdr")
clean_up$sig_tag[clean_up$clin_pval < 0.05] = "V"
clean_up$sig_tag[clean_up$clin_pval > 0.05] = ""

write.csv(clean_up, file="cleaned_clinical_variables_associations_data_sept28_post_cleanup.csv", quote=F, row.names=F)

#post manual cleanup of variables 
clin_results = read.csv("cleaned_clinical_variables_associations_data_sept28_post_cleanup.csv")

clin_results = as.data.table(clin_results)

#which combos are better once lncRNA is used
#look at only those where clinical variable also associated with survival
clinsig = as.data.table(filter(clin_results, sig_tag == "V"))
clinsig_cause_lnc = as.data.table(filter(clinsig, clin_vs_combo_anova < 0.05)) #121/125
clinsig_notcause_lnc = as.data.table(filter(clinsig, clin_vs_combo_anova > 0.05)) #4/125

#36 unique lncRNA-cancer combos 

#get order 
t = as.data.table(table(clin_results$colname))
t = as.data.table(filter(t, N > 0))
t = t[order(-N)]

clin_results$colname = factor(clin_results$colname, levels = t$V1)

#get order 
t = as.data.table(table(clin_results$combo))
t = as.data.table(filter(t, N > 0))
t = t[order(-N)]

clin_results$combo = factor(clin_results$combo, levels = t$V1)
clin_results$concordance_combo_model = as.numeric(clin_results$concordance_combo_model)
clin_results$anova_both_vs_lnc = as.numeric(clin_results$anova_both_vs_lnc)
z = which((clin_results$concordance_combo_model > clin_results$lnc_concordance) & (clin_results$anova_both_vs_lnc < 0.05))
clin_results$better[z] = "V"

pdf("summary_biolinks_subtypes_lncRNA_exp_Sept28.pdf", height=6, width=12)
#make geom_tile plot
ggplot(clin_results, aes(combo, colname)) +
  geom_tile(aes(fill = concordance_combo_model), colour = "grey50") +
  theme_bw() + geom_text(aes(label = better), size=2.5) + 
  theme(legend.title=element_blank(), legend.position="bottom", axis.title.x=element_blank(), 
    axis.text.x = element_text(angle = 45, hjust = 1, size=7)) +
    scale_fill_gradient(low = "tan1", high = "darkred")
    #+ scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

dev.off()



write.csv(clin_results, file="cleaned_clinical_variables_associations_data_sept28_post_cleanup_final_figure_data.csv", quote=F, row.names=F)


#barplot comparing concordance 
#combo 
#y-axis concordance 
#fill bars by lncRNA only/clinical only/combined 
#clinsig = as.data.table(filter(clin_results, sig_tag == "V"))
clinsig = clin_results
clinsig_cause_lnc = as.data.table(filter(clinsig, clin_vs_combo_anova < 0.05)) #121/125
clinsig_cause_lnc = clinsig
clinsig_cause_lnc$combo_clin = paste(clinsig_cause_lnc$combo, clinsig_cause_lnc$colname)
clinsig_cause_lnc$anova_both_vs_lnc = as.numeric(clinsig_cause_lnc$anova_both_vs_lnc)
clinsig_cause_lnc$better_than_clin = ""
z = which((clinsig_cause_lnc$concordance_combo_model > clinsig_cause_lnc$clin_concordance) & (clinsig_cause_lnc$clin_vs_combo_anova < 0.05))
clinsig_cause_lnc$better_than_clin[z] = "V"

z = which(clinsig_cause_lnc$colname %in% c("Spread to Lymph nodes", "treatment_outcome_first_course"))
clinsig_cause_lnc = clinsig_cause_lnc[-z,]
#z = which(clinsig_cause_lnc$lnc_concordance > 0.5)
#clinsig_cause_lnc= clinsig_cause_lnc[z,]

clinsig_cause_lnc = clinsig_cause_lnc[order(-concordance_combo_model)]

lnc_conc = clinsig_cause_lnc[,c("combo", "combo_clin", "lnc_concordance")]
colnames(lnc_conc)[ncol(lnc_conc)] = "Concordance"
lnc_conc$type = "lncRNA"

clin_conc = clinsig_cause_lnc[,c("combo", "combo_clin", "clin_concordance")]
colnames(clin_conc)[ncol(clin_conc)] = "Concordance"
clin_conc$type = "clin"

combo_conc = clinsig_cause_lnc[,c("combo", "combo_clin", "concordance_combo_model")]
colnames(combo_conc)[ncol(combo_conc)] = "Concordance"
combo_conc$type = "combo"

all = rbind(lnc_conc, clin_conc, combo_conc)

all$combo_clin = factor(all$combo_clin, levels=unique(clinsig_cause_lnc$combo_clin))
all$type = factor(all$type, levels=c("clin", "lncRNA", "combo"))

pdf("summary_clinical_concordances_oct11.pdf", width=8, height=10)
g <- ggplot(all, aes(combo_clin, Concordance))
# Number of cars in each class:
g + geom_col(position = "dodge", aes(fill = type)) + scale_fill_manual(values=mypal[c(1,4,3)], name="Predictor type",
                       breaks=c("clin", "lncRNA", "combo"),
                       labels=c("Clinical", "lncRNA", "Combined"))+
theme_bw()+
theme(axis.text.x = element_text(size=6, angle=45, hjust=1),
          axis.text.y = element_text(size=5), legend.position="top") + xlab("lncRNA-Cancer-Clinical") + ylab("Concordance")+
 coord_flip()
dev.off()

#try scatter plot


pdf("summary_clinical_concordances_scatterplot_oct11.pdf",width=6, height=6)
g = ggplot(clinsig_cause_lnc, aes(lnc_concordance, concordance_combo_model, label=combo_clin)) +
 geom_point(aes(color=type, size=-log10(anova_both_vs_lnc)))+
 #scale_size(range = c(0, 3))+
    scale_colour_manual(values = mypal) + 
    xlab("lncRNA Concordance") + ylab("lncRNA & Clinical Combined Concordance") + theme_bw() +
    theme(legend.position = "top", axis.text = element_text(size=12), 
      legend.text=element_text(size=10), legend.title=element_text(size=10)) +
     xlim(0.5,1) + ylim(0.5,1) + geom_abline(intercept=0) + 
     geom_text_repel(data = subset(clinsig_cause_lnc, combo_clin == "MAPT-AS1 BRCA PAM50.mRNA"), size=2, nudge_y = 0.1,
      direction = "x",segment.color = "grey50",
      segment.size = 0.05)
     
     #geom_text_repel(data = subset(prog_pcgs, lncConcordance > 0.7 | pcgConcordance > 0.65), size=2, nudge_y = 0.1,
      #direction = "x",segment.color = "grey50",
      #segment.size = 0.05)

g
dev.off()

pdf("summary_clinical_concordances_vs_lnc_scatterplot_oct11.pdf", width=6, height=6)
g = ggplot(clinsig_cause_lnc, aes(clin_concordance, concordance_combo_model, label=combo_clin)) +
 geom_point(aes(color=type))+
 #scale_size(range = c(0, 3))+
    scale_colour_manual(values = sample(mypal, length(unique(clinsig_cause_lnc$type)))) + 
    xlab("Clinical Concordance") + ylab("lncRNA & Clinical Combined Concordance") + theme_bw() +
    theme(legend.position = "top", axis.text = element_text(size=12), 
      legend.text=element_text(size=10), legend.title=element_text(size=10)) +
     xlim(0.5,0.9) + ylim(0.5,0.9) + geom_abline(intercept=0) + 
     
     geom_text_repel(data = subset(clinsig_cause_lnc, combo_clin == "HOXB-AS2 LGG TERT.promoter.status"), size=2, nudge_y = 0.1,
      direction = "x",segment.color = "grey50",
      segment.size = 1)
g
dev.off()

clinsig_cause_lnc$imp_concord = (clinsig_cause_lnc$concordance_combo_model - clinsig_cause_lnc$clin_concordance)
clinsig_cause_lnc = clinsig_cause_lnc[order(-better_than_clin,-imp_concord )]


#make KM plot for TERT promoter status
lgg = clin_data_lncs[[13]]
lgg = lgg[,which(colnames(lgg) %in% c("patient", "ENSG00000239552", "OS", "OS.time", "TERT.promoter.status"))]
med = median(lgg$ENSG00000239552)
#median = 0
lgg$med = ""
lgg$med[which(lgg$ENSG00000239552 > 0)] = "High"
lgg$med[which(lgg$ENSG00000239552 == 0)] = "Low"

z = which(is.na(lgg$TERT.promoter.status))
lgg = lgg[-z,]

pdf("lgg_tert_vs_HOXBAS2_km_plots.pdf", width=9)

lgg$OS = as.numeric(lgg$OS)
  lgg$OS.time = as.numeric(lgg$OS.time)
  lgg$OS.time = lgg$OS.time/365
  lgg$med = factor(lgg$med, levels = c("Low","High"))
  fit <- survfit(Surv(OS.time, OS) ~ med, data = lgg)
          s <- ggsurvplot(
          fit, 
          xlab = "Time (Years)", 
          surv.median.line = "hv",
          font.main = c(16, "bold", "black"),
          font.x = c(14, "plain", "black"),
          font.y = c(14, "plain", "black"),
          font.tickslab = c(14, "plain", "black"),
          font.legend = 12,
          risk.table.fontsize = 5, 
          legend.labs = c("Low Expression", "High Expression"),             # survfit object with calculated statistics.
          data = lgg,      # data used to fit survival curves. 
          risk.table = TRUE,       # show risk table.
          legend = "right", 
          pval = TRUE,             # show p-value of log-rank test.
          conf.int = FALSE,        # show confidence intervals for 
                            # point estimaes of survival curves.
          xlim = c(0,10),        # present narrower X axis, but not affect
                            # survival estimates.
          break.time.by = 1,     # break X axis in time intervals by 500.
          #palette = colorRampPalette(mypal)(14), 
          palette = mypal[c(4,1)],
          ggtheme = theme_bw(), # customize plot and risk table with a theme.
          risk.table.y.text.col = T, # colour risk table text annotations.
          risk.table.y.text = FALSE # show bars instead of names in text annotations
                            # in legend of risk table
          )
          print(s)

  fit <- survfit(Surv(OS.time, OS) ~ TERT.promoter.status, data = lgg)
          s <- ggsurvplot(
          fit, 
          xlab = "Time (Years)", 
          surv.median.line = "hv",
          font.main = c(16, "bold", "black"),
          font.x = c(14, "plain", "black"),
          font.y = c(14, "plain", "black"),
          font.tickslab = c(14, "plain", "black"),
          font.legend = 8,
          risk.table.fontsize = 5, 
          #legend.labs = c("Low Expression", "High Expression"),             # survfit object with calculated statistics.
          data = lgg,      # data used to fit survival curves. 
          risk.table = TRUE,       # show risk table.
          legend = "right", 
          pval = TRUE,             # show p-value of log-rank test.
          conf.int = FALSE,        # show confidence intervals for 
                            # point estimaes of survival curves.
          xlim = c(0,10),        # present narrower X axis, but not affect
                            # survival estimates.
          break.time.by = 1,     # break X axis in time intervals by 500.
          #palette = colorRampPalette(mypal)(14), 
          palette = mypal[c(4,1)],
          ggtheme = theme_bw(), # customize plot and risk table with a theme.
          risk.table.y.text.col = T, # colour risk table text annotations.
          risk.table.y.text = FALSE # show bars instead of names in text annotations
                            # in legend of risk table
          )
          print(s)


dev.off()









