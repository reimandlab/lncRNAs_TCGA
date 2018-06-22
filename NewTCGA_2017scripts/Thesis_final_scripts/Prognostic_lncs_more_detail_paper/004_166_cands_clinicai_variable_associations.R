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

#------FEATURES-----------------------------------------------------

cands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_May4.rds")

#new cands - final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds

#---new cands file 
#final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds

#cands = filter(cands, data == "PCAWG", pval <=0.05)
cands = filter(cands, AnalysisType == "noFDR")
#colnames(cands)[7] = "canc"
cands$Cancer = NULL
all_cands = cands


#--------This script ------------------------------------------------

#include additional clinical variables from more detailed
#files for each cancer type 
#fit multivariate models using those variables for each candidate
#foresplots?
#correlation plots?

#--------------------------------------------------------------------
#Clinical files
#--------------------------------------------------------------------

clin_files = list.files("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/gdc_clinical_data_june2018")

names(clin_files) = c("BRCA", "KIRC", "LGG", "LIHC", "LUAD", "OV", "PAAD")

#write function that adds tag to whole data group 
#and does survival analysis on whole group

z = which(cancers %in% all_cands$cancer)
cancer_data = canc_datas[z] #cancers list and canc_datas list should be the same 

get_canc_data_for_plot = function(dtt){
  #get cancer specific candidates 
  z = which(colnames(dtt) %in% c(as.character(all_cands$gene[all_cands$cancer == dtt$Cancer[1]]), "age_at_initial_pathologic_diagnosis", 
    "OS.time", "OS", "gender", "race", "patient", "clinical_stage", "histological_grade", "treatment_outcome_first_course", 
    "new_tumor_event_type", "Cancer"))
  dtt = dtt[,z]
  return(dtt)
}

filtered_data = llply(cancer_data, get_canc_data_for_plot)

#--------ADD CLINICAL VARIABLES----------------------------------------

add_clin_vars = function(dtt){
  canc = dtt$Cancer[1]
  canc = rna$type[rna$Cancer == canc][1]
  z = which(names(clin_files) %in% canc)
  if(!(length(z)==0)){

    clin = fread(paste("gdc_clinical_data_june2018/", clin_files[z], sep=""), data.table=F)
    print(length(colnames(clin)))
    print(length(unique(clin$bcr_patient_barcode)))
    for(i in 1:ncol(clin)){
      print(table(clin[,i]))
    }

    #check which columns have enoguh contrasts
    #remove columns where # of NAs is greater than 50% of patietnt cohort
    
    check_nas = function(col){
      check = length(which((col == "[Not Applicable]") | (col == "[Not Available]")))
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

    z = which(colnames(clin) %in% "bcr_patient_barcode")
    colnames(clin)[z] = "patient"
    cols = colnames(clin)[which(colnames(clin) %in% colnames(dtt))]

    dtt = merge(dtt, clin, by=cols)
    return(dtt)

  }#end length(z)==0

} #end add_clin_vars 



clin_data_lncs = llply(filtered_data, add_clin_vars)

#remove Nulls
clin_data_lncs = Filter(Negate(is.null), clin_data_lncs)
#save this file can work on at home 
saveRDS(clin_data_lncs, file="clin_data_lncs_new_variables_June21.rds")

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
    new_dat = dtt[,-z]
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
    z = as.numeric(which((all_cands$cancer %in% canc) & (all_cands$gene %in% lnc) & (all_cands$data == "TCGA")))
    hr = all_cands$HR[z]
    new_dat$risk = ""
    if(hr >1){new_dat$risk = "HighExp"}
    if(hr <1){new_dat$risk = "LowExp"}
    
    #each clinical variable
    canc_col_results = as.data.frame(matrix(ncol=6)) ; colnames(canc_col_results)=c("canc", "lnc", "colname", "cor", "pval", "test")
    for(i in 1:ncol(new_dat)){
      print(i)    
      col = colnames(new_dat)[i]
      print(col)
      if(!(col %in% c("patient", "patient_id", "bcr_patient_uuid", "tissue_source_site", 
        "last_contact_days_to", "days_to_initial_pathologic_diagnosis", "tumor_tissue_site", 
        "form_completion_date"))){

        new_dat_plot = new_dat[,c("patient", col, lnc, "lncRNA_tag", "risk")]
        test = as.numeric(new_dat_plot[,2])  
        
        if(str_detect(col, "year")){
          test[1] = 5
        }

        if(!(is.na(test[1]))){
          new_dat_plot[,2] = as.numeric(new_dat_plot[,2])
          colnames(new_dat_plot)[2] = "Clinical"
          new_dat_plot[,3] = log1p(new_dat_plot[,3])
          colnames(new_dat_plot)[3] = "lncRNA_exp"

          #get correlation results and save results into file
          cor = rcorr(new_dat_plot$lncRNA_exp, new_dat_plot$Clinical, "spearman")$r[2]
          pval_cor = rcorr(new_dat_plot$lncRNA_exp, new_dat_plot$Clinical, "spearman")$P[2]
          row = c(canc, lnc, col, cor, pval_cor, "spearman")
          names(row) = colnames(canc_col_results)
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

        if(is.na(test)[1]){
        #boxplot
        colnames(new_dat_plot)[2] = "Clinical"
        new_dat_plot[,3] = log1p(new_dat_plot[,3])
        med = median(new_dat_plot[,3])
        colnames(new_dat_plot)[3] = "lncRNA_exp"
        #palette
        colourCount = length(unique(new_dat_plot$Clinical))
        getPalette = colorRampPalette(brewer.pal(9, "Set1"))

        check = dim(table(new_dat_plot$Clinical))
        if(check >1){

        m1 = lm(new_dat_plot$lncRNA_exp ~1)
        m2 = lm(new_dat_plot$lncRNA_exp ~ new_dat_plot$Clinical)
        anova = anova(m2, m1)
        anova = anova[2,6]

        row = c(canc, lnc, col, "nocor", anova, "Ftest")
        names(row) = colnames(canc_col_results)
        canc_col_results = rbind(canc_col_results, row)

        new_dat_plot$lncRNA_exp = as.numeric(new_dat_plot$lncRNA_exp)
        z = which(new_dat_plot$Clinical %in% c("[Unknown]", "[Not Available]"))
        if(!(length(z)==0)){
          new_dat_plot = new_dat_plot[-z,]
        }

        p <- ggboxplot(new_dat_plot, x = "Clinical", y = "lncRNA_exp",
          color = "Clinical",
          title = paste(canc, lnc, col), 
          add = "jitter", ylab = "lncRNA expression",  ggtheme = theme_bw()) +
          stat_compare_means() + geom_hline(yintercept=med, linetype="dashed", color = "red") + 
          stat_n_text() + scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Accent"))(colourCount))

        p = ggpar(p,
          font.tickslab = c(10,"plain", "black"),
          xtickslab.rt = 45, legend="none")
        print(p)
          
          } #check >1 
        }
      }
    } #for i in 1:ncol(new_dat)
    canc_col_results = canc_col_results[-1,]
    return(canc_col_results)

  } #end get_cor

  pdf(paste(canc, "clinical_plots.pdf", sep="_"), width=10)
  all_canc_lncs_results = llply(lncs, get_cor)
  all_canc_lncs_results = do.call(rbind.data.frame, all_canc_lncs_results)
  dev.off()
  return(all_canc_lncs_results)

} #end get_clin_lnc_cors


clin_data_lncs_cors = llply(clin_data_lncs, get_clin_lnc_cors)

#--------FDR & Summarize Results-------------------------------------

fdr_sum = function(dtt){

  #first add fdr for p-values
  dtt$fdr = p.adjust(dtt$pval, method="fdr")
  dtt = as.data.table(dtt)
  dtt = filter(dtt, fdr < 0.05)

  #remove OS.time, OS, lncRNA tag... 
  z = which(dtt$colname %in% c("OS", "OS.time", "lncRNA_tag", "vital_status"))
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

clean_up = llply(clin_data_lncs_cors, fdr_sum)
clean_up = ldply(clean_up, data.table)
clean_up = as.data.table(clean_up)
clean_up = clean_up[order(fdr)]

saveRDS(clean_up, file="correlation_results_clinical_lncRNA_exp_June22.rds")



















