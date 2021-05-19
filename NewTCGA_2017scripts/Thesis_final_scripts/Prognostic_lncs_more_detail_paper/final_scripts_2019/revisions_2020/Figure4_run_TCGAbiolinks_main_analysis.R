source("/u/kisaev/lncRNAs_TCGA/NewTCGA_2017scripts/Thesis_final_scripts/Prognostic_lncs_more_detail_paper/final_scripts_2019/revisions_2020/load_data.R")

#------FEATURES-----------------------------------------------------
setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNAs_2019_manuscript")

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = subset(allCands, data == "TCGA") #173 unique lncRNA-cancer combos, #166 unique lncRNAs
allCands$combo = unique(paste(allCands$gene, allCands$cancer, sep="_"))

#library(TCGAbiolinks)
get_median_risk_group = function(PI_centered, PI_lower_thres, epsilon = 0.0001) {

  if (!any(PI_centered > PI_lower_thres)) {
    PI_lower_thres = PI_lower_thres - epsilon
  } else if (all(PI_centered > PI_lower_thres)) {
    PI_lower_thres = PI_lower_thres + epsilon
  }

  risk_group = 1 + (PI_centered > PI_lower_thres)
  risk_group = c("low_risk", "high_risk")[risk_group]
  risk_group = factor(risk_group, levels = c("low_risk", "high_risk"))
  risk_group
}

#--------------------------------------------------------------------
#Clinical files - use TCGAbiolinks via previous script
#--------------------------------------------------------------------

clin = readRDS("clin_data_lncs_new_variables_July19_tcgabiolinks_data.rds")
for(i in 1:length(clin)){
  print(i)
  d = clin[[i]]
  z=which(str_detect(colnames(d), "ENSG"))
  d=d[,-z]

  lncs_keep = filter(allCands, cancer %in% d$Cancer[1])$gene
  gene_exp=as.data.table(filter(rna, Cancer == d$Cancer[1]))
  z=which(colnames(gene_exp) %in% c(lncs_keep, "patient"))
  gene_exp = gene_exp[,..z]
  d = merge(d, gene_exp, by="patient")
  clin[[i]] = d
}

lgg = clin[[1]]
gbm = clin[[12]]
kirc = clin[[3]]

saveRDS(lgg, file="/u/kisaev/TCGA_lgg_wsubtype_info_biolinks.rds")
saveRDS(gbm, file="/u/kisaev/TCGA_gbm_wsubtype_info_biolinks.rds")
saveRDS(kirc, file="/u/kisaev/TCGA_kirc_wsubtype_info_biolinks.rds")

saveRDS(lgg, file="TCGA_lgg_wsubtype_info_biolinks.rds")
saveRDS(gbm, file="TCGA_gbm_wsubtype_info_biolinks.rds")

#--------LOOK AT ASSOCIATIONS BETWEEN EXPRESSION-------------------------------

#For each clinical variable -> xaxis is the clinical variable
#y-axis it each lncRNAs expression
#x-axis is continous if variable is continous such as age...

get_clin_lnc_cors = function(dtt){
  canc = dtt$Cancer[1]
  print(canc)
  print(dim(dtt))
  cancer_type = canc_conv$type[which(canc_conv$Cancer %in% dtt$Cancer)]
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

          z = which(is.na(new_dat_plot[,2]))
          if(!(length(z)==0)){
            new_dat_plot = new_dat_plot[-z,]
          }
          #get correlation results and save results into file
          cor = rcorr(new_dat_plot$lncRNA_exp, new_dat_plot$Clinical, "spearman")$r[2]
          pval_cor = rcorr(new_dat_plot$lncRNA_exp, new_dat_plot$Clinical, "spearman")$P[2]
          chisq_pval = "nochisq"
          kw_pval = "nokw"

          #how good of a predictor of survial is the clinical variable itself?
          z=which(colnames(rna) %in% c("patient", "OS", "OS.time"))
          surv_dat = rna[,..z]
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
            clin_vs_combo_anova = anova(cox_clin, both)[2,4]
          }

          if(length(z)==0){
          anov_pval = anova(cox_lnc, both)[2,4]
          clin_vs_combo_anova = anova(cox_clin, both)[2,4]
          }

          #add lnc, clinical Hazard Ratios and conordance of combined model

          row = c(canc, lnc, col, cor, pval_cor, "Ftest", chisq_pval, kw_pval,
          clin_pval, anov_pval, lnc_concordance, clin_concordance, hr, hr_clin, combo_concordance, clin_vs_combo_anova)
          names(row) = colnames(canc_col_results)

          #print(ggforest(both, main = paste(lnc, col, canc), data=new_dat_plot))

          canc_col_results = rbind(canc_col_results, row)
          #scatter plot
          if(pval_cor < 0.05){
          sp <- ggscatter(new_dat_plot, x = "Clinical", y = "lncRNA_exp",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE # Add confidence interval
          )
          # Add correlation coefficient
          sp = sp + stat_cor(method = "spearman") + theme_bw() + ggtitle(paste(cancer_type, allCands$gene_symbol[allCands$gene==lnc][1], col))
          print(sp)}#only plot if sig

        }
        }

        #if(is.na(test)[1]){
        if(length(which(is.na(test))) == length(test)){
        #boxplot

        #remove catgeories with less than 5 patients
        t = as.data.table(table(new_dat_plot[,2]))
        t = filter(t, N < 5)
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

        if(!(length(z1)==0)){
        new_dat_plot = new_dat_plot[-z1,]}

       z2 = which(new_dat_plot[,which(colnames(new_dat_plot) %in% col)] %in% c("[Unknown]", "[Not Available]",
       "#N/A", "[Not Evaluated]", "[Discrepancy]", "[Not Applicable]", "Unknown", "N/A", "NA", "Not Available", "Not performed",
        "Indeterminate", "[Not Available]", "[Unknown]", "Not Performed Clinical",
       "Performed but Not Available", "[Not submitted]"))

       if(!(length(z2)==0)){
        new_dat_plot = new_dat_plot[-z2,]}

       print(unique(new_dat_plot[,2]))

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
        z = which(colnames(rna) %in% c("patient", "OS", "OS.time"))
        surv_dat = rna[,..z]
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

        if((dim(table(new_dat_plot$lncRNA_tag)) > 1) & lncheck){ #& med_check){

        cox_lnc = coxph(Surv(OS.time, OS) ~ lncRNA_tag, data = new_dat_plot)
        cox_clin = coxph(Surv(OS.time, OS) ~ Clinical, data = new_dat_plot)
        both = coxph(Surv(OS.time, OS) ~ lncRNA_tag + Clinical, data = new_dat_plot)

        clin_pval = glance(cox_clin)[8]

        clin_concordance = glance(cox_clin)$concordance
        lnc_concordance = glance(cox_lnc)$concordance
        combo_concordance = glance(both)$concordance
        hr_clin = summary(cox_clin)$coefficients[2]
        anov_pval = anova(cox_lnc, both)[2,4]
        clin_vs_combo_anova = anova(cox_clin, both)[2,4]
        #print(ggforest(both, main = paste(lnc, col, canc), data=new_dat_plot))

#        check_comp2 = ((col == "clinical_stage") & (lnc == "ENSG00000259641"))
        check_comp2 = TRUE

        if(check_comp2){
          new_dat_plot$OS.time = new_dat_plot$OS.time/365
          cox_lnc = coxph(Surv(OS.time, OS) ~ Clinical, data = new_dat_plot)
          relRisk <- predict(cox_lnc, new_dat_plot, type="risk")   # relative risk
          new_dat_plot$rel_risk_clin_only = relRisk

          # split into two risk groups based on median
          PI_lower_thres = median(new_dat_plot$rel_risk_clin_only)
          PI_max_threshold = summary(new_dat_plot$rel_risk_clin_only)[5]
          new_dat_plot$risk_group_clin_only = get_median_risk_group(
            new_dat_plot$rel_risk_clin_only,
            PI_lower_thres)
          new_dat_plot$risk_group_clin_only = factor(new_dat_plot$risk_group_clin_only, levels=c("low_risk", "high_risk"))

          #get risk based models clinical plus lncRNA only
          cox_lnc = coxph(Surv(OS.time, OS) ~ Clinical + lncRNA_tag, data = new_dat_plot)
          relRisk <- predict(cox_lnc, new_dat_plot, type="risk")   # relative risk
          new_dat_plot$rel_risk_clin_plus_lncRNA = relRisk

          # split into two risk groups based on median
          PI_lower_thres = median(new_dat_plot$rel_risk_clin_plus_lncRNA)
          PI_max_threshold = summary(new_dat_plot$rel_risk_clin_plus_lncRNA)[5]
          new_dat_plot$risk_group_clin_plus_lncRNA = get_median_risk_group(
            new_dat_plot$rel_risk_clin_plus_lncRNA,
            PI_lower_thres)
          new_dat_plot$risk_group_clin_plus_lncRNA = factor(new_dat_plot$risk_group_clin_plus_lncRNA, levels=c("low_risk", "high_risk"))

          file_name = paste("/u/kisaev/", lnc, col, ".pdf", sep="_")

          #pdf(file_name, width=6, height=5)
          fit <- survfit(Surv(OS.time, OS) ~ risk_group_clin_only, data = new_dat_plot)
          cox_lnc = coxph(Surv(OS.time, OS) ~ risk_group_clin_only, data = new_dat_plot)
          hr=round(summary(cox_lnc)$coefficients[2], digits=3)
          pval=round(summary(cox_lnc)$coefficients[5], digits=15)
          lowrisksamps = table(new_dat_plot$risk_group_clin_only)[1]
          highrisksamps = table(new_dat_plot$risk_group_clin_only)[2]

          s1 <- ggsurvplot(fit ,
                  title = paste("HR=", hr, "waldpval=", pval, "riskhigh=", highrisksamps, "risklow=", lowrisksamps),
                  xlab = "Time (Years)",
                  font.main = c(7, "bold", "black"),
                  data = new_dat_plot,      # data used to fit survival curves.
                  pval = TRUE,             # show p-value of log-rank test.
                  conf.int = FALSE,        # show confidence intervals for
                  #xlim = c(0,8),
                  risk.table = FALSE,      # present narrower X axis, but not affect
                  break.time.by = 1,     # break X axis in time intervals by 500.
                  palette =c("#4DBBD5FF", "#E64B35FF"))
          print(s1)
          fit <- survfit(Surv(OS.time, OS) ~ risk_group_clin_plus_lncRNA, data = new_dat_plot)

          cox_lnc = coxph(Surv(OS.time, OS) ~ risk_group_clin_plus_lncRNA, data = new_dat_plot)
          hr=round(summary(cox_lnc)$coefficients[2], digits=3)
          pval=round(summary(cox_lnc)$coefficients[5], digits=15)
          lowrisksamps = table(new_dat_plot$risk_group_clin_plus_lncRNA)[1]
          highrisksamps = table(new_dat_plot$risk_group_clin_plus_lncRNA)[2]

          s2 <- ggsurvplot(fit ,
                  xlab = "Time (Years)",
                  font.main = c(7, "bold", "black"),
                  title = paste("HR=", hr, "waldpval=", pval, "riskhigh=", highrisksamps, "risklow=", lowrisksamps),
                  data = new_dat_plot,      # data used to fit survival curves.
                  pval = TRUE,             # show p-value of log-rank test.
                  conf.int = FALSE,        # show confidence intervals for
                  #xlim = c(0,8),
                  palette =c("#4DBBD5FF", "#E64B35FF"),
                  break.time.by = 1)#,     # break X axis in time intervals by 500.
          print(s2)
          #dev.off()
        }
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

        if(chisq_pval < 0.05){

        #p <- ggboxplot(new_dat_plot, x = "Clinical", y = "lncRNA_exp",
        #  color = "Clinical",
        #  title = paste(cancer_type, get_name(lnc), col),
        #  add = "jitter", ylab = "lncRNA expression",  ggtheme = theme_bw()) +
        #  stat_compare_means() + geom_hline(yintercept=med, linetype="dashed", color = "red") +
        #  stat_n_text(size=5)

        #p = ggpar(p,
        #  font.xtickslab = c(14,"plain", "black"),font.tickslab=c(14,"plain", "black"),
        #  xtickslab.rt = 55, legend="none")
        #print(p)} #only print plot if significant

        #barplot order patients by increasing lncRNA expression
        #plotting the lncRNA values smallest to largest in every subgroup, separate the subgroups by color
        #and facet. like this, except that this is very post-processed. This visual may help us better deal with zeroes than boxplots.

        new_dat_plot = as.data.table(new_dat_plot)
        new_dat_plot = new_dat_plot[order(lncRNA_exp)]
        new_dat_plot$patient = factor(new_dat_plot$patient, levels=new_dat_plot$patient)
        new_dat_plot$lncRNA_tag = factor(new_dat_plot$lncRNA_tag, levels=c(1,0))

         #facet_grid(Clinical ~ , scales = "free", space = "free") +
         g <- ggplot(new_dat_plot, aes(patient, lncRNA_exp)) +  geom_col(aes(fill = Clinical))+
           facet_grid(~lncRNA_tag+Clinical, space="free", scales="free") + theme_bw()+
           theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            strip.text.x = element_text(size = 3, colour = "black"),
            legend.position = "none")+ggtitle(allCands$gene_symbol[allCands$gene==lnc][1])

#           facet_wrap(~ lncRNA_tag+Clinical, scales = "free_x", nrow=1) + theme_bw() +

          print(g)

        }
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

    pdf(paste("/u/kisaev/", cancer_type, "clinical_plots.pdf", sep="_"), width=6, height=5)
    all_canc_lncs_results = llply(lncs, get_cor)
    all_canc_lncs_results = do.call(rbind.data.frame, all_canc_lncs_results)
    dev.off()
    return(all_canc_lncs_results)

  } #end get_clin_lnc_cors

#clin_data_lncs_cors = llply(clin_data_lncs, get_clin_lnc_cors)
clin_data_lncs = clin

d1 = get_clin_lnc_cors(clin_data_lncs[[1]])
d2 = get_clin_lnc_cors(clin_data_lncs[[2]])
d3 = get_clin_lnc_cors(clin_data_lncs[[3]])
d4 = get_clin_lnc_cors(clin_data_lncs[[4]])
d5 = get_clin_lnc_cors(clin_data_lncs[[5]])
d6 = get_clin_lnc_cors(clin_data_lncs[[6]])
d8 = get_clin_lnc_cors(clin_data_lncs[[8]])
d9 = get_clin_lnc_cors(clin_data_lncs[[9]])
d10 = get_clin_lnc_cors(clin_data_lncs[[10]])
d11 = get_clin_lnc_cors(clin_data_lncs[[11]])
d12 = get_clin_lnc_cors(clin_data_lncs[[12]])
d13 = get_clin_lnc_cors(clin_data_lncs[[13]])

#d7 = get_clin_lnc_cors(clin_data_lncs[[7]]) #uterine didnt work not enough data

all_clin = list(d1, d2,d3,d4,d5,d6, d8, d9,d10, d11,d12,d13)
saveRDS(all_clin, file="12_data_sets_biolinks_results.rds")
