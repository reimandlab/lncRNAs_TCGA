source("/u/kisaev/lncRNAs_TCGA/NewTCGA_2017scripts/Thesis_final_scripts/Prognostic_lncs_more_detail_paper/final_scripts_2019/revisions_2020/load_data.R")

#------FEATURES-----------------------------------------------------
setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNAs_2019_manuscript")

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = subset(allCands, data == "TCGA") #173 unique lncRNA-cancer combos, #166 unique lncRNAs
allCands$combo = unique(paste(allCands$gene, allCands$cancer, sep="_"))

library(corrplot)

get_median_risk_group = function(PI_centered, PI_lower_thres, epsilon = 1e-16) {

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

#results from cleaned up analysis
res = readRDS("process_tcga_biolinks_results_for_plotting.rds")
colnames(res)[2] = "lncRNA"
res$pairs = paste(res$lncRNA, res$type, res$colname, sep="_")
pairs = unique(res$pairs)
require(vcd)

get_name=function(g){
  z=which(allCands$gene == g)
  name=allCands$gene_symbol[z]
  name=name[1]
  return(name)
}

#--------LOOK AT ASSOCIATIONS BETWEEN EXPRESSION-------------------------------

#For each clinical variable -> xaxis is the clinical variable
#y-axis it each lncRNAs expression
#x-axis is continous if variable is continous such as age...

get_clin_lnc_plots = function(dtt){
  canc = dtt$Cancer[1]
  print(canc)
  print(dim(dtt))
  cancer_type = canc_conv$type[which(canc_conv$Cancer %in% dtt$Cancer)]
  #get lncs
  z = which(str_detect(colnames(dtt), "ENSG"))
  lncs = colnames(dtt)[z]
  z = which(lncs %in% res$lnc)
  if(!(length(z)==0)){
  lncs=lncs[z]

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

    lnc_id=lnc
    res_dat = as.data.table(filter(res, type==cancer_type, lnc==lnc_id))
    cols_dat = which(colnames(new_dat) %in% res_dat$colname)
    for(i in cols_dat){
      print(i)
      col = colnames(new_dat)[i]
      new_dat_plot = new_dat[,c("patient", col, lnc, "lncRNA_tag", "OS", "OS.time")]
      #palette
      colourCount = length(unique(new_dat_plot$Clinical))
      getPalette = colorRampPalette(brewer.pal(9, "Set1"))
      new_dat_plot[,3] = log1p(new_dat_plot[,3])
      med = median(new_dat_plot[,3])
      colnames(new_dat_plot)[3] = "lncRNA_exp"

      new_dat_plot[,2] = as.character(new_dat_plot[,2])

      z1 = which(is.na(new_dat_plot[,which(colnames(new_dat_plot) %in% col)]))

      if(!(length(z1)==0)){
        new_dat_plot = new_dat_plot[-z1,]}

      z2 = which(new_dat_plot[,which(colnames(new_dat_plot) %in% col)] %in% c("[Unknown]", "[Not Available]",
       "#N/A", "[Not Evaluated]", "[Discrepancy]", "[Not Applicable]", "Unknown", "N/A", "NA", "Not Available", "Not performed",
        "Indeterminate", "[Not Available]", "[Unknown]", "Not Performed Clinical",
       "Performed but Not Available", "[Not submitted]"))

      if(!(length(z2)==0)){
        new_dat_plot = new_dat_plot[-z2,]}

      unq = length(unique(new_dat_plot[,2]))
      colnames(new_dat_plot)[2] = "Clinical"

      print(unique(new_dat_plot[,2]))

      new_dat_plot = as.data.table(new_dat_plot)
      new_dat_plot = new_dat_plot[order(lncRNA_exp)]
      new_dat_plot$patient = factor(new_dat_plot$patient, levels=new_dat_plot$patient)
      new_dat_plot$lncRNA_tag = factor(new_dat_plot$lncRNA_tag, levels=c(1,0))

      #if numerical variable make correlation plot
      if(length(unique(new_dat_plot$Clinical)) > 20){
      new_dat_plot$Clinical = as.numeric(new_dat_plot$Clinical)
      g =  ggscatter(new_dat_plot, x = "Clinical", y = "lncRNA_exp",
          color = "black", shape = 21, size = 3, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "spearman", label.sep = "\n"))+ggtitle(paste(get_name(lnc), cancer_type, colnames(new_dat)[i]))
      print(g)
      }

      #if categorial variable make facetted barplot
      if(!(length(unique(new_dat_plot$Clinical)) > 20)){

        check_comp1 = ((col == "X1p.19q.codeletion") & (lnc == "ENSG00000239552"))
        check_comp2 = ((col == "clinical_stage") & (lnc == "ENSG00000259641"))

        if((check_comp1 | check_comp2)){
          #get risk based models clinical only
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

          pdf(file_name, width=6, height=5)
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
          dev.off()
        }

      g <- ggplot(new_dat_plot, aes(patient, lncRNA_exp)) +  geom_col(aes(fill = Clinical))+
           facet_grid(~lncRNA_tag+Clinical, space="free", scales="free") + theme_bw()+
           theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            strip.text.x = element_text(size = 3, colour = "black"),
            legend.position = "none")+ggtitle(paste(get_name(lnc), cancer_type, colnames(new_dat)[i]))

      print(g)
      tab=table(new_dat_plot$Clinical, new_dat_plot$lncRNA_tag)
      mo = mosaic(~ Clinical + lncRNA_tag, data = new_dat_plot,
      main = paste(get_name(lnc), cancer_type, colnames(new_dat)[i]), shade = TRUE, legend = TRUE)
      print(mo)

      chi=t(table(new_dat_plot$Clinical, new_dat_plot$lncRNA_tag))
      chisq=chisq.test(chi)
      #residuals = corrplot(chisq$residuals, is.cor = FALSE, col=c("blue", "red"))
      #print(residuals)

      #also print boxplot
      box = ggboxplot(new_dat_plot, "Clinical", "lncRNA_exp",
        color = "lncRNA_tag", palette =c("#FC4E07", "#00AFBB"),
        add = "jitter") +ggtitle(paste(get_name(lnc), cancer_type, colnames(new_dat)[i]))
      box = box + stat_compare_means(aes(group = lncRNA_tag), label = "p.format")+stat_n_text()
      print(box)

      #also print KM plot
      new_dat_plot$OS.time = new_dat_plot$OS.time/365

      fit <- survfit(Surv(OS.time, OS) ~ lncRNA_tag + Clinical, data = new_dat_plot)

      s <- ggsurvplot(
          title = paste(get_name(lnc), cancer_type, colnames(new_dat)[i]),
          fit,
          xlab = "Time (Years)",
          font.main = c(14, "bold", "black"),
          font.x = c(12, "plain", "black"),
         font.y = c(12, "plain", "black"),
          font.tickslab = c(11, "plain", "black"),
          font.legend = 4,
          risk.table.fontsize = 2,
          data = new_dat_plot,      # data used to fit survival curves.
          risk.table = TRUE,       # show risk table.
         legend = "right",
          pval = TRUE,             # show p-value of log-rank test.
          conf.int = FALSE,        # show confidence intervals for
                            # point estimaes of survival curves.
          xlim = c(0,10),        # present narrower X axis, but not affect
                            # survival estimates.
          break.time.by = 1,     # break X axis in time intervals by 500.
          palette = "npg",
          risk.table.y.text.col = T, # colour risk table text annotations.
          risk.table.y.text = FALSE # show bars instead of names in text annotations
                            # in legend of risk table
          )

      print(s)

      fit <- survfit(Surv(OS.time, OS) ~ lncRNA_tag, data = new_dat_plot)
      s2 <- ggsurvplot(fit ,
              xlab = "Time (Years)",
              data = new_dat_plot,      # data used to fit survival curves.
              facet.by = "Clinical",
              pval = TRUE,             # show p-value of log-rank test.
              conf.int = FALSE,        # show confidence intervals for
              xlim = c(0,10),
              palette ="npg",
              break.time.by = 1)#,     # break X axis in time intervals by 500.

      print(s2)

      fit <- survfit(Surv(OS.time, OS) ~ Clinical, data = new_dat_plot)
      s3 <- ggsurvplot(fit ,
              xlab = "Time (Years)",
              data = new_dat_plot,      # data used to fit survival curves.
              pval = TRUE,             # show p-value of log-rank test.
              conf.int = FALSE,        # show confidence intervals for
              xlim = c(0,10),
              palette ="npg",
              break.time.by = 1)#,
      print(s3)

      }

    }#all cols dat
  }#get_cor

  llply(lncs, get_cor)

  }#only run if lncs in final clinical res dataset
}#end get_clin_lnc_plots

pdf("/u/kisaev/Jan2021/Figure4_individual_plots_clinical_plots.pdf", width=6, height=5)
llply(clin, get_clin_lnc_plots, .progress="text")
dev.off()
