set.seed(911)

source("/u/kisaev/lncRNAs_TCGA/NewTCGA_2017scripts/Thesis_final_scripts/Prognostic_lncs_more_detail_paper/final_scripts_2019/revisions_2020/load_data.R")
#setWD
setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNAs_2019_manuscript")

#-------------------------------------------------------------------
#------PCAWG DATA---------------------------------------------------
#-------------------------------------------------------------------

tcga_results1 = readRDS("TCGA_results_multivariate_results_Oct3.rds")
tcga_results1$data = "TCGA"
tcga_results1$combo = paste(tcga_results1$gene, tcga_results1$cancer, sep="_")
tcga_results1$clinical_only_concordance = NULL
tcga_results1$num_events = NULL
tcga_results1$perc_wevents = NULL
tcga_results1 = as.data.table(tcga_results1)
saveRDS(tcga_results1, file="final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")

#--------------------------------------------------------------------
#DONT RUN BELOW#
#--------------------------------------------------------------------

pcawg_clin = readRDS("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/top5_cancers_pcawg_thesis_jult2017/raw/Jan26_PCAWG_clinical")

#z = which(tcga_results1$combo %in% robust$combo)
#tcga_results1 = tcga_results1[z,]

#only the robust ones
#robust = readRDS(file="112_combos_robust_internal_validation_survival_lncRNAs_aug8.rds")
#robust = readRDS("148_combos_robust_5perc_increase_internal_validation_survival_lncRNAs_aug9.rds")
head(cands)
cands$combo = paste(cands$gene, cands$canc, sep="_")
#z = which(cands$combo %in% robust$combo)
#cands = cands[z,]

pcawg_data = readRDS("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNA_clinical_data_PCAWG_May2.rds")
pcawg_data$combo = paste(pcawg_data$canc, pcawg_data$histo)

pcawg_data$canc[pcawg_data$combo == "Uterus Adenocarcinoma, endometrioid"] = "Uterine Corpus Endometrial Carcinoma" 
pcawg_data$canc[pcawg_data$combo == "Kidney Adenocarcinoma, clear cell type"] = "Kidney renal clear cell carcinoma"
pcawg_data$canc[pcawg_data$combo == "Breast Infiltrating duct carcinoma"] = "Breast invasive carcinoma"
pcawg_data$canc[pcawg_data$combo == "Ovary Serous cystadenocarcinoma"] = "Ovarian serous cystadenocarcinoma"
pcawg_data$canc[pcawg_data$combo == "Pancreas Pancreatic ductal carcinoma"] = "Pancreatic adenocarcinoma" 
pcawg_data$canc[pcawg_data$combo == "Liver Hepatocellular carcinoma"] = "Liver hepatocellular carcinoma"
pcawg_data$canc[pcawg_data$combo == "Kidney Adenocarcinoma, papillary type"] = "Kidney renal papillary cell carcinoma"
pcawg_data$canc[pcawg_data$combo == "Lung Squamous cell carcinoma"] = "Lung squamous cell carcinoma"
pcawg_data$canc[pcawg_data$combo == "Lung Adenocarcinoma, invasive"] = "Lung adenocarcinoma"
pcawg_data$canc[pcawg_data$combo == "CNS Glioblastoma"] = "Glioblastoma multiforme "
pcawg_data$canc[pcawg_data$combo == "Esophagus Adenocarcinoma"] = "Esophageal carcinoma "
pcawg_data$canc[pcawg_data$canc == "Colon/Rectum"] = "Colon adenocarcinoma"
pcawg_data$canc[pcawg_data$canc == "Stomach"] = "Stomach adenocarcinoma"
pcawg_data$canc[pcawg_data$canc == "Head/Neck"] = "Head and Neck squamous cell carcinoma"
pcawg_data$canc[pcawg_data$canc == "Cervix"] = "Cervical squamous cell carcinoma and endocervical adenocarcinoma"

#remove those already used in TCGA cohort - aka those with a TCGA id 
ids = unique(pcawg_clin[,c("icgc_donor_id", "submitted_donor_id.y")])
colnames(ids) = c("patient", "TCGA_id")
pcawg_data = merge(pcawg_data, ids, by= "patient")
z = which(str_detect(pcawg_data$TCGA_id, "TCGA"))
pcawg_data = pcawg_data[-z,]

#add more pcawg people
cancers_tests = as.list(unique(tcga_results1$cancer[which(tcga_results1$cancer %in% pcawg_data$canc)]))

get_matched_data = function(cancer){
    dtt = subset(pcawg_data, canc == cancer)
    z = which(colnames(dtt) %in% c(as.character(cands$gene[cands$canc == dtt$canc[1]]), "canc", 
    "histo", "time", "status", "sex", "donor_age_at_diagnosis", "patient"))
    dtt = dtt[,z]
    if(nrow(dtt) >= 20){
    return(dtt)}
}

filtered_data = llply(cancers_tests, get_matched_data)
#remove nulls from list
filtered_data = Filter(Negate(is.null), filtered_data) #12 cancers left
getnames = ldply(filtered_data)
getnames =  unique(getnames$canc)
names(filtered_data) = getnames

all_pcawg_cancers = as.data.table(ldply(filtered_data))
all_pcawg_cancers = unique(all_pcawg_cancers[,c("canc", "status", "time", "sex", "histo", "donor_age_at_diagnosis", "patient")])
write.csv(all_pcawg_cancers, file="PCAWG_clinical_data_patients_used_external_survival_analysis.csv", quote=F, row.names=F)

add_tags = function(dtt){
  print(head(dtt))
  z = which(str_detect(colnames(dtt), "ENSG"))
  if(length(z)>1){
  medians = apply(dtt[,z], 2, median)}
  if(length(z)==1){
    medians = median(dtt[,z])
  }
  #add high low tag
  del = c()
  for(k in 1:length(medians)){
    print(k)
    med = medians[k]
    if(med ==0){
    #if median = 0 then anyone greater than zero is 1
    #but first check if non-zero group has a median higher than 0.1  
    l1 = which(dtt[,k+1] > 0)
    l2 = which(dtt[,k+1] ==0)
    check = median(dtt[l1,k+1])
    if((check >= 0.05) & (length(l1) >= 5)){
    dtt[l1,k+1] = 1
    dtt[l2, k+1] = 0
      }
    if((!(check >= 0.05)) | is.na(check)){
      del = c(del, names(medians)[k])
      }
    if(!((check >= 0.05) & (length(l1) >= 5))){
      del = c(del, names(medians)[k])
      }  
    }

    #check if meidan is greater than 0.1 otherwise save gene for deletion
    if(!(med ==0)){
    if(med >= 0.05){
    l1 = which(dtt[,k+1] >= med)
    l2 = which(dtt[,k+1] < med)
    dtt[l1,k+1] = 1
    dtt[l2, k+1] = 0
    }
    if(med < 0.05){
      del = c(del, names(medians)[k])
    }
    }
    }
  z = which(colnames(dtt) %in% del)
  if(!(length(z)==0)){
  dtt = dtt[,-z]}

  z = which(str_detect(colnames(dtt), "ENSG"))
  if(!(length(z)==0)){
  return(dtt)
  }
}

filtered_data_tagged = llply(filtered_data, add_tags, .progress="text")
filtered_data_tagged = Filter(Negate(is.null), filtered_data_tagged) #12 cancers left
getnames = ldply(filtered_data_tagged)
getnames =  unique(getnames$canc)
names(filtered_data_tagged) = getnames

get_survival_models = function(dtt){
  print(head(dtt))
  results_cox1 <- as.data.frame(matrix(ncol=17)) ; colnames(results_cox1) <- c("gene", "coef", "HR", "pval", "low95", "upper95", "cancer", 
    "lnc_test_ph", 'global_test_ph', "num_risk", "perc_risk", "median_nonzero", "sd_nonzero", "min_nonzero", "max_nonzero",
    "multi_model_concordance", "lnc_only_concordance")

  dat = dtt
  dat$canc = NULL
  dat$histo = NULL
  #dat$sex = NULL
  dat$status[dat$status=="alive"] =0
  dat$status[dat$status=="deceased"] =1

  z = which(str_detect(colnames(dat), "ENSG"))
  check_contrasts = function(col){
        check = dim(table(col))
        if(check >1){
          return("keep")
        }
  }
  keep = unlist(apply(dat, 2, check_contrasts))
  print(keep)

  dat = dat[,c(which(colnames(dat) %in% names(keep)))]

  dat$status = as.numeric(dat$status)
  dat$time = as.numeric(dat$time)
  dat$donor_age_at_diagnosis = as.numeric(dat$donor_age_at_diagnosis)
  num_genes = which(str_detect(colnames(dat), "ENSG"))

  for(i in 1:length(num_genes)){
  print(i)
  gene = num_genes[i]  
  k = which(!(str_detect(colnames(dat), "ENSG")))

  newdat = dat[,c(gene,k)]
  gene = colnames(newdat)[1]
  gene_name = colnames(newdat)[1]

  colnames(newdat)[1] = "gene"

  lncs = coxph(Surv(time, status)  ~ gene, data = newdat)
  test.ph <- cox.zph(lncs)
  lnc_test_ph = test.ph$table[1,3]
  global = test.ph$table[nrow(test.ph$table),3]

  if(!(is.na(lnc_test_ph))){

  #mutlvariate concordance 
  cmulti = "not_avail"
  clnconly = as.numeric(glance(lncs)[11])

  #calculate power
  #k is ratio of participants in group E (experimental group) compared to group C (controlgroup).
  kval = as.numeric(table(newdat[,1])[1]/table(newdat[,1])[2])
  #m = expected total number of events over both groups.
  mval = length(which(newdat$status == 1))

  #power = powerCT.default0(k=kval,m=mval, RR=2, alpha=0.05)

  hr = summary(lncs)$coefficients[1,c(1,2,5)][2]
  if(hr >1){
    risk_num = length(which(newdat[,1] == 1))
    perc = risk_num/nrow(newdat)
  }

  if(hr < 1){
    risk_num = length(which(newdat[,1] == 0))
    perc = risk_num/nrow(newdat)
  }
  

  #check that HR Confidence Interval is not infiinite 
  #check that at least 5 patients have expression greater than 0.1 
  lower95 = summary(lncs)[8][[1]][3]
  upper95 = summary(lncs)[8][[1]][4]

  check = ((!(lower95 == "Inf")) & (!(upper95=="Inf")))
  if(check){
  colnames(newdat)[1] = "gene"
  newdat$time = newdat$time/365
  newdat$gene = factor(newdat$gene, levels=c(1,0))
  
  fit <- survfit(Surv(time, status) ~ gene, data = newdat)
          s <- ggsurvplot(
          title = paste(gene_name, canc_conv$type[canc_conv$Cancer == dtt$canc[1]][1], "HR =", round(hr, digits=2)),
          fit, 
          xlab = "Time (Years)", 
          #surv.median.line = "hv",
          font.main = c(14, "bold", "black"),
          font.x = c(12, "plain", "black"),
          font.y = c(12, "plain", "black"),
          font.tickslab = c(11, "plain", "black"),
          font.legend = 10,
          risk.table.fontsize = 5, 
          legend.labs = c("High Expression", "Low Expression"),             # survfit object with calculated statistics.
          data = newdat,      # data used to fit survival curves. 
          risk.table = TRUE,       # show risk table.
          legend = "right", 
          pval = TRUE,             # show p-value of log-rank test.
          conf.int = FALSE,        # show confidence intervals for 
                            # point estimaes of survival curves.
          xlim = c(0,10),        # present narrower X axis, but not affect
                            # survival estimates.
          break.time.by = 1,     # break X axis in time intervals by 500.
          #palette = colorRampPalette(mypal)(14), 
          #palette = mypal[c(4,1)],
          palette = "npg", 
          #ggtheme = theme_bw(), # customize plot and risk table with a theme.
          risk.table.y.text.col = T, # colour risk table text annotations.
          risk.table.y.text = FALSE # show bars instead of names in text annotations
                            # in legend of risk table
          )
          #print(s)      


      #generate boxplot 
      z = which(names(filtered_data) == dtt$canc[1])
      exp_data = filtered_data[[z]]
      exp_data$patient = rownames(exp_data)
      exp_data = exp_data[,which(colnames(exp_data) %in% c(gene, "patient"))]
      newdat$patient = rownames(newdat)
      exp_data = merge(exp_data, newdat, by="patient")
      colnames(exp_data)[2] = "geneexp"
      
      #get median non-zero expression and SD 
      if(perc < 0.45){
      median_nonzero = median(exp_data[which(!(exp_data[,2] == 0)),2])
      sd_nonzero = sd(exp_data[which(!(exp_data[,2] == 0)),2])
      min_nonzero = min(exp_data[which(!(exp_data[,2] == 0)),2])
      max_nonzero = max(exp_data[which(!(exp_data[,2] == 0)),2])
      } else if (perc > 0.55){
      median_nonzero = median(exp_data[which(!(exp_data[,2] == 0)),2])
      sd_nonzero = sd(exp_data[which(!(exp_data[,2] == 0)),2])
      min_nonzero = min(exp_data[which(!(exp_data[,2] == 0)),2])
      max_nonzero = max(exp_data[which(!(exp_data[,2] == 0)),2])
      } else {
      median_nonzero = "notavail"
      sd_nonzero = "notavail"
      min_nonzero = "notavail"
      max_nonzero = "notavail"
      }

    #check at least 5 pateints witih exp > 5
    num_pats = length(which(exp_data$geneexp >= 0.1))
    pval = as.numeric(glance(lncs)[4])

    if(num_pats >= (0.01 * nrow(exp_data))){
    print(s)
    row <- c(gene_name, summary(lncs)$coefficients[1,c(1,2)], pval,  summary(lncs)$conf.int[1,c(3,4)], dtt$canc[1], 
    lnc_test_ph, global, risk_num, perc, median_nonzero,
    sd_nonzero,
    min_nonzero,
    max_nonzero, cmulti, clnconly)

    names(row) <- names(results_cox1) 
    results_cox1 = rbind(results_cox1, row)

    #make density plot using FPKM-UQ values and logged values
    #visualize bimodality 
    exp_data$median=""
    exp_data$median[exp_data$gene==0] = "Low"
    exp_data$median[exp_data$gene==1] = "High"
    exp_data$median = factor(exp_data$median, levels=c("High", "Low"))
   
    gg <- ggplot(exp_data)
    gg <- gg + geom_density(aes(x=geneexp, y=..scaled.., fill=median), alpha=1/2)
    gg <- gg + theme_bw() + ggtitle(paste(gene, "Expression", dtt$Cancer[1] , sep=" "))
    #print(gg)

        p <- ggboxplot(exp_data, x = "gene", y = "geneexp",
          color = "gene",
         palette = mypal[c(4,1)], title = paste(gene_name, "Expression", dtt$canc[1] , sep=" "), 
          add = "jitter", ylab = "FPKM-UQ",  ggtheme = theme_bw())
        # Change method
       p = p + stat_compare_means(method = "wilcox.test")
       #print(p)


   #get measure of bimodality 
   #SIBER(y=exp_data$geneexp, model='LN')

    exp_data$geneexp = log1p(exp_data$geneexp)
  
    #gg <- ggplot(exp_data)
    #gg <- gg + geom_density(aes(x=geneexp, y=..scaled.., fill=median), alpha=1/2)
    #gg <- gg + theme_bw() + ggtitle(paste(gene, "Expression", dtt$Cancer[1] , sep=" "))
    #print(gg)

      #exp_data$geneexp = log1p(exp_data$geneexp)
      p <- ggboxplot(exp_data, x = "gene", y = "geneexp",
          color = "gene",
         palette = mypal[c(4,1)], title = paste(gene_name, "Expression", dtt$canc[1] , sep=" "), 
          add = "jitter", ylab = "log1p(FPKM-UQ)",  ggtheme = theme_bw())
        # Change method
       p = p + stat_compare_means(method = "wilcox.test")
       #print(p)

  p <- ggboxplot(exp_data, x = "median", y = "geneexp",
          color = "median",
         palette = mypal[c(4,1)], title = paste(gene_name, "Expression", canc_conv$type[canc_conv$Cancer == dtt$canc[1]][1], sep=" "), 
          add = "jitter", ylab = "log1p(FPKM-UQ)",  ggtheme = theme_classic())
        # Change method
  p = p + stat_compare_means(method = "wilcox.test") + stat_n_text() + scale_color_npg() 
  print(p)


}
}
}
}

results_cox1 = results_cox1[-1,]
if(!(dim(results_cox1)[1] == 0)){

#get clinical model
z = which(!(str_detect(colnames(dat), "ENSG")))
dat = dat[,z]
clinmodel = coxph(Surv(time, status)  ~ ., data = dat)
conc = unlist(glance(clinmodel)[11])

row = results_cox1[1,]
row = c("clin", rep("na", 15), conc, "na")

results_cox1 = rbind(results_cox1, row)
results_cox1$cancer = results_cox1$cancer[1]
print(results_cox1)

return(results_cox1)
}
}

#pdf("PCAWG_validating_individual_TCGA_candidates_survival_plots_final_cands_May3rd_full_lifespan.pdf")
pdf("PCAWG_validating_individual_TCGA_candidates_survival_plots_Nov262019_five_year_survival.pdf", width=6, height=5)
pcawg_results = llply(filtered_data_tagged, get_survival_models, .progress="text")
dev.off()

#all coxph results for lcnRNAs in TCGA (these p-values came from including clinical variables in the models)
pcawg_results1 = ldply(pcawg_results, data.frame)

#z = which(pcawg_results1$lnc_only_concordance == "na")
#clinical_concs = pcawg_results1[z,]
#pcawg_results1 = pcawg_results1[-z,]
pcawg_results1$pval = as.numeric(pcawg_results1$pval)
pcawg_results1$fdr_pval = p.adjust(pcawg_results1$pval, method="fdr")
z = which(pcawg_results1$gene == "clin")
clin= pcawg_results1[z,]
pcawg_results1 = pcawg_results1[-z,]

#combine results from TCGA and PCAWG
pcawg_results1$data = "PCAWG"
pcawg_results1$num_risk = as.numeric(pcawg_results1$num_risk)
#pcawg_results1 = filter(pcawg_results1, num_risk >=5)
pcawg_results1 = as.data.table(pcawg_results1)
filter(pcawg_results1, pval < 0.05)

#z = which(pcawg_results1$upper95 == "Inf")
#pcawg_results1 = pcawg_results1[-z,]

pcawg_results1$combo = paste(pcawg_results1$gene, pcawg_results1$cancer, sep="_")

saveRDS(pcawg_results1, file="PCAWG_external_validation_lncCands_feb19.rds")

#plot power versus Hazard Ratio 
pcawg_results1$pval = as.numeric(pcawg_results1$pval)
write.table(pcawg_results1, file="all_lncRNAs_evaluated_by_PCAWG_KI_Nov262019.txt", quote=F, row.names=F, sep="\t")

#all-results
#tcga_results1$lnc_test_ph =NULL
#tcga_results1$global_test_ph = NULL
tcga_results1$groupy = NULL
tcga_results1$gene_name = NULL
tcga_results1$lnc_better = NULL
pcawg_results1[,1] = NULL

all_results = as.data.table(rbind(tcga_results1, pcawg_results1))
all_results = all_results[order(gene, cancer)]

#add more info on lncRNAs
lnc_info = read.csv("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/fantom_genebased_evidence_supp_table_17.csv")
lnc_info = lnc_info[which(str_detect(lnc_info$CAT_geneID, "ENSG")),]

#shorten the gene names 
extract3 <- function(row){
  gene <- as.character(row[[1]])
  ens <- gsub("\\..*","",gene)
  return(ens)
}
lnc_info[,1] <- apply(lnc_info[,1:2], 1, extract3) #5049 lncRNAs 

colnames(lnc_info)[1] = "gene"

all_results = merge(all_results, lnc_info, by="gene")
all_results = all_results[order(data, pval)]
all_results$best_cell_trait_fdr = NULL
all_results$CAT_browser_link = NULL

write.table(all_results, file="TCGA_PCAWG_prognostic_lncRNAs.txt", quote=F, row.names=F, sep="\t")

#-----check which actually match---------------------------------------------------------------------------

all_results_orig = fread("TCGA_PCAWG_prognostic_lncRNAs.txt")

all_results = all_results_orig[!duplicated(all_results_orig), ]
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
matches$combo = paste(matches$gene, matches$cancer, sep="_")

ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
#z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
#ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,6]))
ucsc <- ucsc[-z,]
colnames(ucsc)[6] = "gene"
matches = merge(matches, ucsc, by=c("gene"))

#remove weird KIRC one (already did above by filtering)
#z = which(matches$gene == "ENSG00000250360")
#matches = matches[-z,]
write.table(matches, file="5_unique_lncNRAs_validate_PCAWG.txt", quote=F, row.names=F, sep=";")

#write.table(matches, file="4_unique_lncNRAs_validate_PCAWG.txt", quote=F, row.names=F, sep=";")
matches = as.data.frame(matches)
matches = matches[,1:13]
pdf("5_unique_lncNRAs_validate_PCAWG.pdf", width=24)
p<-tableGrob(matches)
grid.arrange(p)
dev.off()
matches$combo = paste(matches$gene, matches$cancer, sep="_")

#Clean up and save as spreadsheet
all_results_orig$combo = paste(all_results_orig$gene, all_results_orig$cancer, sep="_")
z = which(all_results_orig$combo %in% matches$combo)
all_results_orig$top_pcawg_val = ""
all_results_orig$top_pcawg_val[z] = "YES"
all_results_orig$combo = NULL
all_results_orig$best_eQTL_mRNA_fdr = NULL

all_results_orig$coef = as.numeric(all_results_orig$coef)
all_results_orig$coef = round(all_results_orig$coef, digits=4)

all_results_orig$HR = as.numeric(all_results_orig$HR)
all_results_orig$HR = round(all_results_orig$HR, digits=4)

all_results_orig$pval = as.numeric(all_results_orig$pval)
all_results_orig$pval = round(all_results_orig$pval, digits=6)

all_results_orig$low95 = as.numeric(all_results_orig$low95)
all_results_orig$low95 = round(all_results_orig$low95, digits=4)

all_results_orig$upper95 = as.numeric(all_results_orig$upper95)
all_results_orig$upper95 = round(all_results_orig$upper95, digits=4)

all_results_orig$lnc_test_ph = as.numeric(all_results_orig$lnc_test_ph)
all_results_orig$lnc_test_ph = round(all_results_orig$lnc_test_ph, digits=4)

all_results_orig$global_test_ph = as.numeric(all_results_orig$global_test_ph)
all_results_orig$global_test_ph = round(all_results_orig$global_test_ph, digits=4)

all_results_orig$perc_risk = as.numeric(all_results_orig$perc_risk)
all_results_orig$perc_risk = round(all_results_orig$perc_risk, digits=4)

all_results_orig$fdr_pval = as.numeric(all_results_orig$fdr_pval)
#all_results_orig$fdr_pval = round(all_results_orig$fdr_pval, digits=4)

all_results_orig$combo = paste(all_results_orig$gene, all_results_orig$cancer, sep="_")

#add info on number of events and % of events in each cancer type
tcga = readRDS("TCGA_results_multivariate_results_Oct3.rds")
#tcga = unique(tcga[,c("cancer", "perc_risk")])
#all_results_orig = merge(all_results_orig, tcga, by="cancer")
write.csv(all_results_orig, file="168_lncRNA_cancers_combos_22_cancer_types_feb19.csv", quote=F, row.names=F)
#write.csv(all_results_orig, file="112_lncRNA_cancers_combos_22_cancer_types_aug8.csv", quote=F, row.names=F)
saveRDS(all_results_orig, file="final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")

write.csv(all_results_orig, file="pcawg_results_tcga_results_file_supp_file_KI.csv", quote=F, row.names=F)



