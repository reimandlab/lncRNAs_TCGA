#------FUNCTIONS----------------------------------------------------

#######
##[1]##-------------------------------------------------------------
#######

#get gene ID from name
#input: GeneName
#output: GeneID

get_name_pcg = function(pcg){
        z = which(hg38$ensgene == pcg)
  if(length(z)>1){
    z = z[1]
  }
  return(hg38$symbol[z])
}

get_ensg_pcg = function(pcg){
  z = which(ucsc$hg19.ensemblToGeneName.value == pcg)
  if(length(z)>1){
    z = z[1]
  }
  return(ucsc$hg19.ensGene.name2[z])
}

get_ensg = function(lnc){
  z = which(fantom$CAT_geneName == lnc)
  return(fantom$CAT_geneID[z])
}

get_name = function(ensg){
    z = which(fantom$CAT_geneID == ensg)
    return(fantom$CAT_geneName[z][1])
}

#######
##[2]##-------------------------------------------------------------
#######

#get lncRNA expression across cancers plot
#input: lncRNA id
#output: boxplot of expression across cancer types

#test case: pca3 - ENSG00000225937
#mypal = readRDS("best_pal.rds")

get_exp_plots = function(lnc){
        dat = rna[,which(colnames(rna) %in% c("type", lnc))]
        colnames(dat)[1] = "lncRNAExp"
        dat$lncRNAExp = log1p(dat$lncRN)
        dat = as.data.table(dat)
        dat = dat[order(lncRNAExp)]

        #get order of cancer types by median
        sum = as.data.table(dat %>% dplyr::group_by(type) %>%
                dplyr::summarize(median=quantile(lncRNAExp, 0.75)))
        sum = sum[order(median)]

        dat$type = factor(dat$type, levels=unique(sum$type))
        #p = ggboxplot(dat, x="type", y = "lncRNAExp", title=lnc, error.plot = "linerange")+
        #rotate_x_text(45)

        #try density plot
        p = ggdensity(dat, x = "lncRNAExp", title=paste(lnc, get_name(lnc)), color="type", palette=mypal5)+
  xlab("log1p(FPKM)")
        print(p)
        print("done plot")

}


get_num_peaks_den = function(lnc){
  dat = rna[,which(colnames(rna) %in% c("type", lnc))]
  colnames(dat)[1] = "lncRNAExp"
  dat$lncRNAExp = log1p(dat$lncRN)
  dat = as.data.table(dat)
  dat = dat[order(lncRNAExp)]

  #get order of cancer types by median
  sum = as.data.table(dat %>% dplyr::group_by(type) %>%
    dplyr::summarize(median=median(lncRNAExp)))
  sum = sum[order(median)]

  dat$type = factor(dat$type, levels=unique(sum$type))
  #p = ggboxplot(dat, x="type", y = "lncRNAExp", title=lnc, error.plot = "linerange")+
  #rotate_x_text(45)

  #how many peaks within each cancer?
  types = as.character(unique(dat$type))
  get_peaks = function(canc){
     D = dat$lncRNAExp[dat$type == canc]

     firstquantile = summary(D)[2]
     thirdquantile = summary(D)[5]
     med = median(D)

     #some data
     d <- density(D)

     #make it a time series
     ts_y<-ts(d$y)

     #calculate turning points (extrema)
     require(pastecs)
     tp<-turnpoints(ts_y)
     print(canc)
     row = c(lnc, canc, firstquantile, med, thirdquantile, length(d$x[tp$tppos]))
     names(row) = c("lnc", "cancer", "firstQ", "median", "thirdQ", "numPeaks")
     return(row)
  }
  canc_peaks = llply(types, get_peaks)
  canc_peaks = ldply(canc_peaks)
  canc_peaks = as.data.table(canc_peaks)
  canc_peaks$firstQ = as.numeric(canc_peaks$firstQ)
  canc_peaks$median = as.numeric(canc_peaks$median)
  canc_peaks$thirdQ = as.numeric(canc_peaks$thirdQ)
  canc_peaks = canc_peaks[order(-firstQ, -median, -thirdQ)]
  canc_peaks$med_count[canc_peaks$median == 0] = 0
  canc_peaks$med_count[canc_peaks$median > 0] = 1
  canc_peaks$firstQ_count[canc_peaks$firstQ == 0] = 0
  canc_peaks$firstQ_count[canc_peaks$firstQ > 0] = 1
  canc_peaks$thirdQ_count[canc_peaks$thirdQ == 0] = 0
  canc_peaks$thirdQ_count[canc_peaks$thirdQ > 0] = 1
  s = as.data.table(table(canc_peaks$firstQ_count, canc_peaks$med_count, canc_peaks$thirdQ_count))
  s = as.data.table(filter(s, N >0))
  if(dim(s)[1] > 1){
    lnc_type = "bimodal"
  }
   if(dim(s)[1] == 1){
    lnc_type = "unimodal"
  }
  return(c(lnc, lnc_type))
}

#pdf("pca3_expression_across_cancers.pdf")
#get_exp_plots("ENSG00000225937")
#get_num_peaks_den("ENSG00000225937")
#dev.off()


#######
##[x]##-------------------------------------------------------------
#######

#identify which lncRNAs are bimodal and which ones are
#highly expressed in all cancers

#get peaks
#density(data$V2)$x[which.max(density(data$V2)$y)]



#######
##[3]##-------------------------------------------------------------
#######

#check how many patients have expression value X FPKM
#(because we want at least 15 patients with 100FPKM expression)

get_num_pats = function(lnc, canc, exp_cut){
        z = which(colnames(rna) %in% c("Cancer", lnc))
        dat = rna[,..z]
        dat = dat[which(dat$Cancer == canc),]
        colnames(dat)[1] = "lncRNAExp"
        num_pats = length(which(dat[,1] > exp_cut))
        if(num_pats >= 15){
                return("great success")
        }
        if(num_pats < 15){
                return("problem")
        }
}

#######
##[4]##-------------------------------------------------------------
#######

#check expression of PCG how it's different between a lncRNA's risk group
#input: lncRNA, pcg, cancer

get_pcg_enrich = function(lnc, pcg, canc){
        dat = all[,which(colnames(all) %in% c("Cancer", lnc, pcg, "patient"))]
        dat = dat[which(dat$Cancer == canc),]
        colnames(dat)[which(colnames(dat)==lnc)] = "lncRNAExp"
        colnames(dat)[which(colnames(dat)==pcg)] = "pcgExp"

        #add risk group
        lnc_hr = as.numeric(allCands$HR[which(allCands$combo == paste(lnc, canc, sep="_"))])

        #get med
        med = median(dat$lncRNAExp)
        dat$lnc_median = ""
       if(med ==0){
        #if median = 0 then anyone greater than zero is 1
        l1 = which(dat$lncRNAExp > 0)
        l2 = which(dat$lncRNAExp ==0)
        dat$lnc_median[l1] = 1
        dat$lnc_median[l2] = 0
        }

      if(!(med ==0)){
        l1 = which(dat$lncRNAExp >= med)
        l2 = which(dat$lncRNAExp < med)
        dat$lnc_median[l1] = 1
        dat$lnc_median[l2] = 0
    }

   if(lnc_hr > 1){
        dat$lnc_median[dat$lnc_median == 1] = "risk"
        dat$lnc_median[dat$lnc_median == 0] = "less_risk"
   }

    if(lnc_hr < 1){
        dat$lnc_median[dat$lnc_median == 0] = "risk"
        dat$lnc_median[dat$lnc_median == 1] = "less_risk"
   }

   mean_diff = round(as.numeric(coexp$mean_diff[which(coexp$combo2 == paste(lnc, canc, pcg, sep="_"))]), digits=4)
   dat$pcgExp = log1p(dat$pcgExp)
   g = ggboxplot(dat, x = "lnc_median", y="pcgExp", title=paste(lnc, pcg, mean_diff), fill="lnc_median") +
        stat_compare_means()
   print(g)
}


#######
##[5]##-------------------------------------------------------------
#######

#gene = "ENSG00000225937" #PCA3
#cancer = "PAAD"
pal = pal_npg("nrc")(10)

###EASY WAY TO MAKE KM PLOT
get_km_plot_os = function(gene, cancer){
  all_g = all
  all_g = as.data.frame(all_g)
  z = which(colnames(all) %in% c("type", gene, "OS", "OS.time"))
  dat = all_g[,z]
  z = which(str_detect(colnames(dat), "ENSG"))
  if(!(length(z)==0)){
  colnames(dat)[z] = "gene"
  dat = subset(dat, type == cancer)
  #split patients
  med = median(dat$gene)
  #add high low tag
    if(med ==0){
    #if median = 0 then anyone greater than zero is 1
    l1 = which(dat[,z] > 0)
    l2 = which(dat[,z] ==0)
    dat[l1,z] = 1
    dat[l2, z] = 0
    }

    if(!(med ==0)){
    l1 = which(dat[,z] >= med)
    l2 = which(dat[,z] < med)
    dat[l1,z] = 1
    dat[l2, z] = 0
    }

  dat$OS = as.numeric(dat$OS)
  dat$OS.time = as.numeric(dat$OS.time)
  dat$OS.time = dat$OS.time/365
  dat$gene = factor(dat$gene, levels = c(0,1))
  gene_name = get_name(gene)
  if(is.na(gene_name)){
    gene_name = get_name_pcg(gene)
  }
  cox_mod = coxph(Surv(OS.time, OS) ~ gene, data = dat)
  print(glance(cox_mod)$concordance)
  conc = round(glance(cox_mod)$concordance, digits=2)
  hr = summary(cox_mod)$coefficients[2]
  dtt = dat

  fit <- survfit(Surv(OS.time, OS) ~ gene, data = dtt)

  s <- ggsurvplot(
          title = paste(get_name(gene), cancer, "HR =", round(hr, digits=2)),
          fit,
          ylab = "Survival Probability" ,
          xlab = "Time (Years)",
          #surv.median.line = "hv",
          font.main = c(4, "bold", "black"),
          font.x = c(3, "plain", "black"),
          font.y = c(3, "plain", "black"),
          font.tickslab = c(3, "plain", "black"),
          font.legend = 3,
          risk.table.fontsize = 1,
          #legend.labs = c("High Expression", "Low Expression"),             # survfit object with calculated statistics.
          data = dtt,      # data used to fit survival curves.
          risk.table = FALSE,       # show risk table.
          legend = "right",
          pval = TRUE,             # show p-value of log-rank test.
          conf.int = FALSE,        # show confidence intervals for
          xlim = c(0,10),        # present narrower X axis, but not affect
          break.time.by = 1,     # break X axis in time intervals by 500.
          palette = pal[c(2,1)])
          #risk.table.y.text.col = T, # colour risk table text annotations.
          #risk.table.y.text = FALSE )
          return(s)
}
}

get_km_plot_pfi = function(gene, cancer){
  all_g = all
  all_g = as.data.frame(all_g)
  z = which(colnames(all) %in% c("type", gene, "PFI", "PFI.time"))
  dat = all_g[,z]
  z = which(str_detect(colnames(dat), "ENSG"))
  if(!(length(z)==0)){
  colnames(dat)[z] = "gene"
  dat = subset(dat, type == cancer)
  #split patients
  med = median(dat$gene)

  #add high low tag
    if(med ==0){
    #if median = 0 then anyone greater than zero is 1
    l1 = which(dat[,z] > 0)
    l2 = which(dat[,z] ==0)
    dat[l1,z] = 1
    dat[l2, z] = 0
    }

    if(!(med ==0)){
    l1 = which(dat[,z] >= med)
    l2 = which(dat[,z] < med)
    dat[l1,z] = 1
    dat[l2, z] = 0
    }

  dat$PFI = as.numeric(dat$PFI)
  dat$PFI.time = as.numeric(dat$PFI.time)
  dat$PFI.time = dat$PFI.time/365
  dat$gene = factor(dat$gene, levels = c(0,1))
  gene_name = get_name(gene)
  if(is.na(gene_name)){
    gene_name = get_name_pcg(gene)
  }
  cox_mod = coxph(Surv(PFI.time, PFI) ~ gene, data = dat)
  print(glance(cox_mod)$concordance)
  conc = round(glance(cox_mod)$concordance, digits=2)
  hr = summary(cox_mod)$coefficients[2]
  dtt = dat

  fit <- survfit(Surv(PFI.time, PFI) ~ gene, data = dtt)

  s <- ggsurvplot(
          title = paste(get_name(gene), cancer, "HR =", round(hr, digits=2)),
          fit,
          xlab = "Time (Years)",
          ylab = "Progression Free Probability" ,
          #surv.median.line = "hv",
          font.main = c(4, "bold", "black"),
          font.x = c(4, "plain", "black"),
          font.y = c(4, "plain", "black"),
          font.tickslab = c(3, "plain", "black"),
          font.legend = 3,
          risk.table.fontsize = 1,
          #legend.labs = c("High Expression", "Low Expression"),             # survfit object with calculated statistics.
          data = dtt,      # data used to fit survival curves.
          risk.table = FALSE,       # show risk table.
          legend = "right",
          pval = TRUE,             # show p-value of log-rank test.
          conf.int = FALSE,        # show confidence intervals for
          xlim = c(0,10),        # present narrower X axis, but not affect
          break.time.by = 1,     # break X axis in time intervals by 500.
          palette = pal[c(2,1)])
          #risk.table.y.text.col = T, # colour risk table text annotations.
          #risk.table.y.text = FALSE )
          return(s)
}
}

#######
##[7]##-------------------------------------------------------------
#######

###EAST WAY TO FOREST PLOT (1 variable)
get_forest_plot = function(gene, cancer){
  dat = all[,c(which(colnames(all) %in% c("type", gene, "OS", "OS.time")))]
  z = which(str_detect(colnames(dat), "ENSG"))
  if(!(length(z)==0)){
  colnames(dat)[z] = "gene"
  dat = subset(dat, type == cancer)
  #split patients
  med = median(dat$gene)
  #add high low tag
    if(med ==0){
    #if median = 0 then anyone greater than zero is 1
    l1 = which(dat[,z] > 0)
    l2 = which(dat[,z] ==0)
    dat[l1,z] = 1
    dat[l2, z] = 0
    }

    if(!(med ==0)){
    l1 = which(dat[,z] >= med)
    l2 = which(dat[,z] < med)
    dat[l1,z] = 1
    dat[l2, z] = 0
    }

  dat$OS = as.numeric(dat$OS)
  dat$OS.time = as.numeric(dat$OS.time)

  #cox model using both
  cox_lnc = coxph(Surv(OS.time, OS) ~ gene, data = dat)

  #compare coefficeints
  print(ggforest(cox_lnc, data=exp_dat))
  print("done")
}
}


surv_test = function(gene){

  #print(gene)
  results_cox <- as.data.frame(matrix(ncol=6)) ; colnames(results_cox) <- c("gene", "coefficient", "HR", "pval", "low95", "upper95")

  #1. Subset lnc_rna to those patients in cancer
  z =  which(colnames(train) %in% gene)

  if(!(length(z)==0)){

  df = as.data.frame(train)
  cols_keep = c("age_at_initial_pathologic_diagnosis", "gender", "race", "clinical_stage", "histological_grade", "OS", "OS.time")
  z2 = which(colnames(train) %in% cols_keep)

  df = df[,c(z,z2)]

  gene <- colnames(df)[1]
  df$OS <- as.numeric(df$OS)
  df$OS.time <- as.numeric(df$OS.time)
  #df$stage = as.numeric(df$stage)
  #df$grade = as.numeric(df$grade)
  df$age_at_initial_pathologic_diagnosis = as.numeric(df$age_at_initial_pathologic_diagnosis)

  colnames(df)[1] = "median"

  if((!(table(df$OS)[2] <5)) & (!(dim(table(df$OS)) == 1))){

  #cox regression
  res.cox <- coxph(Surv(OS.time, OS) ~ median + age_at_initial_pathologic_diagnosis, data = df)

  row <- c(gene, summary(res.cox)$coefficients[1,c(1,2,5)],  summary(res.cox)$conf.int[1,c(3,4)])
  names(row) <- names(results_cox)
  print(row)

  return(row)

}
}
}
