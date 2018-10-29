#check GBM and LGG IDH status influence on lncRNA expression and clinical outcomes 

lgg = as.data.table(filter(clean_up, canc == "Brain Lower Grade Glioma"))
gbm = as.data.table(filter(clean_up, canc == "Glioblastoma multiforme"))

#get IDH variables
z = unique(lgg$colname[which(str_detect(lgg$colname, 'IDH'))])
col = "IDH.status"

lgg = as.data.table(filter(lgg, colname == col))

#hox genes 
 g=subset(allCands, cancer== "Glioblastoma multiforme")
 g=g[1,1]

#plot their expression in gbm 
library(TCGAbiolinks)
gbm_subtype <- TCGAquery_subtype(tumor = "lgg")

pdf("lncRNA_cands_fromGBM_inLGG.pdf")
for(i in 1:length(g$gene)){
    gene=g[i]
    gene=unlist(gene)
    exp = which(colnames(all) %in% c(gene, "type", "patient", "OS", "OS.time", "age_at_initial_pathologic_diagnosis.y", "gender.y"))
    exp = all[,exp]
    exp = subset(exp, type=="LGG")
    #which have IDH mutation 
    z = which(exp$patient %in% gbm_subtype$patient)
    exp = exp[z,]
    z = which(str_detect(colnames(gbm_subtype), "IDH.status"))
    z1=which(colnames(gbm_subtype) %in% "patient")
    gbm_subtype = gbm_subtype[,c(z,z1)]
    exp=merge(exp, gbm_subtype, by = "patient")
    med = median(exp[,5])
    z=which(exp[,5] >= med)
    exp$tag = ""
    exp$tag[z] = "High"
    exp$tag[-z] = "Low"
    colnames(exp)[5] = "gene_exp"
    exp[,5] = log1p(exp[,5]) 
    exp$combo = paste(exp$tag, exp$IDH.status, sep="_")
    exp$combo = factor(exp$combo, levels = c("Low_Mutant", "Low_WT", "High_Mutant", "High_WT"))
    z = which(is.na(exp$IDH.status))
    if(!(length(z)==0)){
        exp = exp[-z,]
    }

    p <- ggboxplot(exp, x = "combo", y = "gene_exp",
          color = "combo", title = paste(get_name(gene), "Expression in LGG"),
          add = "jitter", ylab = "lncRNA exp log1p(FPKM-UQ)",  ggtheme = theme_bw()) +
          stat_compare_means() + geom_hline(yintercept=log1p(med), linetype="dashed", color = "red") + 
          stat_n_text(size=5)

        p = ggpar(p,
          font.xtickslab = c(14,"plain", "black"),font.tickslab=c(14,"plain", "black"), 
          xtickslab.rt = 55, legend="none")
        print(p)
    
    exp$OS = as.numeric(exp$OS)
    exp$OS.time = as.numeric(exp$OS.time)
    exp$OS.time = exp$OS.time/365
    exp$tag = factor(exp$tag, levels = c("Low","High"))
    
    #lnc model
    lnc_model = coxph(Surv(OS.time, OS) ~ tag, data = exp)
    hr_lnc = round(summary(lnc_model)$coefficients[2], digits=3)

    #idh model
    idh_model = coxph(Surv(OS.time, OS) ~ IDH.status, data = exp)
    hr_idh = round(summary(idh_model)$coefficients[2], digits=3)

    #lnc + idh model
    both = coxph(Surv(OS.time, OS) ~ tag + IDH.status, data = exp)

    #lnc vs both
    lnc_vs_both = anova(lnc_model, both)[2,4]

    #idh vs both
    idh_vs_both = anova(idh_model, both)[2,4]


    fit <- survfit(Surv(OS.time, OS) ~ tag + IDH.status, data = exp)
          s <- ggsurvplot(
          fit, 
          title=paste(get_name(gene), "Expression in LGG"), 
          xlab = "Time (Years)", 
          surv.median.line = "hv",
          font.main = c(16, "bold", "black"),
          font.x = c(14, "plain", "black"),
          font.y = c(14, "plain", "black"),
          font.tickslab = c(14, "plain", "black"),
          font.legend = 8,
          risk.table.fontsize = 5, 
          #legend.labs = c("Low Expression", "High Expression"),             # survfit object with calculated statistics.
          data = exp,      # data used to fit survival curves. 
          risk.table = TRUE,       # show risk table.
          legend = "right", 
          pval = TRUE,             # show p-value of log-rank test.
          conf.int = FALSE,        # show confidence intervals for 
                            # point estimaes of survival curves.
          xlim = c(0,10),        # present narrower X axis, but not affect
                            # survival estimates.
          break.time.by = 1,     # break X axis in time intervals by 500.
          #palette = colorRampPalette(mypal)(14), 
          palette = mypal[c(4,1,2,3)],
          ggtheme = theme_bw(), # customize plot and risk table with a theme.
          risk.table.y.text.col = T, # colour risk table text annotations.
          risk.table.y.text = FALSE # show bars instead of names in text annotations
                            # in legend of risk table
          )
          print(s)
}
dev.off()


#all 9 lncRNAs are associated with IDH status and the model with both predictors
#is stronger than either one alone 

#add name 
lgg$name = unlist(llply(lgg$lnc, get_name))

#get IDH variables - GBM

z = unique(gbm$colname[which(str_detect(gbm$colname, 'IDH'))])
col = "IDH.specific.DNA.Methylation.Cluster"

gbm = as.data.table(filter(gbm, colname == col))

dtt = clin_data_lncs[[13]]

#idh vs lncRNAs 

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

    	i = which(colnames(new_dat) == "IDH.status")

      	print(i)    
      	col = colnames(new_dat)[i]
    

     	new_dat_plot = new_dat[,c("patient", col, lnc, "lncRNA_tag", "risk")]
       
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

  
        #oncoprint 
        new_dat_plot = as.data.table(new_dat_plot)
        new_dat_plot = new_dat_plot[order(Clinical, lncRNA_tag)]

        exp = new_dat_plot[,c("patient", "lncRNA_tag")]
        exp$type = "lncRNA_exp"
        mut = new_dat_plot[,c("patient", "Clinical")]
        mut$type = "IDH_status"
        mut$Clinical = as.character(mut$Clinical)
        mut$Clinical[mut$Clinical == "WT"] = 0
        mut$Clinical[mut$Clinical == "Mutant"] = 1
        colnames(mut)[2] = "lncRNA_tag"
        all = rbind(exp, mut)
        all$patient = factor(all$patient, levels=new_dat_plot$patient)
        all$lncRNA_tag = factor(all$lncRNA_tag, levels=c(1,0))
        name = get_name(lnc)
		g = ggplot(all, aes(patient, type)) +
		  geom_tile(aes(fill = lncRNA_tag), colour = "grey50") +
		  theme(axis.text.x = element_blank(), axis.ticks = element_blank())+
    		scale_fill_manual(values = c("red", "blue"), name="Value",
                       breaks=c(1, 0),
                       labels=c("HighExp (top) or IDH mut (bottom)", "LowExp (top) or IDH WT (bottom)"))+     		
    		xlab("patient") + ylab("Variable") +
 		ggtitle(paste(name, "IDH status", "\nRisk = ", new_dat_plot$risk[1]))+theme(legend.position="top")
		print(g)
          }
          } #check >1 
        }
      }
    }

    pdf("oncoplot.pdf", height=5, width=8)
    llply(lncs, get_cor)
    dev.off()






