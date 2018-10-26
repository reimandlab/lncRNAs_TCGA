
#for GBM 
#add extra patient samples frome xternal data set
ext = readRDS("all_genes_external_tcga_all_cancers_March13_wclinical_data.rds")

#canc_conv
canc_conv = unique(rna[,c("type", "Cancer")])

#check if cands are significant using data from ext 
pats = as.data.table(table(ext$type))
pats = as.data.table(filter(pats, N >= 15))
colnames(pats)[1] = "type"

pats = merge(pats, canc_conv, by="type")

###EASY WAY TO MAKE KM PLOT
get_km_plot = function(gene, cancer){
  all_g = ext
  all_g = as.data.frame(all_g)
  dat = ext[,c(which(colnames(ext) %in% c("type", gene, "OS", "OS.time")))]
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
  
  #balance check
  bal_check = (table(dat[,1])[1] >= 5) & (table(dat[,1])[2] >= 5)
  if(bal_check){

  cox_mod = coxph(Surv(OS.time, OS) ~ gene, data = dat)
  print(glance(cox_mod)$concordance)
  conc = round(glance(cox_mod)$concordance, digits=2)
  fit <- survfit(Surv(OS.time, OS) ~ gene, data = dat)
          s <- ggsurvplot(
          title = paste(gene_name, dat$type[1], "\nConcordance=", conc),
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
          data = dat,      # data used to fit survival curves. 
          risk.table = TRUE,       # show risk table.
          legend = "right", 
          pval = TRUE,             # show p-value of log-rank test.
          conf.int = FALSE,        # show confidence intervals for 
                            # point estimaes of survival curves.
          xlim = c(0,5),        # present narrower X axis, but not affect
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
}
}	
}


for(i in 1:length(unique(allCands$combo))){
	
	z = which(pats$Cancer %in% allCands$cancer[i])
	if(!(length(z)==0)){
	canc = pats$type[z]
	gene = allCands$gene[i]

	print(get_km_plot(gene, canc))
}
}

dev.off()




#get gbm
gbm = subset(ext, type=="GBM")

z = which(colnames(all) %in% colnames(gbm))
cols = colnames(all)[z]
all = all[,z]

z = which(colnames(gbm) %in% colnames(all))
cols = colnames(gbm)[z]
gbm = gbm[,z]

r = rbind(all, gbm)

all=r


pdf("four_GBM_LGG_ionchannels_of_interest.pdf")

get_km_plot(get_ensg_pcg("CATSPER1"), "GBM")
get_km_plot(get_ensg_pcg("CATSPER1"), "LGG")

get_km_plot(get_ensg_pcg("SCN9A"), "GBM")
get_km_plot(get_ensg_pcg("SCN9A"), "LGG")

get_km_plot(get_ensg_pcg("AQP9"), "GBM")
get_km_plot(get_ensg_pcg("AQP9"), "LGG")

get_km_plot(get_ensg_pcg("KCNN4"), "GBM")
get_km_plot(get_ensg_pcg("KCNN4"), "LGG")

get_km_plot(get_ensg_pcg("GJB2"), "GBM")
get_km_plot(get_ensg_pcg("GJB2"), "LGG")

dev.off()