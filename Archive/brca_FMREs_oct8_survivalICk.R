#get breast cancer data

dtt = clin_data_lncs[[6]]
z = which(colnames(dtt) %in% c("patient", "Cancer", "PAM50.mRNA"))
dtt= dtt[,z]

#idh vs lncRNAs 

canc = dtt$Cancer[1]
#Combined into one dataframe because need to get ranks 
com = colnames(pcg)[which(colnames(pcg) %in% colnames(rna))]
all <- merge(rna, pcg, by = com)  

mypal = readRDS("best_pal.rds")

get_ensg_pcg = function(pcg){
  z = which(ucsc$hg19.ensemblToGeneName.value == pcg)
  if(length(z)>1){
    z = z[1]
  }
  return(ucsc$hg19.ensGene.name2[z])
}

get_subtype_km_plots = function(gene){
    
    gene_id = get_ensg_pcg(gene)
    z = which(colnames(all) %in% c("patient", "type", "Cancer", gene_id, "OS", "OS.time"))

    dat = all[,z]
    dat = merge(dat, dtt, by=c("patient", "Cancer"))

    #first is survival overall associated with ICK expression
    #when subtype is accounted for
    colnames(dat)[6] = "gene"

    dat$OS = as.numeric(dat$OS)
    dat$OS.time=as.numeric(dat$OS.time)
    dat$PAM50.mRNA = as.character(dat$PAM50.mRNA)
    z = which(is.na(dat$PAM50.mRNA))
    dat=dat[-z,]

    med = median(dat$gene)
    if(!(med == 0)){
    z1 = which(dat$gene < med)
    z2 = which(dat$gene >= med)
    dat$gene_tag = ""
    dat$gene_tag[z1] = "low"
    dat$gene_tag[z2] = "high"}

    if(med==0){
    z1 = which(dat$gene == 0)
    z2 = which(dat$gene > 0)
    dat$gene_tag = ""
    dat$gene_tag[z1] = "low"
    dat$gene_tag[z2] = "high" 
    }

    surv_gene= coxph(Surv(OS.time, OS) ~ gene_tag, data = dat)
    surv_combo = coxph(Surv(OS.time, OS) ~ gene_tag + PAM50.mRNA, data = dat)

    out <- split(dat , f = dat$PAM50.mRNA)
    print(gene)

    get_km = function(dato){
        if(dim(dato)[1] >20){
            med = median(dato$gene)
            
            if(!(med == 0)){
            z1 = which(dato$gene < med)
            z2 = which(dato$gene >= med)
            dato$gene_tag = ""
            dato$gene_tag[z1] = "low"
            dato$gene_tag[z2] = "high"}

            if(med==0){
            z1 = which(dato$gene == 0)
            z2 = which(dato$gene > 0)
            dato$gene_tag = ""
            dato$gene_tag[z1] = "low"
            dato$gene_tag[z2] = "high" 
            }
            
            print(dato$PAM50.mRNA[1])
            dato$gene_tag = factor(dato$gene_tag, levels=c("low", "high"))

            surv_gene= coxph(Surv(OS.time, OS) ~ gene_tag, data = dato)
            #plot KM 
            dato$OS.time = dato$OS.time/365
            fit <- survfit(Surv(OS.time, OS) ~ gene_tag, data = dato)
              s <- ggsurvplot(
              title = paste(gene, dato$PAM50.mRNA[1], "\nHR =" , round(summary(surv_gene)$coefficients[2], digits=3),
              "pval =" , round(summary(surv_gene)$coefficients[5], digits=3)), 
              fit, 
              xlab = "Time (Years)", 
              surv.median.line = "hv",
              font.main = c(16, "bold", "black"),
              font.x = c(14, "plain", "black"),
              font.y = c(14, "plain", "black"),
              font.tickslab = c(14, "plain", "black"),
              font.legend = 12,
              risk.table.fontsize = 5, 
              #legend.labs = c("low", "high"),             # survfit object with calculated statistics.
              data = dato,      # data used to fit survival curves. 
              risk.table = TRUE,       # show risk table.
              legend = "right", 
              #pval = TRUE,             # show p-value of log-rank test.
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
              t = (as.data.table(glance(surv_gene)))
              t$gene = gene
              t$HR = summary(surv_gene)$coefficients[2]
              t$pval = summary(surv_gene)$coefficients[5]
              return(t)
        }    
    }

    genes_results = llply(out, get_km)
    genes_results = ldply(genes_results)
    print(paste("done", gene))
    return(genes_results)
}

genes = c("MDM4", "ICK", "ZKSCAN3", "NBPF1", "CLPTM1L", "AC104794.4")

pdf("brca_subtypes_KM_plots.pdf")

surv_models_genes = llply(genes, get_subtype_km_plots, .progress="text")

dev.off()


surv_models_genes1 = ldply(surv_models_genes)
surv_models_genes1 = as.data.table(surv_models_genes1)
surv_models_genes1 = surv_models_genes1[order(p.value.log)]
surv_models_genes1$pval = NULL
write.csv(surv_models_genes1, file="BRCA_diff_exp_genes_subtype_survival_results.csv", quote=F, row.names=F)






