dim(lgg_idh)
head(lgg_idh)
lgg = as.data.table(subset(rna, type == "LGG"))
lgg = lgg[,c("ENSG00000229140", "patient", "type", "OS", "OS.time")]
lgg = merge(lgg, lgg_idh, by="patient")

med = median(lgg$ENSG00000229140)

z1 = which(lgg$ENSG00000229140 >= med)
z2 = which(lgg$ENSG00000229140 < med)

lgg$ENSG00000229140[z1] = "High"
lgg$ENSG00000229140[z2] = "Low"

colnames(lgg)[2] = "CCDC26"

#get km plots

lgg$OS.time = lgg$OS.time

gene_name = colnames(lgg)[2]

lgg$IDH.status = factor(lgg$IDH.status, levels=c("WT", "Mutant"))
lgg$OS.time = lgg$OS.time/365

fit <- survfit(Surv(OS.time, OS) ~ CCDC26 + IDH.status, data = lgg)

pdf("lgg_CCDC26_lncRNAs_cands_figure6d.pdf", width=8, height=7)

s <- ggsurvplot(
          title = gene_name, 
          fit, 
          xlab = "Time (Years)", 
          #surv.median.line = "hv",
          font.main = c(14, "bold", "black"),
          font.x = c(12, "plain", "black"),
          font.y = c(12, "plain", "black"),
          font.tickslab = c(12, "plain", "black"),
          font.legend = 10,
          risk.table.fontsize = 5, 
          #legend.labs = c("High Expression", "Low Expression"),             # survfit object with calculated statistics.
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
          #palette = mypal[c(4,1)],
          palette = "npg", 
           legend.labs = c("High exp & IDH WT", "High exp & IDH Mutant", "Low exp & IDH WT", "Low exp & IDH Mutant"), 
          #ggtheme = theme_bw(), # customize plot and risk table with a theme.
          risk.table.y.text.col = T, # colour risk table text annotations.
          risk.table.y.text = FALSE # show bars instead of names in text annotations
                            # in legend of risk table
          )

print(s)
dev.off()




head(lgg_idh)
lgg = as.data.table(subset(rna, type == "GBM"))
lgg = lgg[,c("ENSG00000229140", "patient", "type", "OS", "OS.time")]
lgg = merge(lgg, gb_idh, by="patient")

med = median(lgg$ENSG00000229140)

z1 = which(lgg$ENSG00000229140 >= med)
z2 = which(lgg$ENSG00000229140 < med)

lgg$ENSG00000229140[z1] = "High"
lgg$ENSG00000229140[z2] = "Low"

colnames(lgg)[2] = "CCDC26"

#get km plots

lgg$OS.time = lgg$OS.time

gene_name = colnames(lgg)[2]

lgg$IDH.status = factor(lgg$IDH.status, levels=c("WT", "Mutant"))
lgg$OS.time = lgg$OS.time/365

fit <- survfit(Surv(OS.time, OS) ~ CCDC26 + IDH.status, data = lgg)

pdf("gbm_CCDC26_lncRNAs_cands_figure6d.pdf", width=8, height=7)

s <- ggsurvplot(
          title = gene_name, 
          fit, 
          xlab = "Time (Years)", 
          #surv.median.line = "hv",
          font.main = c(14, "bold", "black"),
          font.x = c(12, "plain", "black"),
          font.y = c(12, "plain", "black"),
          font.tickslab = c(12, "plain", "black"),
          font.legend = 10,
          risk.table.fontsize = 5, 
          #legend.labs = c("High Expression", "Low Expression"),             # survfit object with calculated statistics.
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
          #palette = mypal[c(4,1)],
          palette = "npg", 
           #legend.labs = c("High exp & IDH WT", "High exp & IDH Mutant", "Low exp & IDH WT", "Low exp & IDH Mutant"), 
          #ggtheme = theme_bw(), # customize plot and risk table with a theme.
          risk.table.y.text.col = T, # colour risk table text annotations.
          risk.table.y.text = FALSE # show bars instead of names in text annotations
                            # in legend of risk table
          )

print(s)
dev.off()



pdf("CCD26_LGG_GBM_no_IDH.pdf")
get_km_plot("ENSG00000229140", "LGG")
get_km_plot("ENSG00000229140", "GBM")
dev.off()




