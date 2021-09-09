lgg = as.data.table(filter(all, type == "LGG"))
lgg$OS.time = lgg$OS.time/365

fit <- survfit(Surv(OS.time, OS) ~ 1, data = lgg)
# Drawing curves
#ggsurvplot(fit, color = "#2E9FDF")

pdf("lgg_gbm_km_plots.pdf", width=8, height=7)

s <- ggsurvplot(
          title = "LGG Overall Survival", 
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
          # legend.labs = c("High exp & IDH WT", "High exp & IDH Mutant", "Low exp & IDH WT", "Low exp & IDH Mutant"), 
          #ggtheme = theme_bw(), # customize plot and risk table with a theme.
          risk.table.y.text.col = T, # colour risk table text annotations.
          risk.table.y.text = FALSE # show bars instead of names in text annotations
                            # in legend of risk table
          )
print(s)

gbm = as.data.table(filter(all, type == "GBM"))
gbm$OS.time = gbm$OS.time/365

fit <- survfit(Surv(OS.time, OS) ~ 1, data = gbm)
# Drawing curves
s <- ggsurvplot(
          title = "GBM Overall Survival", 
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
          # legend.labs = c("High exp & IDH WT", "High exp & IDH Mutant", "Low exp & IDH WT", "Low exp & IDH Mutant"), 
          #ggtheme = theme_bw(), # customize plot and risk table with a theme.
          risk.table.y.text.col = T, # colour risk table text annotations.
          risk.table.y.text = FALSE # show bars instead of names in text annotations
                            # in legend of risk table
          )
print(s)
dev.off()





brain = as.data.table(filter(all, type %in% c("LGG", "GBM")))
brain$OS.time = brain$OS.time/365

fit <- survfit(Surv(OS.time, OS) ~ type, data = brain)
s <- ggsurvplot(
          title = "Gliomas Overall Survival", 
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
          data = brain,      # data used to fit survival curves. 
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
          # legend.labs = c("High exp & IDH WT", "High exp & IDH Mutant", "Low exp & IDH WT", "Low exp & IDH Mutant"), 
          #ggtheme = theme_bw(), # customize plot and risk table with a theme.
          risk.table.y.text.col = T, # colour risk table text annotations.
          risk.table.y.text = FALSE # show bars instead of names in text annotations
                            # in legend of risk table
          )

pdf("all_gliomas_km_plots.pdf", width=8, height=7)
print(s)
dev.off()










