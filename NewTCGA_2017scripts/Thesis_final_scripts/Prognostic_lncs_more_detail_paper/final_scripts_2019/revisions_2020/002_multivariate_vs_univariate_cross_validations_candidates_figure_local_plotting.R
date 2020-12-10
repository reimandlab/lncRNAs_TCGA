library(ggpubr)
library(data.table)
library(dplyr)

setwd("/Users/kisaev/Documents/lncRNAs")
input_plot = readRDS("multivaraite_vs_univariate_cindices_input_for_plots.rds")

input_plot$type = factor(input_plot$type, levels=c("lncRNA", "lncRNA&clin",
	"all lncRNAs", "all lncRNAs + clinical", "clinical"))

#get order of cancers by decreasing
lncs = as.data.table(filter(input_plot, type == "lncRNA"))
lncs = as.data.table(lncs %>% group_by(type_canc) %>% summarize(mean=median(cindex)))
lncs = lncs[order(-mean)]
input_plot$type_canc = factor(input_plot$type_canc, levels=unique(lncs$type_canc))

pdf("multivaraite_vs_univariate_cindices_input_for_plots.pdf", width=11, height=5)

 g = ggerrorplot(input_plot, x = "type", y = "cindex",
            desc_stat = "median_mad", palette=c("darkred", "darksalmon", "dodgerblue4", "skyblue2", "grey"),
            error.plot = "pointrange",            # Change error plot type
            add = "median" , color="type" ,
             add.params = list(size = 0.0001, shape=19)              # Add mean points
            ) + theme_bw()

g = facet(g, facet.by = "type_canc", panel.labs.font = list(size =4), nrow=1)

g=ggpar(g, legend="none", font.xtickslab=c(6, "plain", "black"), font.ytickslab=c(4, "plain", "black")) +
ylim(c(0.4,1))+
	geom_hline(yintercept=0.5, linetype="dashed", color = "red") +
	  geom_text(aes(label = lncRNA), size=1, y=1)+rremove("x.text") +rremove("x.ticks")+xlab("")+ylab("c-index")
	  print(g)


dev.off()


pdf("multivaraite_vs_univariate_cindices_input_for_plots_barplot.pdf", width=8, height=5)

g = ggbarplot(input_plot, x = "type", y = "cindex",
            palette=c("darkred", "darksalmon", "dodgerblue4", "skyblue2", "grey"),
            add = c("median_mad"),            # Change error plot type
            fill="type", facet.by="type_canc") + theme_bw()

#g = facet(g, facet.by = "type_canc", panel.labs.font = list(size =4), nrow=1)
#g=ggpar(g, font.xtickslab=c(6, "plain", "black"), font.ytickslab=c(4, "plain", "black")) +
#ylim(c(0.4,1))+
#	geom_hline(yintercept=0.5, linetype="dashed", color = "red") +
#	  rremove("x.text") +xlab("")+ylab("c-index")
#	  print(g)

print(g)


dev.off()
