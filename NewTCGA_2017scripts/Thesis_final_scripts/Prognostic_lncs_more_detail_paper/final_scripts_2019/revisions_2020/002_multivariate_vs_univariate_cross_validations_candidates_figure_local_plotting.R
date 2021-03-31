library(ggpubr)
library(data.table)
library(dplyr)

setwd("/Users/kisaev/Documents/lncRNAs/Jan2021")

#color palette
colours_palette=readRDS("23_cancers_color_palette.rds")

myColors=colours_palette$color
names(myColors)=colours_palette$cancer
myColors
colScale <- scale_colour_manual(name = "type",values = myColors)

#get colours for 30 cancers
row=c("PRAD", "gray26")
colours_palette=rbind(colours_palette, row)
row=c("TGCT", "darkblue")
colours_palette=rbind(colours_palette, row)
row=c("THYM", "yellow")
colours_palette=rbind(colours_palette, row)
row=c("READ", "darkgreen")
colours_palette=rbind(colours_palette, row)
row=c("KICH", "sienna1")
colours_palette=rbind(colours_palette, row)
row=c("PCPG", "hotpink")
colours_palette=rbind(colours_palette, row)
row=c("UVM", "maroon1")
colours_palette=rbind(colours_palette, row)

myColors_full=colours_palette$color
names(myColors_full)=colours_palette$cancer
myColors_full
colScale_full <- scale_colour_manual(name = "cancer",values = myColors_full)

input_plot = readRDS("multivaraite_vs_univariate_cindices_input_for_plots.rds")

input_plot$type = factor(input_plot$type, levels=c("lncRNA", "all lncRNAs","lncRNA&clin",
	"all lncRNAs + clinical", "clinical"))

#get order of cancers by decreasing
lncs = as.data.table(filter(input_plot, type == "lncRNA"))
lncs = as.data.table(lncs %>% group_by(type_canc) %>% summarize(mean=median(cindex)))
lncs = lncs[order(-mean)]
input_plot$type_canc = factor(input_plot$type_canc, levels=unique(lncs$type_canc))

#one plot color for each cancer
g = ggline(input_plot, x = "type", y = "cindex",
					 add = c("median") , color="type_canc" ,
						add.params = list(size = 0.0001, shape=19)) + theme_classic()
g=ggpar(g, legend="bottom", font.legend=c(3, "plain", "black"),
font.xtickslab=c(6, "plain", "black"), font.ytickslab=c(8, "plain", "black")) +
ylim(c(0.4,1))+ colScale+
 geom_hline(yintercept=0.5, linetype="dashed", color = "red") +
	  xlab("Model Type")+ylab("c-index")
#print(g)
ggsave("multivaraite_vs_univariate_cindices_one_plot_median_shown.pdf",g, width=4, height=6, units="in")

g = ggline(input_plot, x = "type", y = "cindex",
					 add = c("median_mad") , color="type_canc" ,
						add.params = list(size = 0.0001, shape=19)) + theme_classic()
g=ggpar(g, legend="bottom", font.legend=c(3, "plain", "black"),
font.xtickslab=c(6, "plain", "black"), font.ytickslab=c(8, "plain", "black")) +
ylim(c(0.4,1))+ colScale+
 geom_hline(yintercept=0.5, linetype="dashed", color = "red") +
	  xlab("Model Type")+ylab("c-index")
#print(g)
ggsave("multivaraite_vs_univariate_cindices_one_plot_median_mad_shown.pdf",g, width=4, height=6, units="in")


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


pdf("multivaraite_vs_univariate_cindices_input_for_plots_barplot.pdf", width=8, height=7)

g = ggbarplot(input_plot, x = "type", y = "cindex",
            palette=c("darkred", "darksalmon", "dodgerblue4", "skyblue2", "grey"),
            add = c("median_mad"),            # Change error plot type
            fill="type", facet.by="type_canc") + theme_bw()+rremove("x.text")+ylab("c-index")+
						theme(strip.text.x = element_text(size = 4, colour = "black"))+xlab("Model type")+
						theme(strip.background = element_rect(colour="black", fill="white",
                                       size=0.5, linetype="solid"))

#g = facet(g, facet.by = "type_canc", panel.labs.font = list(size =4), nrow=1)
#g=ggpar(g, font.xtickslab=c(6, "plain", "black"), font.ytickslab=c(4, "plain", "black")) +
#ylim(c(0.4,1))+
#	geom_hline(yintercept=0.5, linetype="dashed", color = "red") +
#	  rremove("x.text") +xlab("")+ylab("c-index")
#	  print(g)

print(g)
dev.off()
