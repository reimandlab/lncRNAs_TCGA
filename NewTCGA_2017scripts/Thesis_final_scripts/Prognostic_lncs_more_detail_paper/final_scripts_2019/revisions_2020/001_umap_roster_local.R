library(ggpubr)
library(data.table)
library(dplyr)
library(ggrastr)
library(ggrepel)

setwd("/Users/kisaev/Documents/lncRNAs/Jan2021")

#colours
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

#umap data
umap_dat = fread("umap_29_cancers_data_for_figure.txt")

#plot
g = ggplot(umap_dat,aes(x, y, label = label)) +
geom_point_rast(aes(x=x, y=y, colour=cancer),alpha=0.4, stroke=0) +
#scale_colour_manual(values=mypal)+
colScale_full+theme_classic()+
geom_text_repel(data = subset(umap_dat, !(label == "no")))

ggsave("UMAP_30_cancers_ggrastr.pdf",g, width=5, height=3, units="in")
