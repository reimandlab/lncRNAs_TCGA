setwd("/Users/kisaev/Documents/lncRNAs/Jan2021")
library(ActivePathways)
library(stringr)
library(data.table)

dat = readRDS("LGG_all_up_down_genes_fActivepathways_March25_updFC.rds")
dat = as.data.frame(dat)
rownames(dat) = dat$ID
dat$ID = NULL
dat = as.matrix(dat)
dat_esng = dat

colnames(dat)[which(colnames(dat)=="ENSG00000177337")] = "AC108134.2"
colnames(dat)[which(colnames(dat)=="ENSG00000215196")] = "AC026790.1"
colnames(dat)[which(colnames(dat)=="ENSG00000223768")] = "LINC00205"
colnames(dat)[which(colnames(dat)=="ENSG00000229267")] = "AC016708.1"
colnames(dat)[which(colnames(dat)=="ENSG00000239552")] = "HOXB-AS2"
colnames(dat)[which(colnames(dat)=="ENSG00000248693")] = "LINC02100"
colnames(dat)[which(colnames(dat)=="ENSG00000250237")] = "AC010273.2"
colnames(dat)[which(colnames(dat)=="ENSG00000250360")] = "AC008808.2"
colnames(dat)[which(colnames(dat)=="ENSG00000253187")] = "HOXA10-AS"
colnames(dat)[which(colnames(dat)=="ENSG00000254271")] = "AC022390.1"
colnames(dat)[which(colnames(dat)=="ENSG00000254635")] = "WAC-AS1"
colnames(dat)[which(colnames(dat)=="ENSG00000255020")] = "AF131216.3"
colnames(dat)[which(colnames(dat)=="ENSG00000261889")] = "AC108134.2"

gmt.file <- "/Users/kisaev/Documents/lncRNAs/gprofiler_hsapiens.ENSG/gmt_file.gmt"
gmt <- read.GMT(gmt.file)

#Run using gene symbols for enrichment map
lgg_res = ActivePathways(dat, merge.method  = "Brown", gmt, cytoscape.file.tag ="cytoscape")
lgg_res[which(str_detect(lgg_res$term.name, "brain")),]

#RUN on ENSG ids
gmt <- read.GMT(gmt.file)
lgg_res = ActivePathways(dat_esng, merge.method  = "Brown", gmt, cytoscape.file.tag ="cytoscape_esng")
lgg_res[which(str_detect(lgg_res$term.name, "brain")),]

saveRDS(lgg_res, file = "LGG_all_up_down_genes_activepathways.rds")

#Edit legend for gene symbol based analysis
cols = c("red", "#72E1E0", "#DD9955", "#73E590", "#E180E3" ,"#D4DC8D", "#A6E6BB", "#A4C9E9",
 "#DC5FA5" ,"#D8DB50" ,"#6891DB" ,"#633CE1" ,"#DDDDCF" ,"#DAAADA" ,"#D849DE",
"#7C9F6A" ,"#E6565F" ,"#8C68D2" ,"#63859B")

f=fread("cytoscapesubgroups.txt")
tests = colnames(f)[2:14]
col.colors <- cols[1:length(tests)]
instruct.str <- paste('piechart:',
											' attributelist="',
											paste(tests, collapse=','),
											'" colorlist="',
											paste(col.colors, collapse=','),
											'" showlabels=FALSE', sep='')
f$instruct=instruct.str
write.table(f, file="cytoscapesubgroups_new_cols.txt", quote=F, row.names=F, sep="\t")

# Making a Legend
library(ggplot2)
dummy_plot = ggplot(data.frame("tests" = factor(tests, levels = tests),
																		"value" = 1), aes(tests, fill = tests)) +
			 geom_bar() +
			 scale_fill_manual(name = "Contribution", values=col.colors)

		 grDevices::pdf(file = NULL) # Suppressing Blank Display Device from ggplot_gtable
		 dummy_table = ggplot_gtable(ggplot_build(dummy_plot))
		 grDevices::dev.off()

		 legend = dummy_table$grobs[[which(sapply(dummy_table$grobs, function(x) x$name) == "guide-box")]]

		 # Estimating height & width
		 legend_height = ifelse(length(tests) > 20,
														5.5,
														length(tests)*0.25+1)
		 legend_width = ifelse(length(tests) > 20,
													 ceiling(length(tests)/20)*(max(nchar(tests))*0.05+1),
													 max(nchar(tests))*0.05+1)
		 ggsave(legend,
						device = "pdf",
						filename = "new_legend.pdf",
						height = legend_height,
						width = legend_width,
						scale = 1)
