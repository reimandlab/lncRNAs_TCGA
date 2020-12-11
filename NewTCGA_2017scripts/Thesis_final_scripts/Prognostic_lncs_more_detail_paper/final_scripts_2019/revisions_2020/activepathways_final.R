setwd("/Users/kisaev/Documents/lncRNAs")
library(ActivePathways)

	dat = readRDS("LGG_all_up_down_genes_fActivepathways_March25_updFC.rds")
	dat = as.data.frame(dat)
	rownames(dat) = dat$ID
	dat$ID = NULL
	dat = as.matrix(dat)
	dat_esng = dat

	colnames(dat)[which(colnames(dat)=="ENSG00000215196")] = "AC091878.1"
	colnames(dat)[which(colnames(dat)=="ENSG00000224950")] = "RP5-1086K13.1"
	colnames(dat)[which(colnames(dat)=="ENSG00000225511")] = "LINC00475.1"
	colnames(dat)[which(colnames(dat)=="ENSG00000228021")] = "RP11-383C5.31"
	colnames(dat)[which(colnames(dat)=="ENSG00000231265")] = "TRERNA1"
	colnames(dat)[which(colnames(dat)=="ENSG00000239552")] = "HOXB-AS2"
	colnames(dat)[which(colnames(dat)=="ENSG00000253187")] = "HOXA10-AS"
	colnames(dat)[which(colnames(dat)=="ENSG00000254271")] = "RP11-131N11.4"
	colnames(dat)[which(colnames(dat)=="ENSG00000254635")] = "WAC-AS1"
	colnames(dat)[which(colnames(dat)=="ENSG00000255020")] = "AF131216.5"
	colnames(dat)[which(colnames(dat)=="ENSG00000261889")] = "RP11-473M20.16"

	gmt.file <- "/Users/kisaev/Documents/lncRNAs/gprofiler_hsapiens.ENSG/gmt_file.gmt"

	gmt <- read.GMT(gmt.file)
#	gmt <- Filter(function(term) length(term$genes) >= 10, gmt)
#	gmt <- Filter(function(term) length(term$genes) <= 500, gmt)
	lgg_res = ActivePathways(dat,
		 gmt,
	cutoff = 0.05,
	significant=0.05, cytoscape.file.tag ="cytoscape", geneset.filter=c(10,500))


	gmt <- read.GMT(gmt.file)
#	gmt <- Filter(function(term) length(term$genes) >= 10, gmt)
#	gmt <- Filter(function(term) length(term$genes) <= 500, gmt)
	lgg_res = ActivePathways(dat_esng,
		 gmt,
	cutoff = 0.05,
	significant=0.05, cytoscape.file.tag ="cytoscape_esng", geneset.filter=c(10,500))
saveRDS(lgg_res, file = "LGG_all_up_down_genes_activepathways.rds")

cols = c("#72E1E0" ,"#7CE147", "#73E590", "#E180E3" ,"#D39794", "#D4DC8D", "#A6E6BB", "#A4C9E9",
 "#DC5FA5" ,"#D8DB50" ,"#6891DB" ,"#633CE1" ,"#DDDDCF" ,"#DAAADA" ,"#D849DE",
"#7C9F6A" ,"#E6565F" ,"#8C68D2" ,"#63859B" ,"#DD9955")

f=fread("cytoscapesubgroups.txt")
tests = colnames(f)[2:13]
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
