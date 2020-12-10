setwd("/Users/kisaev/Documents/lncRNAs")
library(ActivePathways)

	dat = readRDS("LGG_all_up_down_genes_fActivepathways_March25_updFC.rds")
	dat = as.data.frame(dat)
	rownames(dat) = dat$ID
	dat$ID = NULL
	dat = as.matrix(dat)

	gmt.file <- "/Users/kisaev/Documents/lncRNAs/gprofiler_hsapiens.ENSG/gmt_file.gmt"

	gmt <- read.GMT(gmt.file)
#	gmt <- Filter(function(term) length(term$genes) >= 10, gmt)
#	gmt <- Filter(function(term) length(term$genes) <= 500, gmt)

	lgg_res = ActivePathways(dat,
		 gmt,
	cutoff = 0.05,
	significant=0.05, cytoscape.file.tag ="cytoscape", geneset.filter=c(10,500))

saveRDS(lgg_res, file = "LGG_all_up_down_genes_activepathways.rds")
