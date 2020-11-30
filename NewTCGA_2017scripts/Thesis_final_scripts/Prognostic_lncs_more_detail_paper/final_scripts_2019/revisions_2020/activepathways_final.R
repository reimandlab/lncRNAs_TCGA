	library(ActivePathways)

	dat = readRDS("LGG_all_up_down_genes_fActivepathways_March25_updFC.rds")

	lgg_res = ActivePathways(dat,
	"all_paths.txt", cutoff = 0.05,
	significant=0.01, cytoscape.file.dir="cytoscape", geneset.filter=c(10,500))


	saveRDS(lgg_res, file = "LGG_all_up_down_genes_activepathways.rds")
