source("/u/kisaev/lncRNAs_TCGA/NewTCGA_2017scripts/Thesis_final_scripts/Prognostic_lncs_more_detail_paper/final_scripts_2019/revisions_2020/load_data.R")

#------FEATURES-----------------------------------------------------
setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNAs_2019_manuscript")

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = subset(allCands, data == "TCGA") #173 unique lncRNA-cancer combos, #166 unique lncRNAs
allCands$combo = unique(paste(allCands$gene, allCands$cancer, sep="_"))

mypal = c("#EAD286" ,"#D1EB7B", "#96897F" ,"#E5C0A6" ,
  "#72A93B", "#74DAE3" ,"#49B98D" ,"#D97B8F" ,"#70A2A4", "#64709B" ,"#DFBF38" ,"#61EA4F" ,
  "#C7CBE7", "#786DDA",
"#CFA0E0" ,"#67E9D0" ,"#7C9BE1", "#D94753" ,
"#AAE6B0", "#D13BDF" ,"#DEAEC7" ,"#BBE6DF" ,"#B2B47A" ,"#E6ECBA", "#C86ED7",
 "#7BEE95" ,"#6F46E6" ,"#65B9E0", "#C0EC3E",
"#DE8D54" ,"#DF4FA6")

mypal = "purple"

#library(M3C)
library(umap)

lgg=readRDS("TCGA_lgg_wsubtype_info_biolinks.rds")
lgg_wsurv = unique(lgg[,c("patient", "IDH.status", "OS", "OS.time")])
lgg = unique(lgg[,c("patient", "IDH.status")])
lgg_rna = as.data.table(filter(rna, type =="LGG"))
lgg_pfi = unique(lgg_rna[,c("patient", "PFI", "PFI.time")])
lgg_pfi = merge(lgg, lgg_pfi, by="patient")
saveRDS(lgg_pfi, "/u/kisaev/TCGA_LGG_PFS_time.rds")

#-------------------------------------------------------------------
#-----------------PCA using just all expressed lncRNA --------------
#-------------------------------------------------------------------

dim(rna)
rownames(rna) = rna$patient

z1 = which(str_detect(colnames(rna), "ENSG"))
z2 = which(colnames(rna) %in% c("type", "patient", "OS", "OS.time"))
rna = as.data.frame(rna)
rna = rna[,c(z1, z2)]

get_umap_lgg = function(canc){

	canc_dat = as.data.table(filter(rna, type==canc))

	if(canc == "LGG"){
		canc_dat = as.data.table(merge(lgg, canc_dat, by="patient"))
    z = which(is.na(canc_dat$IDH.status))
    canc_dat = canc_dat[-z,]
		cancer.labels = canc_dat$IDH.status
		canc_dat$IDH.status = NULL
	}

	#1. remove those not expressed at all
	z1 = which(str_detect(colnames(canc_dat), "ENSG"))
	canc_dat[,z1] = log1p(canc_dat[,..z1])

	set.seed(100)

	if(!(canc %in% c("LGG", "KIRC"))){
		cancer.labels = canc_dat$type
	}
	df = canc_dat
	df$type = NULL
	df$patient = NULL
  df$OS = NULL
  df$OS.time = NULL
	embedding = umap::umap(df)
#	head(embedding)

	layout = as.data.table(embedding$layout) ; colnames(layout)=c("x", "y")
	layout$cols_names = cancer.labels

	z = which(duplicated(layout$col))
	layout$label[z] ="no"
	layout$label[-z] = "yes"

	z = (which(layout$label == "yes"))
	for(i in 1:length(z)){
		canc = layout$col[z[i]]
		layout$label[z[i]] = canc
	}


  g = ggplot(layout, aes(x=x, y=y, colour=cols_names)) +
    geom_point() +
    geom_point(data=layout %>%
                 group_by(cols_names) %>%
                 summarise_at(c("x", "y"), mean),
               size=5, shape=3) +
    theme_classic()+
    scale_colour_manual(values=c("purple", "orange"))

  #geom_point(aes(x=x, y=y, color=cols_names),alpha=0.5, stroke=0) +
  #	geom_text_repel(data = subset(layout, !(label == "no")))+
  #theme_classic()
	print(g)

}

pdf("/u/kisaev/Jan2021/LGG_IDH_umap_plots.pdf")
get_umap_lgg("LGG")

#make survival plot just for IDH mutation status

z = which(is.na(lgg_wsurv$IDH.status))
lgg_wsurv = lgg_wsurv[-z,]
lgg_wsurv$OS = as.numeric(lgg_wsurv$OS)
lgg_wsurv$OS.time = as.numeric(lgg_wsurv$OS.time)
lgg_wsurv$OS.time = lgg_wsurv$OS.time/365

fit <- survfit(Surv(OS.time, OS) ~ IDH.status, data = lgg_wsurv)

s <- ggsurvplot(fit ,
        xlab = "Time (Years)",
        data = lgg_wsurv,      # data used to fit survival curves.
        pval = TRUE,             # show p-value of log-rank test.
        conf.int = FALSE,        # show confidence intervals for
        xlim = c(0,10),
        risk.table = TRUE,      # present narrower X axis, but not affect
        break.time.by = 1,     # break X axis in time intervals by 500.
        palette = c("purple", "orange"))
print(s)
dev.off()
