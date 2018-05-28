source("source_code_Cox_MonteCarlo_CV_April12.R")
require(caTools)

#start with only lncRNA_intergenic
#lincs = subset(fantom, CAT_geneClass == "lncRNA_intergenic")
#z = which(colnames(rna) %in% lincs$gene)
#rna = as.data.frame(rna)
#rna = rna[,c(z, (ncol(rna)-5):ncol(rna))]

###[2.] Data splitting 

###---------------------------------------------------------------
###PCA using lncRNA expression 
#can then compare how using all genes compared to just using
#the ones chosen by LASSO at the end 
###---------------------------------------------------------------

#check if this person is in my analysis: TCGA-61-2095
library(glmnet)
library(survcomp)
library(caret)
library(stringr)
library(factoextra)
library(plyr)
library(doParallel)

rna = as.data.frame(rna)
norm = as.data.frame(norm)


#only look at cancers with at least 50 patients in them 
cancers_dist = filter(as.data.table(table(rna$Cancer)), N >=50)
rna = subset(rna, rna$Cancer %in% cancers_dist$V1)
norm = subset(norm, norm$Cancer %in% rna$Cancer)

ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
#z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
#ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,6]))
ucsc <- ucsc[-z,]

#add location to lncRNAs 

#how many lncRNAs  --> 5,785

rownames(rna) = rna$patient
z1 = which(str_detect(colnames(rna), "ENSG"))	
z2 = which(colnames(rna) %in% "Cancer")
rna = rna[,c(z1,z2)]

#normal lncRNAs
z = which(colnames(norm) %in% colnames(rna))
norm = norm[,z]
#make sure looking at same lncRNAs 
z = which(colnames(rna) %in% colnames(norm))
rna = rna[,z]
#only keep normal tissues with at least 10 patient samples
pats = as.data.table(table(norm$Cancer))
pats = pats[order(N)]
pats = filter(pats, N >=10)
norm = subset(norm, norm$Cancer %in% pats$V1)

#subset rna file to just cancers in normal file
rna = subset(rna, rna$Cancer %in% norm$Cancer)

#how many detectable lncRNAs in normal tissues versus tumours 

cancers = unique(rna$Cancer)

get_canc_data = function(cancer){
	norm_data = subset(norm, Cancer == cancer)
	norm_data$type = "normal" 
	canc_data = subset(rna, Cancer == cancer)
	canc_data$type = "cancer"
	all_data = rbind(norm_data, canc_data)
	return(all_data)
}

canc_datas = llply(cancers, get_canc_data)

#get median value for each lncRNA
colnames(ucsc)[6] = "gene"
ucsc = ucsc[,c(6, 2, 4,5)]
colnames(ucsc) = c("gene", "chr", "start", "end")

#need to figure out % of people that have expression greater than 1 FPKM 
#in the cohort 

get_medians = function(dtt){
	
	norm = subset(dtt, type == "normal")
	canc = subset(dtt, type=="cancer")

	z1 = which(str_detect(colnames(dtt), "ENSG"))	
	#dtt[,z1] = log1p(dtt[,z1])
	#meds 
	lncs = colnames(dtt)[1:5781]
		calc_freq_canc = function(lnc){
		print(lnc)
		newdat = canc[,which(colnames(canc) %in% lnc)]
		madd = mad(newdat)
		#l = (length(which(newdat >=80)))/nrow(dtt)
		if(madd > 0){stat = "detectable"}
		if(madd <= 0){stat = "NOTdetectable"}
		return(stat)
		}


		calc_freq_norm = function(lnc){
		print(lnc)
		newdat = norm[,which(colnames(norm) %in% lnc)]
		madd = mad(newdat)
		#l = (length(which(newdat >=80)))/nrow(dtt)
		if(madd > 0){stat = "detectable"}
		if(madd <= 0){stat = "NOTdetectable"}
		return(stat)
		}


	summary_dat_canc = as.data.frame(matrix(ncol = 3, nrow=length(lncs))) ; colnames(summary_dat_canc) = c("cancer", "lncRNA", "status")
	summary_dat_canc$cancer = canc$Cancer[1]
	summary_dat_canc$lncRNA = lncs
	summary_dat_canc$status = unlist(llply(lncs, calc_freq_canc))
	summary_dat_canc$num_patient = nrow(canc)
	summary_dat_canc$type = "cancer"

	summary_dat_norm = as.data.frame(matrix(ncol = 3, nrow=length(lncs))) ; colnames(summary_dat_norm) = c("cancer", "lncRNA", "status")
	summary_dat_norm$cancer = norm$Cancer[1]
	summary_dat_norm$lncRNA = lncs
	summary_dat_norm$status = unlist(llply(lncs, calc_freq_norm))
	summary_dat_norm$num_patient = nrow(norm)
	summary_dat_norm$type = "normal"

	all_results = rbind(summary_dat_canc, summary_dat_norm)
	return(all_results)

}

meds_cancers = llply(canc_datas, get_medians, .progress="text")
meds_cancers1 = ldply(meds_cancers, data.frame)

#justcancerand num patients
pats = meds_cancers1[,c(1, 4:5)]
pats = pats[!duplicated(pats), ]
colnames(pats)[1] = "Cancer"
colnames(pats)[3] ="tissuetype"

all_detectable = subset(meds_cancers1, status=="detectable")

lncspercancer = as.data.table(table(meds_cancers1$cancer, meds_cancers1$type, meds_cancers1$status))
lncspercancer = filter(lncspercancer, V3 == "detectable")
colnames(lncspercancer) = c("Cancer", "tissuetype", "lncRNAStatus", "NumDetected")
lncspercancer = merge(lncspercancer, pats, by= c("Cancer", "tissuetype"))

lncspercancer = as.data.table(lncspercancer)
lncspercancer = lncspercancer[order(NumDetected)]
lncspercancer$NumDetected = lncspercancer$NumDetected/5781
#assign levels 
order = unique(lncspercancer$Cancer)

lncspercancer$Cancer <- factor(lncspercancer$Cancer, levels = lncspercancer$Cancer[order(lncspercancer$NumDetected)])
lncspercancer$Cancer  # notice the changed order of factor levels

pdf("Figure1A_lncspercancer_norms.pdf", height=5, width=7)

g = ggbarplot(lncspercancer, x="Cancer", y="NumDetected", color = "grey", fill = "tissuetype", palette = "Paired",
  position = position_dodge(0.9)) + theme_light() 

ggpar(g,
 font.tickslab = c(7,"plain", "black"),
 xtickslab.rt = 65)

dev.off()


######look at visual of overall lncRNA expression in cancer compared to normal

#median expression tumours
meds_lncs_canc = as.data.frame(apply(rna[,1:5781], 2, mad))
colnames(meds_lncs_canc)[1] = "median"
meds_lncs_canc$type = "cancer"

#median expression normals 
meds_lncs_norm = as.data.frame(apply(norm[,1:5781], 2, mad))
colnames(meds_lncs_norm)[1] = "median"
meds_lncs_norm$type = "normal"

all_meds = rbind(meds_lncs_norm, meds_lncs_canc)
all_meds$median = as.numeric(all_meds$median)
all_meds$median = log1p(all_meds$median)

summary(all_meds$median)
pdf("mad_dis_tum_vs_normal.pdf")
ggviolin(all_meds, x = "type", y = "median", color="type", palette = mypal, xlab="MAD") +  theme_light() +  stat_compare_means()    
dev.off()


########calculate number of lncRNAs up/downregulated for each cancer type 

#for each cancer types, get number of lncRNAs "differentailly expressed"

get_diff_exp = function(dtt){
	library(stringr)
	z1 = which(str_detect(colnames(dtt), "ENSG"))	
	lncs = as.list(colnames(dtt)[z1]) #5781 
	#for each lncRNA calculate mean difference in expression
	#plot expression
	pdf(file=paste(dtt$Cancer[1], "diff_expressed_lncRNAs.pdf", sep="_"))

	each_lnc = function(lnc, dtt){
		#subset data to lncRNA
		z = which(colnames(dtt) %in% c(lnc, "type"))
		justlnc= dtt[,z]
		justlnc[,1] = log1p(justlnc[,1])
		colnames(justlnc)[1] ="lncRNA"
		w = wilcox.test(justlnc$lncRNA[justlnc$type=="cancer"], justlnc$lncRNA[justlnc$type=="normal"])$p.value
		median_diff = median(justlnc$lncRNA[justlnc$type=="cancer"]) / median(justlnc$lncRNA[justlnc$type=="normal"])
		mean_diff = mean(justlnc$lncRNA[justlnc$type=="cancer"]) / mean(justlnc$lncRNA[justlnc$type=="normal"])
		library(ggpubr)
		#print(ggdensity(justlnc, x="lncRNA", color="type", title=paste(dtt$Cancer[1], lnc, "wilcoxon-p=", round(w, digits=2))))
		row = c(dtt$Cancer[1], lnc, w, median_diff, mean_diff)
		names(row) = c("Cancer", "lnc", "Wilcoxon_p", "median_diff", "mean_diff")
		return(row)
	}
	library(parallel)
	cl <- makeCluster(5)
	library(doParallel)
	registerDoParallel(cl)
	opts <- list(preschedule=TRUE)
	clusterSetRNGStream(cl, 123)
	#lncs = lncs[1:100]
	lncs_summary = llply(lncs, dtt, 
           .fun = each_lnc,
           .parallel = TRUE,
           .paropts = list(.options.snow=opts))
	dev.off()
	#combine
	library(data.table)
	lncs_summary1 = do.call("rbind", lncs_summary)
	lncs_summary1 = as.data.frame(lncs_summary1)
	lncs_summary1$fdr = ""
	lncs_summary1$fdr = p.adjust(as.numeric(lncs_summary1$Wilcoxon_p), method="fdr")
	lncs_summary1 = as.data.table(lncs_summary1)
	lncs_summary1 = lncs_summary1[order(fdr)]
	return(lncs_summary1)
}

#cl <- makeCluster(5)
#registerDoParallel(cl)
#opts <- list(preschedule=TRUE)
#clusterSetRNGStream(cl, 123)

#canc_datas = canc_datas[1:2]

lncs_summary_all = llply(canc_datas,
           .fun = get_diff_exp,
          	.progress = "text")

lncs_summary_all = ldply(lncs_summary_all, data.frame)

saveRDS(lncs_summary_all, "cancer_normals_lncRNA_diff_expression_analysis_May28.rds")


#keep only fdr significant 
lncs_summary_all = as.data.table(lncs_summary_all)
lncs_summary_all = as.data.table(filter(lncs_summary_all, fdr <=0.05))
lncs_summary_all = lncs_summary_all[order(mean_diff)]
z = which(lncs_summary_all$median_diff == "Inf")
lncs_summary_all = lncs_summary_all[-z,]
z = which(lncs_summary_all$mean_diff == "Inf")
lncs_summary_all = lncs_summary_all[-z,]
lncs_summary_all$fdr = as.numeric(lncs_summary_all$fdr)
lncs_summary_all$mean_diff = as.numeric(lncs_summary_all$mean_diff)
lncs_summary_all$mean_diff = log2(lncs_summary_all$mean_diff)


pdf("summary_tumour_normal_diff_expression_May28th.pdf", width=9, height=7)
g = ggviolin(lncs_summary_all, x="Cancer", y="mean_diff", ylab="log2(Difference in mean expression)", fill="fdr", add = "mean_sd") + theme_light() + 
geom_jitter(aes(colour = fdr), position=position_jitter(width=0.05,height=0),
         alpha=0.7,
         size=0.5) + scale_color_gradient(low = "black", high = "red")

ggpar(g,
 font.tickslab = c(7,"plain", "black"),
 xtickslab.rt = 65) + geom_hline(yintercept=c(-(log2(2/3)), log2(0.5)), linetype="dashed", color = "red") 

dev.off()

#summarize how many significant "diff expressed" lncRNAs for each cancer type

lncs_summary_all$diff = ""
lncs_summary_all$diff[lncs_summary_all$mean_diff > log2(1.5)] = "Upregulated"
lncs_summary_all$diff[lncs_summary_all$mean_diff <= log2(2/3)] = "Downregulated"
lncs_summary_all$diff[lncs_summary_all$diff == ""] = "NotDiff"

summary_diff = as.data.table(table(lncs_summary_all$Cancer, lncs_summary_all$diff))
summary_diff = subset(summary_diff, V2 %in% c("Downregulated", "Upregulated"))
colnames(summary_diff) = c("Cancer", "TypeRegulated", "NumberofLNCRNAS")


pdf("summary_barplot_tumour_normal_diff_expression_May28th.pdf", width=9, height=7)
ggbarplot(summary_diff, "Cancer", "NumberofLNCRNAS", fill = "TypeRegulated", color = "grey", palette = "Paired",
  label = TRUE,
  position = position_dodge(0.9))
dev.off()

saveRDS(summary_diff, "summary_detectable_cancer_normals_lncRNA_diff_expression_analysis_May28.rds")





