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

rna = as.data.frame(rna)
norm = as.data.frame(norm)

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

all_detectable = subset(meds_cancers1, status="detectable")

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
ggdensity(all_meds, x = "median", color="type", palette = mypal, xlab="MAD") +  theme_light()
dev.off()


########calculate number of lncRNAs up/downregulated for each cancer type 







































