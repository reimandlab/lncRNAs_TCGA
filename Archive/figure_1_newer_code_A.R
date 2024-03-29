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

#------FEATURES-----------------------------------------------------

#check if this person is in my analysis: TCGA-61-2095
library(glmnet)
library(survcomp)
library(caret)
library(stringr)
library(factoextra)

rna = as.data.frame(rna)

ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
#z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
#ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,6]))
ucsc <- ucsc[-z,]

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = filter(allCands, data=="TCGA") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 

cands_dups = unique(allCands$gene[which(duplicated(allCands$gene))])
cands = allCands

#add location to lncRNAs 

#how many lncRNAs  --> 5,785

rownames(rna) = rna$patient

z1 = which(str_detect(colnames(rna), "ENSG"))	
z2 = which(colnames(rna) %in% c("type", "OS", "OS.time"))
rna = rna[,c(z1,z2)]

#need to get median expression for each gene for each cancer type 
cancers = unique(rna$type)

#remove cancer types 
#with less than 50 patients
pats_num = as.data.table(table(rna$type))
pats_num = filter(pats_num, N <50)
canc_rm = pats_num$V1

#remove those ones
cancers = cancers[which(!(cancers %in% canc_rm))]

get_canc_data = function(cancer){
	canc_data = subset(rna, type == cancer)
	return(canc_data)
}

canc_datas = llply(cancers, get_canc_data)

#get median value for each lncRNA
colnames(ucsc)[6] = "gene"
ucsc = ucsc[,c(6, 2, 4,5)]
colnames(ucsc) = c("gene", "chr", "start", "end")

#need to figure out % of people that have expression greater than 1 FPKM 
#in the cohort 

get_medians = function(dtt){
	z1 = which(str_detect(colnames(dtt), "ENSG"))	
	#dtt[,z1] = log1p(dtt[,z1])
	#meds 
	lncs = colnames(dtt)[1:5785]
	calc_freq = function(lnc){
		print(lnc)
		canc_data_genes_analyze = dtt
		z = which(colnames(canc_data_genes_analyze) == lnc)
  		dat = canc_data_genes_analyze[,c(z,(ncol(canc_data_genes_analyze)-2):ncol(canc_data_genes_analyze))]
  		dat$OS.time = as.numeric(dat$OS.time)
  		dat$OS = as.numeric(dat$OS)
  		#remove NAs
  		z = which(is.na(dat$OS.time))
  		if(!(length(z) ==0)){
  		  dat = dat[-z,]}
		med_gene = median(as.numeric(dat[,1]))  	
		dat$med = ""
		if(med_gene ==0){
        #if median = 0 then anyone greater than zero is 1 
        l1 = which(dat[,1] > 0)
        l2 = which(dat[,1] ==0)
        dat$med[l1] = 1
        dat$med[l2] = 0
        }

    	if(!(med_gene ==0)){
        l1 = which(dat[,1] >= med_gene)
        l2 = which(dat[,1] < med_gene)
        dat$med[l1] = 1
        dat$med[l2] = 0
    	}

    if(dim(table(dat$med)) >1){
    check1 = table(dat$med)[1] >= 15
    check2 = table(dat$med)[2] >= 15
    med = med_gene
    if(check1 & check2){
	stat = "detectable"}
	if(!(check1 & check2)){
	stat = "NOTdetectable"}
	}

	if(dim(table(dat$med)) <= 1){
    med = med_gene
	stat = "NOTdetectable"}
	
	med = round(med, digits=2)
	nonzero = length(which(dat[,1] > 100))
	return(paste(stat, med, nonzero, sep="_"))
	}

	summary_dat = as.data.frame(matrix(ncol = 3, nrow=length(lncs))) ; colnames(summary_dat) = c("cancer", "lncRNA", "status")
	summary_dat$cancer = dtt$type[1]
	summary_dat$lncRNA = lncs
	summary_dat$status = unlist(llply(lncs, calc_freq))
	summary_dat$num_patient = nrow(dtt)

	return(summary_dat)

}

meds_cancers = llply(canc_datas, get_medians, .progress="text")
meds_cancers1 = ldply(meds_cancers, data.frame)

meds_cancers1$stat = sapply(meds_cancers1$status, function(x){unlist(strsplit(x, "_"))[1]})
meds_cancers1$pats_det = sapply(meds_cancers1$status, function(x){unlist(strsplit(x, "_"))[3]})

z = which(meds_cancers1$stat == "detectable")
all_detectable = meds_cancers1[z,]
saveRDS(all_detectable, file="all_TCGA_cancers_lncRNAs_detectable_July20.rds")

z = which((meds_cancers1$stat == "detectable") & (meds_cancers1$pats_det >= 15))
meds_cancers1$real_detz[z] = "ReallyDet"
meds_cancers1$real_detz[-z] = "NotDetected"
lncspercancer = as.data.table(table(meds_cancers1$cancer, meds_cancers1$real_detz))
lncspercancer = filter(lncspercancer, V2 == "ReallyDet")
colnames(lncspercancer) = c("Cancer", "lncRNAStatus", "NumDetected")

lncspercancer = as.data.table(lncspercancer)
lncspercancer = lncspercancer[order(NumDetected)]
lncspercancer$NumDetected = lncspercancer$NumDetected/5785
#assign levels 
order = lncspercancer$Cancer

lncspercancer$Cancer <- factor(lncspercancer$Cancer, levels = lncspercancer$Cancer[order(lncspercancer$NumDetected)])
lncspercancer$Cancer  # notice the changed order of factor levels

pdf("Figure1A_new_ratio_split_july20.pdf", height=5, width=7)

g = ggbarplot(lncspercancer, x="Cancer", y="NumDetected", fill = "steelblue")

ggpar(g,
 font.tickslab = c(7,"plain", "black"),
 xtickslab.rt = 45)

dev.off()


#now how many cancers per lncRNA 
meds_cancers1 = as.data.table(meds_cancers1)
det_lncs = as.data.table(filter(meds_cancers1, real_detz=="ReallyDet"))
det_lncs = as.data.table(table(det_lncs$lncRNA, det_lncs$cancer))
det_lncs =filter(det_lncs, N >=1)
det_lncs = table(det_lncs$V1)
det_lncs = as.data.table(det_lncs)
det_lncs = det_lncs[order(N)]

#414 not detectable in any cancer --> are they detectable in normal tissues 

highlncs = fantom$gene[which(fantom$CAT_geneName %in% c("MALAT1", "NEAT1", "HOTAIR", "XIST"))]
colnames(det_lncs)[1] = "gene"
det_lncs  = as.data.frame(det_lncs)
det_lncs$gene = as.factor(det_lncs$gene)
labels = det_lncs$gene

det_lncs$N = as.numeric(det_lncs$N)
det_lncs$group = ""
det_lncs$group[det_lncs$N %in% c(1:4)] = "1:4"
det_lncs$group[det_lncs$N %in% c(5:9)] = "5:9"
det_lncs$group[det_lncs$N %in% c(10:14)] = "10:14"
det_lncs$group[det_lncs$N %in% c(15:19)] = "15:19"
det_lncs$group[det_lncs$N %in% c(20:24)] = "20:24"
det_lncs$group[det_lncs$N %in% c(25:28)] = "25:28"

cancersperlncs = as.data.table(table(det_lncs$group))
order = c("1:4",
"5:9",
"10:14",
"15:19",
"20:24",
"25:28")

cancersperlncs$N = cancersperlncs$N/5785

pdf("Figure1A_part2_new_ratio_split_july20.pdf", height=5, width=7)

g = ggbarplot(cancersperlncs, x = "V1", y = "N",
            xlab = "Number of detectable cancers",
            ylab = "Number of lncRNAs", order=order, 
            fill = "steelblue")
#ggpar(g, xticks.by =2)
g
dev.off()

head(fantom)
det_lncs = merge(det_lncs, fantom, by="gene")
det_lncs = as.data.table(det_lncs)
all_cancers = as.data.table(filter(det_lncs, N ==32))

#within those expressed in all cancer types 
#40% of them, what kinds of lcnRNAs are they? 
all_cancers_types = as.data.table(table(all_cancers$CAT_geneCategory))
all_cancers_types = as.data.table(all_cancers_types[order(N)])

all_cancers_types$V1 <- factor(all_cancers_types$V1, levels = all_cancers_types$V1[order(all_cancers_types$N)])
all_cancers_types$V1  # notice the changed order of factor levels
all_cancers_types$N = all_cancers_types$N/nrow(all_cancers)

pdf("Figure1A_part3_FPKM5.pdf", height=5, width=7)

g = ggbarplot(all_cancers_types, x="V1", y="N", fill = "steelblue", color = "steelblue", 
	lab.pos = "in", lab.col = "white", lab.size=2.2, #label = TRUE, 
	 xlab = "lncRNA Classification",
            ylab = "Proportion of lncRNAs detected in all cancers")

ggpar(g,
 font.tickslab = c(10,"plain", "black"),
 xtickslab.rt = 45)

dev.off()


#is there a signficinat enrichment of divergent? 
expected = table(det_lncs$CAT_geneCategory)/nrow(det_lncs)
all_cancers_types = as.data.table(table(all_cancers$CAT_geneCategory))
all_cancers_types = as.data.table(all_cancers_types[order(N)])
observed = all_cancers_types$N 
names(observed) = all_cancers_types$V1
expected = expected[c(4, 1, 2, 3)]


observed = as.data.frame(observed)
expected = table(det_lncs$CAT_geneCategory)/nrow(det_lncs)*sum(observed$observed)
expected = expected[c(4, 1, 2, 3)]

observed$expected = expected
observed$type = rownames(observed)
observed = melt(observed)

observed$value = round(observed$value, digits=1)
pdf("expected_vs_observed_fractions_lncRNA_types.pdf", width=9, height=5)
ggbarplot(observed, x="type", y="value", fill="variable", palette = "Paired",
  label = TRUE,
  position = position_dodge(0.9)) + theme_light() 
dev.off()


#2. What is the overlap of detectable lncRNAs? 
#use only lncRNAs detectable in at least 1 cancer 
dets = as.character(unique(det_lncs$gene))
rna = rna[,which(colnames(rna) %in% c(dets, "Cancer"))]

#rows need to be cancer types and columns are lncRNAs 
#so need to get a summary metric for each lncRNA
#try mean and median 

######CLUSTERING##################################################################################################

hier = as.data.frame(matrix(ncol=ncol(rna)-1, nrow=length(unique(rna$Cancer))))
rownames(hier) = unique(rna$Cancer)
colnames(hier) = colnames(rna)[1:(ncol(rna)-1)]

#inefficient but whatever
for(i in 1:nrow(hier)){
	canc = rownames(hier)[i]
	canc_data = subset(rna, Cancer == canc)
	canc_data = canc_data[,1:(ncol(canc_data)-1)]
	for(y in 1:ncol(hier)){
		lnc = colnames(hier)[y]
		lnc_med = median(canc_data[,which(colnames(canc_data) %in% lnc)])
		hier[i,y] = lnc_med
	}
}

hier = scale(hier)
df = hier
# Dissimilarity matrix
d <- dist(df, method = "euclidean")
# Hierarchical clustering using Complete Linkage
hc1 <- hclust(d, method = "ward" )
# Plot the obtained dendrogram
plot(hc1, cex = 0.6, hang = -1)
# Compute with agnes
hc2 <- agnes(df, method = "ward")
# Agglomerative coefficient
hc2$ac
## [1] 0.8531583
# methods to assess
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")
# function to compute coefficient
ac <- function(x) {
  agnes(df, method = x)$ac
}
map_dbl(m, ac)
# average    single  complete      ward 
#0.4083263 0.3567615 0.4963445 0.5877997 
#i think this one is the best one 
hc3 <- agnes(df, method = "ward")
pltree(hc3, cex = 0.6, hang = -1, main = "Dendrogram of agnes") 
# compute divisive hierarchical clustering
hc4 <- diana(df)
# Divise coefficient; amount of clustering structure found
hc4$dc
## [1] 0.8514345
# plot dendrogram
pltree(hc4, cex = 0.6, hang = -1, main = "Dendrogram of diana")
# Ward's method
hc5 <- hclust(d, method = "ward.D2" )
# Cut tree into 4 groups
sub_grp <- cutree(hc5, k = 4)
pdf("Figure1B_FPKM5.pdf", width=7)
plot(hc5, cex = 0.6)
#rect.hclust(hc5, k = 4)
dev.off()
# Number of members in each cluster
table(sub_grp)
#fviz_cluster(list(data = df, cluster = sub_grp))

######EXPRESSION######################################################################################################

#show that even though lncRNAs are detectable in all cancer types, their expression still varies between cancer types 
highlncs

highlncs = rna[,which(colnames(rna) %in% c(highlncs, "Cancer"))]
highlncs[,1:4] = log1p(highlncs[,1:4])

pdf("Figure1B_part3.pdf", height=5, width=8)

for(i in 1:4){
	
	df = highlncs[,c(i, 5)]
	meds = as.data.table(aggregate(df[,1], list(df$Cancer), median))		
	meds = meds[order(x)]
	order= meds$Group.1

	df$Cancer <- factor(df$Cancer, levels = order)
	df$Cancer  # notice the changed order of factor levels

	lnc = colnames(df)[1]
	colnames(df)[1] = "lncRNA"

	g = ggboxplot(df, x = "Cancer", y = "lncRNA", fill = "grey", color = "steelblue", title=paste(lnc, "Expression"))
	g = ggpar(g,
 		font.tickslab = c(7,"plain", "black"),
 		xtickslab.rt = 60, legend = "none") + 

	    stat_compare_means(method = "anova", label.y = 20, label.x=5)       # Add global annova p-value
  		#stat_compare_means(label = "p.signif", method = "t.test",
        #             ref.group = ".all.", hide.ns = TRUE)      # Pairwise comparison against all
	
	print(g)

}


dev.off()


######PCGS###########################################################################################################

#get median expression of protein coding genes 
#get median expression of lncRNA expression 

#x-axis median values 
#y-axis frequency 

#get all lncRNA median expression across all cancers 
meds_lncs = apply(rna[,1:(ncol(rna)-1)], 2, median)
meds_lncs = as.data.frame(meds_lncs)
meds_lncs$type = "lncRNA"
colnames(meds_lncs)[1] = "medianFPKM"
meds_lncs$gene = rownames(meds_lncs)
meds_lncs = merge(meds_lncs, fantom, by="gene")
meds_lncs = meds_lncs[,c(1,2,3,5)]
colnames(meds_lncs)[4] = "genetype"

#get all pcg median expression across all cancers
rownames(pcg) = pcg$patient
meds_pcgs = apply(pcg[,2:19351], 2, median)
meds_pcgs = as.data.frame(meds_pcgs)
meds_pcgs$gene = rownames(meds_pcgs)
meds_pcgs$type = "pcg"
colnames(meds_pcgs)[1] = "medianFPKM"
meds_pcgs$genetype = "pcg"
meds_pcgs = meds_pcgs[,c("gene", "medianFPKM", "type", "genetype")]

all_meds = rbind(meds_lncs, meds_pcgs)
rownames(all_meds) = 1:nrow(all_meds)

all_meds$medianFPKM = as.numeric(all_meds$medianFPKM)
all_meds$medianFPKM = log1p(all_meds$medianFPKM)
all_meds = as.data.frame(all_meds)

pdf("figure1supllement.pdf", width=8, height=5)
ggdensity(all_meds, x = "medianFPKM", color="genetype", palette = mypal)
dev.off()


#how much bigger is the pcg median 
meds = as.data.table(aggregate(all_meds[,2], list(all_meds$type), median))	

pdf("figure1supllementB.pdf", width=8, height=5)
ggboxplot(all_meds, x = "type", y = "medianFPKM", color="type", palette = mypal)+
stat_compare_means(method = "wilcox.test") 
dev.off()


saveRDS(det_lncs, file="summary_detectable_lncs_howmanycancers_typesLNCRNAS.rds")













