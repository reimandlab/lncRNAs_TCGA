###---------------------------------------------------------------
###Load libraries and data 
###---------------------------------------------------------------

source("source_code_Cox_MonteCarlo_CV_Mar13.R")
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

detectable = readRDS("PCAWG_detectable_genes_4cancers_March20.rds")

rna = subset(rna, canc %in% c("Kidney renal clear cell carcinoma", "Liver hepatocellular carcinoma", 
	"Ovarian serous cystadenocarcinoma", "Pancreatic adenocarcinoma"))

z = which(colnames(rna) %in% detectable$gene)
rna = as.data.frame(rna)
rna = rna[,c(z, 5786:ncol(rna))]

#how many detectable genes/cancer
cancgenes = as.data.table(table(detectable$canc))
cancgenes = cancgenes[order(N)]
pdf("num_detectacle_genesperCAncer.pdf", width=5, height=5)
ggbarplot(cancgenes, x = "V1", y = "N",
          main = "Number of detectable genes in each cancer type",
          fill = "V1",               # change fill color by cyl
          color = "white",            # Set bar border colors to white
          palette = "Dark2",            # jco journal color palett. see ?ggpar
          sort.val = "desc",          # Sort the value in dscending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 65, legend="none",
          ggtheme = theme_minimal()           # Rotate vertically x axis texts
          )
dev.off()


#remove 0 sums
sums = apply(rna[,1:(ncol(rna)-5)], 2, sum)
z = which(sums==0)
#rna = rna[,-(which(colnames(rna) == names(z)))]
groups <- as.factor(rna$canc)
res.pca <- prcomp(rna[,1:(ncol(rna)-5)],  scale = TRUE)
pdf("4cancers_PCA_using_intergenic_lncRNAs.pdf")
# Change title and axis labels
p = fviz_pca_ind(res.pca, geom="point", label="none", habillage=rna$canc, addEllipses=TRUE, ellipse.level=0.95) +
  labs(title ="PCA using expression of PCAWG detectable lncRNAs") #+ 
   #xlim(-25, 15) + ylim (-10, 20)
# Change group colors using RColorBrewer color palettes
p + scale_color_brewer(palette="Dark2") +
     theme_minimal()
dev.off()

####PCA using only cands = readRDS("chosen_features_wFANTOM_data_Mar22_1000CVs_8020splits.rds")
cands = readRDS("chosen_features_wFANTOM_data_Mar22_1000CVs_8020splits.rds")
z = which(colnames(rna) %in% cands$gene)
rna = rna[,c(z, 899:ncol(rna))]
groups <- as.factor(rna$canc)

#how many candidate lncRNAs per cancer in TCGA
candsgenes = as.data.table(table(cands$Cancer))
candsgenes = candsgenes[order(N)]
candsgenes$validated = c(1, 2, 2, 3)
colnames(candsgenes) = c("Cancer", "NumTCGAcandidates", "validated")
candsgenes$validated = as.factor(candsgenes$validated)

pdf("num_candidate_lncRNAs_perTCGAcancer.pdf", width=5, height=5)
ggbarplot(candsgenes, x = "Cancer", y = "NumTCGAcandidates",
          main = "Number of candidate genes in each cancer type",
          fill = "validated",               # change fill color by cyl
          color = "white",            # Set bar border colors to white
          palette = "Dark2",            # jco journal color palett. see ?ggpar
          sort.val = "desc",          # Sort the value in dscending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 65, legend="right",
          ggtheme = theme_minimal()           # Rotate vertically x axis texts
          )
dev.off()



res.pca <- prcomp(rna[,1:(ncol(rna)-5)],  scale = TRUE)
pdf("PCA_inTCGA_using_only_25cands.pdf")
# Change title and axis labels
p = fviz_pca_ind(res.pca, geom="point", label="none", habillage=rna$canc, addEllipses=TRUE, ellipse.level=0.95) +
  labs(title ="PCA using expression of top 25 candidate lncRNAs") #+ 
   #xlim(-25, 15) + ylim (-10, 20)
# Change group colors using RColorBrewer color palettes
p + scale_color_brewer(palette="Dark2") +
     theme_minimal()
dev.off()


####Summary of how many patients in each cancer 
pats = as.data.table(table(rna$canc))
pats = pats[order(N)]
pdf("num_patients_eachTCGAcancer.pdf", width=5, height=5)
ggbarplot(pats, x = "V1", y = "N",
          main = "Number of patients in TCGA",
          fill = "V1",               # change fill color by cyl
          color = "white",            # Set bar border colors to white
          palette = "Dark2",            # jco journal color palett. see ?ggpar
          sort.val = "desc",          # Sort the value in dscending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 65, legend="none",
          ggtheme = theme_minimal()           # Rotate vertically x axis texts
          )
dev.off()

#####-----------------------------------------------------------------------------------------------
####PCAWG-------------------------------------------------------------------------------------------
#####-----------------------------------------------------------------------------------------------

pcawg_data = readRDS("lncRNA_clinical_data_PCAWG_March20.rds")
z = which(colnames(pcawg_data) %in% detectable$gene)
pcawg_data = pcawg_data[,c(z, 6014:ncol(pcawg_data))]
pcawg_data = subset(pcawg_data, canc %in% c("Kidney Adenocarcinoma, clear cell type", "Liver Hepatocellular carcinoma", 
	"Ovary Serous cystadenocarcinoma", "Pancreas Pancreatic ductal carcinoma"))
groups <- as.factor(pcawg_data$canc)
#remove 0 sums
sums = apply(pcawg_data[,1:(ncol(pcawg_data)-5)], 2, sum)
z = which(sums==0)
res.pca <- prcomp(pcawg_data[,1:(ncol(pcawg_data)-5)],  scale = TRUE)

pdf("4cancers_PCA_using_detectable_PCAWGdata_lncRNAs.pdf")
# Change title and axis labels
p = fviz_pca_ind(res.pca, geom="point", label="none", habillage=pcawg_data$canc, addEllipses=TRUE, ellipse.level=0.95) +
  labs(title ="PCA using expression of PCAWG detectable lncRNAs") #+ 
   #xlim(-25, 15) + ylim (-10, 20)
# Change group colors using RColorBrewer color palettes
p + scale_color_brewer(palette="Dark2") +
     theme_minimal()
dev.off()


####Summary of how many patients in each cancer 
pats = as.data.table(table(pcawg_data$canc))
pats = pats[order(-N)]
pdf("num_patients_eachPCAWGcancer.pdf", width=5, height=5)
ggbarplot(pats, x = "V1", y = "N",
          main = "Number of patients in PCAWG",
          fill = "V1",               # change fill color by cyl
          color = "white",            # Set bar border colors to white
          palette = "Dark2",            # jco journal color palett. see ?ggpar
          sort.val = "desc",          # Sort the value in dscending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 65, legend="none",
          ggtheme = theme_minimal()           # Rotate vertically x axis texts
          )
dev.off()


#####-----------------------------------------------------------------------------------------------
####Overview of detectable lncRNAs------------------------------------------------------------------
#####-----------------------------------------------------------------------------------------------

z = which(fantom$gene %in% detectable$gene)
fantom = fantom[z,]
####Summary of how lncRNAs of each class
CAT_geneClass = as.data.table(table(fantom$CAT_geneClass))
CAT_geneClass = CAT_geneClass[order(-N)]
pdf("detectable_lncs_class.pdf", width=5, height=5)
ggbarplot(CAT_geneClass, x = "V1", y = "N",
          main = "Fantom gene classes for detectable lncRNAs",
          fill = "V1",               # change fill color by cyl
          color = "white",            # Set bar border colors to white
          palette = "Dark2",            # jco journal color palett. see ?ggpar
          sort.val = "desc",          # Sort the value in dscending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90, legend="none",
          ggtheme = theme_minimal()           # Rotate vertically x axis texts
          )
dev.off()
####Summary of how lncRNAs of each category
CAT_geneCategory = as.data.table(table(fantom$CAT_geneCategory))
CAT_geneCategory = CAT_geneCategory[order(-N)]
pdf("detectable_lncs_categories.pdf", width=5, height=5)
ggbarplot(CAT_geneCategory, x = "V1", y = "N",
          main = "Fantom gene categories for detectable lncRNAs",
          fill = "V1",               # change fill color by cyl
          color = "white",            # Set bar border colors to white
          palette = "Dark2",            # jco journal color palett. see ?ggpar
          sort.val = "desc",          # Sort the value in dscending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90, legend="none",
          ggtheme = theme_minimal()           # Rotate vertically x axis texts
          )
dev.off()

#both category and class
lncs = as.data.table(table(fantom$CAT_geneCategory, fantom$CAT_geneClass))
lncs = lncs[order(N)]
lncs = filter(lncs, N >0)

colnames(lncs) = c("CATgeneCategory", "CATgeneClass", "Frequency")

pdf("detectable_lncs_categories_classes.pdf", width=8, height=8)
ggbarplot(lncs, "CATgeneCategory", "Frequency",
  fill = "CATgeneClass", color = "CATgeneClass", palette = "Dark2",
  label = TRUE,
  position = position_dodge(0.9))
dev.off()

#####-----------------------------------------------------------------------------------------------
####Overview of all lncRNAs in PCAWG----------------------------------------------------------------
#####-----------------------------------------------------------------------------------------------

pcawg_data = readRDS("lncRNA_clinical_data_PCAWG_March20.rds")

fantom <- fread("lncs_wENSGids.txt", data.table=F) #6088 lncRNAs 
extract3 <- function(row){
	gene <- as.character(row[[1]])
	ens <- gsub("\\..*","",gene)
	return(ens)
}
fantom[,1] <- apply(fantom[,1:2], 1, extract3)
#remove duplicate gene names (gene names with multiple ensembl ids)
z <- which(duplicated(fantom$CAT_geneName))
rm <- fantom$CAT_geneName[z]
z <- which(fantom$CAT_geneName %in% rm)
fantom <- fantom[-z,]
colnames(fantom)[1] = "gene"

z = which(colnames(pcawg_data) %in% fantom$gene)
pcawg_data = pcawg_data[,c(z, 6014:ncol(pcawg_data))]
pcawg_data = subset(pcawg_data, canc %in% c("Kidney Adenocarcinoma, clear cell type", "Liver Hepatocellular carcinoma", 
	"Ovary Serous cystadenocarcinoma", "Pancreas Pancreatic ductal carcinoma"))
groups <- as.factor(pcawg_data$canc)

sums = apply(pcawg_data[,1:(ncol(pcawg_data)-5)], 2, sum)
z = which(sums==0)
pcawg_data = pcawg_data[,-(which(colnames(pcawg_data) == names(z)))]

z = which(fantom$gene %in% colnames(pcawg_data))
fantom = fantom[z,]

lncs = as.data.table(table(fantom$CAT_geneCategory, fantom$CAT_geneClass))
lncs = lncs[order(N)]
lncs = filter(lncs, N >0)

colnames(lncs) = c("CATgeneCategory", "CATgeneClass", "Frequency")

pdf("detectable_lncs_categories_classes_ALL_pcawg_lncs.pdf", width=8, height=8)
ggbarplot(lncs, "CATgeneCategory", "Frequency",
  fill = "CATgeneClass", color = "CATgeneClass", palette = "Dark2",
  label = TRUE,
  position = position_dodge(0.9))
dev.off()

#####-----------------------------------------------------------------------------------------------
####Plot distribution plot comparing to PCGs--------------------------------------------------------
#####-----------------------------------------------------------------------------------------------

lnc_rna <- readRDS("6028_pcawg_lncRNAs_RNASeq_data.rds")
lnc_rna <- as.data.frame(lnc_rna)
lnc_rna$patient <- rownames(lnc_rna)

#PCGs
pcg_rna <- readRDS("20166_pcawg_PCGs_RNASeq_data.rds")
pcg_rna <- as.data.frame(pcg_rna)
pcg_rna$patient <- rownames(pcg_rna)

#remove duplicated column names 
dups <- colnames(pcg_rna)[which(duplicated(colnames(pcg_rna)))]   
#save them in a list for future reference 
#pcg_rna <- pcg_rna[,-(which(colnames(pcg_rna) %in% dups))]

#Clinical file - available only for 485/497 patients 
clin <- readRDS("Jan26_PCAWG_clinical")
z <- which(clin$icgc_donor_id %in% rownames(lnc_rna))
clin <- clin[z,]

lnc_rna <- lnc_rna[which(rownames(lnc_rna) %in% clin$icgc_donor_id),] #485 patients remain
pcg_rna <- pcg_rna[which(rownames(pcg_rna) %in% clin$icgc_donor_id),] #485 patients remain 

#only look at lncRNAs included in fantom
z = which(colnames(lnc_rna) %in% fantom$gene)
lnc_rna = lnc_rna[,z]
#which are detectable in all cancers?
#meds = apply(lnc_rna, 2, median)
#det = meds[meds>=1]

#For each patient add survival status and days since last seen 
lnc_rna$canc = ""
lnc_rna$status = ""
lnc_rna$time = ""
lnc_rna$sex = ""

#lncs
for(i in 1:nrow(lnc_rna)){
  pat <- rownames(lnc_rna)[i]
  z <- which(clin$icgc_donor_id %in% pat)
  lnc_rna$canc[i] <- clin$histology_abbreviation[z]
  lnc_rna$status[i] <- clin$donor_vital_status[z]
  lnc_rna$sex[i] <- clin$donor_sex[z]
  t <- clin$donor_survival_time[z]
  if(is.na(t)){
        t <- clin$donor_interval_of_last_followup[z]
        }
        lnc_rna$time[i] <- t
}

#1. remove all genes that are not expressed at all
lnc_rna = subset(lnc_rna, canc %in% c("Kidney-RCC", "Liver-HCC", "Ovary-AdenoCA", "Panc-AdenoCA"))
pcg_rna = pcg_rna[which(rownames(pcg_rna) %in% rownames(lnc_rna)),]

sums = apply(lnc_rna[,1:(ncol(lnc_rna)-4)], 2, sum)
z = which(sums==0)

sums = apply(pcg_rna[,1:(ncol(pcg_rna)-2)], 2, sum)
z = which(sums==0)
pcg_rna = pcg_rna[,-(which(colnames(pcg_rna) == names(z)))]


####################################################
#2. get median FPKM for each gene 
#####PLOTTING#######################################

pdf("dist_lncRNA_medians_PCAWG_ind_cancers.pdf", width=5, height=5)
for(i in 1:length(unique(lnc_rna$canc))){
canc = subset(lnc_rna, canc == unique(lnc_rna$canc)[i])
meds_lncs = apply(canc[,1:(ncol(canc)-4)], 2, median)
meds_lncs = as.data.frame(meds_lncs)
meds_lncs$gene = rownames(meds_lncs)
meds_lncs$type = ""
for(i in 1:nrow(meds_lncs)){
	z = which(fantom$gene == meds_lncs$gene[i])
	meds_lncs$type[i] = fantom$CAT_geneCategory[z]
}
colnames(meds_lncs)[1] = "median"
meds_lncs = as.data.table(meds_lncs)

pcg_canc = pcg_rna[which(rownames(pcg_rna) %in% rownames(canc)),]
meds_pcgs = apply(pcg_canc[,1:(ncol(pcg_canc)-2)], 2, median)
meds_pcgs = as.data.frame(meds_pcgs)
meds_pcgs$gene = rownames(meds_pcgs)
meds_pcgs$type = "pcg"
colnames(meds_pcgs)[1] = "median"
meds_pcgs = as.data.table(meds_pcgs)

all_meds = rbind(meds_lncs, meds_pcgs)
#all_meds$median = floor(all_meds$median)
all_meds$type = as.factor(all_meds$type)

#remove super highly expressed outliers
z = which(all_meds$median >= 1000)
if(!(length(z)==0)){
all_meds = all_meds[-z,]
}

#zoomed in 
all_meds$gene = NULL
gg <- ggplot(all_meds) + labs(title=canc$canc[1], y = "density", x="median FPKM")
gg <- gg + geom_density(aes(x=median, y=..scaled.., colour=type, fill=type), alpha=0.01)
gg <- gg + theme_bw()
print(ggpar(gg, xlim=c(0,100)))

#zoomed out
gg <- ggplot(all_meds) + labs(title=canc$canc[1], y = "density", x="median FPKM")
gg <- gg + geom_density(aes(x=median, y=..scaled.., colour=type, fill=type), alpha=0.01)
gg <- gg + theme_bw()
#ggpar(p, xlim=c(0,250))
print(gg)

#3. log medians first 
all_meds$median = log1p(all_meds$median)
#zoomed in 
gg <- ggplot(all_meds) + labs(title=canc$canc[1], y = "density", x="log1p(median FPKM)")
gg <- gg + geom_density(aes(x=median, y=..scaled.., colour=type, fill=type), alpha=0.01)
gg <- gg + theme_bw()
print(gg)
#print(ggpar(gg, xlim=c(0,100)))
}
dev.off()


#####-----------------------------------------------------------------------------------------------
####Looking at only detectable lncRNAs--------------------------------------------------------------
#####-----------------------------------------------------------------------------------------------

#detectable only lncs
lnc_rna = lnc_rna[,c(which(colnames(lnc_rna) %in% detectable$gene), 6013:ncol(lnc_rna))]
#PLOTTING - get medians 
pdf("dist_lncRNA_medians_PCAWG_ind_cancers_ONLY_detectable_ones.pdf", width=5, height=5)
for(i in 1:length(unique(lnc_rna$canc))){
canc = subset(lnc_rna, canc == unique(lnc_rna$canc)[i])
meds_lncs = apply(canc[,1:(ncol(canc)-4)], 2, median)
meds_lncs = as.data.frame(meds_lncs)
meds_lncs$gene = rownames(meds_lncs)
meds_lncs$type = ""
for(i in 1:nrow(meds_lncs)){
  z = which(fantom$gene == meds_lncs$gene[i])
  meds_lncs$type[i] = fantom$CAT_geneCategory[z]
}
colnames(meds_lncs)[1] = "median"
meds_lncs = as.data.table(meds_lncs)

pcg_canc = pcg_rna[which(rownames(pcg_rna) %in% rownames(canc)),]
meds_pcgs = apply(pcg_canc[,1:(ncol(pcg_canc)-2)], 2, median)
meds_pcgs = as.data.frame(meds_pcgs)
meds_pcgs$gene = rownames(meds_pcgs)
meds_pcgs$type = "pcg"
colnames(meds_pcgs)[1] = "median"
meds_pcgs = as.data.table(meds_pcgs)

all_meds = rbind(meds_lncs, meds_pcgs)
#all_meds$median = floor(all_meds$median)
all_meds$type = as.factor(all_meds$type)

#remove super highly expressed outliers
z = which(all_meds$median >= 1000)
if(!(length(z)==0)){
all_meds = all_meds[-z,]
}

#zoomed in 
all_meds$gene = NULL
gg <- ggplot(all_meds) + labs(title=canc$canc[1], y = "density", x="median FPKM")
gg <- gg + geom_density(aes(x=median, y=..scaled.., colour=type, fill=type), alpha=0.01)
gg <- gg + theme_bw()
print(ggpar(gg, xlim=c(0,100)))

#zoomed out
gg <- ggplot(all_meds) + labs(title=canc$canc[1], y = "density", x="median FPKM")
gg <- gg + geom_density(aes(x=median, y=..scaled.., colour=type, fill=type), alpha=0.01)
gg <- gg + theme_bw()
#ggpar(p, xlim=c(0,250))
print(gg)

#3. log medians first 
all_meds$median = log1p(all_meds$median)
#zoomed in 
gg <- ggplot(all_meds) + labs(title=canc$canc[1], y = "density", x="log1p(median FPKM)")
gg <- gg + geom_density(aes(x=median, y=..scaled.., colour=type, fill=type), alpha=0.01)
gg <- gg + theme_bw()
print(gg)
#print(ggpar(gg, xlim=c(0,100)))
}
dev.off()






































