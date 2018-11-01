library(survAUC)
#source("source_code_Cox_MonteCarlo_CV_Mar13.R")
require(caTools)
#check if this person is in my analysis: TCGA-61-2095
library(glmnet)
library(survcomp)
library(caret)
library(stringr)

source("universal_LASSO_survival_script.R")

library(ggpubr)
library(ggrepel)
library(viridis)
library(patchwork)
library(caret)  
library(Rtsne)
library("FactoMineR")


#TO-DO:
#clutser all 6,000 lncRNAs in each cancer type
#then look if those selected by elastic net meaningfully represnent 
#clusters 


#------FEATURES-----------------------------------------------------

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = filter(allCands, data=="TCGA") #175 unique lncRNA-cancer combos, #166 unique lncRNAs 

canc_conv = rna[,which(colnames(rna) %in% c("Cancer", "type"))]
canc_conv = canc_conv[!duplicated(canc_conv), ]

mypal = c("#E5DFD9","#EAD286" ,"#D1EB7B", "#96897F" ,"#E5C0A6" ,
  "#72A93B", "#74DAE3" ,"#49B98D" ,"#D97B8F" ,"#70A2A4", "#64709B" ,"#DFBF38" ,"#61EA4F" ,
  "#C7CBE7", "#786DDA",
"#CFA0E0" ,"#67E9D0" ,"#7C9BE1", "#D94753" ,
"#AAE6B0", "#D13BDF" ,"#DEAEC7" ,"#BBE6DF" ,"#B2B47A" ,"#E6ECBA", "#C86ED7",
 "#7BEE95" ,"#6F46E6" ,"#65B9E0", "#C0EC3E",
"#DE8D54" ,"#DF4FA6")

#-------------------------------------------------------------------
#-----------------PCA using just all expressed lncRNA --------------
#-------------------------------------------------------------------

#remove cancer types with less than 50 patients 
pats_num = as.data.table(table(rna$Cancer))
pats_num = filter(pats_num, N <50)
canc_rm = pats_num$V1

#remove those ones
cancers = unlist(cancers[which(!(cancers %in% canc_rm))])
rna = subset(rna, Cancer %in% cancers)

dim(rna)
rownames(rna) = rna$patient
z1 = which(str_detect(colnames(rna), "ENSG"))
z2 = which(colnames(rna) %in% "type")
rna = rna[,c(z1, z2)]

#1. remove those not expressed at all
z1 = which(str_detect(colnames(rna), "ENSG"))
sums = apply(rna[,z1], 2, sum)
z = which(sums ==0)
if(!(length(z)==0)){
rm = names(sums)[z]}

#2. remove those with MAD < 0? 
#1. remove those not expressed at all
z1 = which(str_detect(colnames(rna), "ENSG"))
#sums = apply(rna[,z1], 2, mad)
#z = which(sums <= 0)
#rna = rna[,-z]

#2. cluster cancer types, ie -> each dot should be cancer type 
library("FactoMineR")
z1 = which(str_detect(colnames(rna), "ENSG"))
logged_rna = rna
logged_rna[,z1] = log1p(logged_rna[,z1])


#logged - DONE! 
#pdf("logged_nonMAD0lncRNAs_cands_PCA_plots_Aug27.pdf", width=12)
#autoplot(prcomp(logged_rna[,z1]), data = logged_rna, colour = 'type')+ 
#scale_colour_manual(values = mypal)
#dev.off()

#try t-SNE 
library(Rtsne)
iris_unique <- unique(iris) # Remove duplicates
iris_matrix <- as.matrix(iris_unique[,1:4])
set.seed(42) # Set a seed if you want reproducible results
tsne_out <- Rtsne(logged_rna[,z1]) # Run TSNE
saveRDS(tsne_out, file="tsne_out_28_cancers_Oct1.rds")

tsne_plot = readRDS("tsne_out_28_cancers_Oct1.rds")
# Show the objects in the 2D tsne representation
tsne_plot <- data.frame(x = tsne_out$Y[,1], y = tsne_out$Y[,2], col = logged_rna$type)

z = which(duplicated(tsne_plot$col))
tsne_plot$label[z] ="no"
tsne_plot$label[-z] = "yes"

z = (which(tsne_plot$label == "yes"))
for(i in 1:length(z)){
	canc = tsne_plot$col[z[i]]
	tsne_plot$label[z[i]] = canc
}


#assign color to variable ******




pdf("tSNE_28_cancers_oct1_all_lncs.pdf", width=9)
ggplot(tsne_plot,aes(x, y, label = label)) + geom_point(aes(x=x, y=y, color=col)) + scale_colour_manual(values=mypal)+
geom_text_repel(data = subset(tsne_plot, !(label == "no")))
dev.off()


#tSNA just candidates 
z1 = which(colnames(logged_rna) %in% allCands$gene)
tsne_out <- Rtsne(logged_rna[,z1]) # Run TSNE
saveRDS(tsne_out, file="just_cands_tsne_out_28_cancers_Oct1.rds")

tsne_plot <- data.frame(x = tsne_out$Y[,1], y = tsne_out$Y[,2], col = logged_rna$type)

z = which(duplicated(tsne_plot$col))
tsne_plot$label[z] ="no"
tsne_plot$label[-z] = "yes"

z = (which(tsne_plot$label == "yes"))
for(i in 1:length(z)){
	canc = tsne_plot$col[z[i]]
	tsne_plot$label[z[i]] = canc
}


pdf("tSNE_28_cancers_oct1_just_can_lncs.pdf", width=9)
ggplot(tsne_plot,aes(x, y, label = label)) + geom_point(aes(x=x, y=y, color=col)) + scale_colour_manual(values=mypal)+
geom_text_repel(data = subset(tsne_plot, !(label == "no")))
dev.off()


#NOT RUN ---------------------------------------------------------------
#not logged & scaled 
#pdf("scaled_nonMAD0lncRNAs_cands_PCA_plots_June25.pdf", width=12)
#z1 = which(str_detect(colnames(rna), "ENSG"))
#autoplot(prcomp(rna[,z1], scale. = TRUE), data = rna, colour = 'type')
#dev.off()
#PCA--------------------------------------------------------------------
#NOT RUN END------------------------------------------------------------


#-------------------------------------------------------------------
#-----------------PCA using mean/lncRNA/cancer to get --------------
#-----------------		32 points in the end          --------------
#-------------------------------------------------------------------

#1. get median for each gene within each cancer - Logged version

z = which(str_detect(colnames(logged_rna), "ENSG"))
means_rna = aggregate(logged_rna[, z], list(logged_rna$type), mean)
colnames(means_rna)[1] = "Cancer"
z = which(str_detect(colnames(means_rna), "ENSG"))
rownames(means_rna) = means_rna$Cancer

color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
z1 = which(str_detect(color, c("yellow")))
z2 = which(str_detect(color, c("snow")))
z3 = which(str_detect(color, c("lemon")))
z4 = which(str_detect(color, c("light")))
z5 = which(str_detect(color, c("pale")))
z6 = which(str_detect(color, c("white")))
z7 = which(str_detect(color, c("sandy")))
z8 = which(str_detect(color, c("blush")))
z9 = which(str_detect(color, c("corn")))
z10 = which(str_detect(color, c("darksea")))
z11 = which(str_detect(color, c("darkolive")))
z12 = which(str_detect(color, c("gainsboro")))
z13 = which(str_detect(color, c("lavender")))
z15 = which(str_detect(color, c("seashell")))
z16 = which(str_detect(color, c("honeydew")))

color = color[-c(z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12, z13, z15, z16)]

mypal = sample(color, 32)
z = which(mypal == "ivory2")
mypal[z] = "black"

#saveRDS(mypal, file="palette_32_cancer_types.rds")

mypal = readRDS("palette_32_cancer_types.rds")

pdf("cancer_logged_means_PCA_plots_july9.pdf", width=6, height=6)
g = autoplot(prcomp(means_rna[,z], scale. = TRUE), data = means_rna, label=TRUE, colour = 'Cancer', label.size = 3, label.repel=TRUE) +
theme_minimal() + scale_fill_manual(values = mypal) + 
scale_color_manual(values =mypal) + ggtitle("PCA using cancer mean log1p(lncRNA expression)")
ggpar(g, legend="none")
dev.off()

#2. get median for each gene within each cancer - Non-Logged version

z = which(str_detect(colnames(rna), "ENSG"))
means_rna = aggregate(rna[, z], list(rna$type), mean)
colnames(means_rna)[1] = "Cancer"
z = which(str_detect(colnames(means_rna), "ENSG"))
rownames(means_rna) = means_rna$Cancer

pdf("cancer_NOTloggd_means_PCA_plots_july9.pdf", width=6, height=6)
g = autoplot(prcomp(means_rna[,z]), data = means_rna, label=TRUE, colour = 'Cancer', label.size = 3, label.repel=TRUE) +
theme_minimal() + scale_fill_manual(values = mypal) + 
scale_color_manual(values =mypal) + ggtitle("PCA using cancer mean lncRNA expression")
ggpar(g, legend="none")
dev.off()


#-------------------------------------------------------------------
#-----------------PCA using median/lncRNA/cancer to get ------------
#-----------------		32 points in the end          --------------
#-------------------------------------------------------------------

z = which(str_detect(colnames(logged_rna), "ENSG"))
means_rna = aggregate(logged_rna[, z], list(logged_rna$type), median)
colnames(means_rna)[1] = "Cancer"
z = which(str_detect(colnames(means_rna), "ENSG"))
rownames(means_rna) = means_rna$Cancer

#remove columns with zero variance
vars = apply(means_rna[,z], 2, var)
rm = names(vars)[which(vars ==0)]
means_rna = means_rna[,-(which(colnames(means_rna) %in% rm))]
z = which(str_detect(colnames(means_rna), "ENSG"))

pdf("cancer_loggef_median_PCA_plots_June29.pdf", width=6, height=6)
g = autoplot(prcomp(means_rna[,z]), data = means_rna, label=TRUE, colour = 'Cancer', label.size = 3, label.repel=TRUE) +
theme_minimal() + scale_fill_manual(values = mypal) + 
scale_color_manual(values =mypal) + ggtitle("PCA using cancer median log1p(lncRNA expression)")
ggpar(g, legend="none")
dev.off()


#-------------------------------------------------------------------
#-----------------PCA using MAD/lncRNA/cancer to get ---------------
#-----------------		32 points in the end          --------------
#-------------------------------------------------------------------


z = which(str_detect(colnames(logged_rna), "ENSG"))
means_rna = aggregate(logged_rna[, z], list(logged_rna$type), mad)
colnames(means_rna)[1] = "Cancer"
z = which(str_detect(colnames(means_rna), "ENSG"))
rownames(means_rna) = means_rna$Cancer

#remove columns with zero variance
vars = apply(means_rna[,z], 2, var)
rm = names(vars)[which(vars ==0)]
means_rna = means_rna[,-(which(colnames(means_rna) %in% rm))]
z = which(str_detect(colnames(means_rna), "ENSG"))

pdf("cancer_loggef_mad_PCA_plots_June29.pdf", width=6, height=6)
g = autoplot(prcomp(means_rna[,z]), data = means_rna, label=TRUE, colour = 'Cancer', label.size = 3, label.repel=TRUE) +
theme_minimal() + scale_fill_manual(values = mypal) + 
scale_color_manual(values =mypal) + ggtitle("PCA using cancer mad log1p(lncRNA expression)")
ggpar(g, legend="none")
dev.off()



#-------------------------------------------------------------------
#-----------------PCA using just 166 lncRNA candidates--------------
#-------------------------------------------------------------------

#subset to cancers,  n = 23
colnames(canc_conv)[2] = "cancer"
allCands = merge(allCands, canc_conv, by="cancer")

z = which(logged_rna$type %in% allCands$type)

#subset to lncRNAs, keep cancer type 
cands_logged_rna = logged_rna[z,]

z1 = which(colnames(cands_logged_rna) %in% allCands$gene)
z2 = which(colnames(cands_logged_rna) %in% c("type", "patient", "type"))

cands_logged_rna = cands_logged_rna[,c(z1, z2)]
rownames(cands_logged_rna) = rna$patient

z = which(str_detect(colnames(cands_logged_rna), "ENSG"))
means_rna = aggregate(cands_logged_rna[, z], list(cands_logged_rna$type), mean)
colnames(means_rna)[1] = "Cancer"
z = which(str_detect(colnames(means_rna), "ENSG"))
rownames(means_rna) = means_rna$Cancer

pdf("cancer_loggef_mean_just_166_cands_PCA_plots_June29.pdf", width=6, height=6)
g = autoplot(prcomp(means_rna[,z]), data = means_rna, label=TRUE, colour = 'Cancer', label.size = 3, label.repel=TRUE) +
theme_minimal() + scale_fill_manual(values = mypal) + 
scale_color_manual(values =mypal) + ggtitle("PCA using 166 cands cancer \nmean log1p(lncRNA expression)")
ggpar(g, legend="none")
dev.off()






