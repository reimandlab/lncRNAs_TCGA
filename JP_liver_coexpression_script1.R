#---------------------------------------------------------
#JP_liver_coexpression_script1.R
#---------------------------------------------------------

#Author: Karina_Isaev
#Date_started: June 26th 2017
#dir: Thesis/pcawg_liver_JP

#------------
#Description:
#------------

#Have two RNA-Seq files, one for lncRNAs (all from UCSC including low confidence) 
#and one for PCGs. Conduct here LM and NGB regression to identify 
#signficantly co-expressed lncRNA-PCGs in liver cancer patients 
#using an array of confounders 
#------------
#confounders:
#------------
#[1]. PCG CNA
#[2]. lncRNA CNA
#[3]. Clinical features 
#[4]. is methylation available?

#SCRIPT1 - PROCESS GENES, VISUALIZE EXPRESSION OF GENES AND 
#SET UP DATA FRAME FOR REGRESSION ANALYSIS 

#---------------------------------------------------------
#Preamble
#---------------------------------------------------------

options(stringsAsFactors=F)

#---------------------------------------------------------
#Libraries
#---------------------------------------------------------
library(data.table)
library(survival)
library(ggplot2)
library(ggfortify)
library(cluster)
library(qqman)
library(dplyr)
library(rafalib)
library(RColorBrewer) 
library(genefilter)
library(gplots) ##Available from CRAN
library(survminer)
library(MASS)
library(Hmisc)
library(gProfileR)
library(wesanderson)
library(ggsci)
library(gridExtra)
library(ggpubr)

#---------------------------------------------------------
#Data
#---------------------------------------------------------

#These are gene expression files based on the original RNA-Seq
#PCAWG file
lnc <- readRDS("liver_jp_lncRNA_expression_6028.rds")
all <- readRDS("liver_jp_pcg_expression.rds")

#change feature column with gene names so that they are the rownames
rownames(lnc) <- lnc[,1] ; lnc <- lnc[,-1] #13479 lncRNAs as defined by ensembl "lincRNA", "antisense", "sense_intronic", "sense_overlapping"
rownames(all) <- all[,1] ; all <- all[,-1] #18039 PCGs as defined by ensembl "protein_coding"

#---------------------------------------------------------
#Processing lncRNAs and PCGs 
#---------------------------------------------------------

#1. REMOVE ALL THOSE THAT ARE 0 ACROSS THE BOARD BUT SAVE THEM
##lncs
sums <- apply(lnc, 1, sum)
z <- which(sums ==0) #845 have 0 expression in all patients 
zero_e <- rownames(lnc)[z] ; write.table(zero_e, file="jp_liver_0expressing_lncs.txt", quote=F, row.names=F)
#remove from expression file
lnc <- lnc[-z,]
##pcgs
sums <- apply(all, 1, sum)
z <- which(sums ==0) #845 have 0 expression in all patients 
zero_e <- rownames(all)[z] ; write.table(zero_e, file="jp_liver_0expressing_pcgs.txt", quote=F, row.names=F)
#remove from expression file
all <- all[-z,]

#2. CHECK HOW MANY ARE EXPRESSED ABOVE 0 IN 50%, 60%, 70%, 80%, 90% of patients 
check_e <- function(row, perc){
	z <- length(which(!(row<0.5)))
	check <- (z > (perc * length(row)))
	if(check){
		return("yes")
	}
	if(!(check)){
		return("no")
	}
}
lncs_res <- c()
e_50 <- apply(lnc, 1, function(x) check_e(x, perc=0.5)) ; lncs_res <- c(lncs_res, length(which(e_50 == "yes")))
e_60 <- apply(lnc, 1, function(x) check_e(x, perc=0.6)) ; lncs_res <- c(lncs_res, length(which(e_60 == "yes")))
e_70 <- apply(lnc, 1, function(x) check_e(x, perc=0.7)) ; lncs_res <- c(lncs_res, length(which(e_70 == "yes")))
e_80 <- apply(lnc, 1, function(x) check_e(x, perc=0.8)) ; lncs_res <- c(lncs_res, length(which(e_80 == "yes")))
e_90 <- apply(lnc, 1, function(x) check_e(x, perc=0.9)) ; lncs_res <- c(lncs_res, length(which(e_90 == "yes")))

pcg_res <- c()
e_50 <- apply(all, 1, function(x) check_e(x, perc=0.5)) ; pcg_res <- c(pcg_res, length(which(e_50 == "yes")))
e_60 <- apply(all, 1, function(x) check_e(x, perc=0.6)) ; pcg_res <- c(pcg_res, length(which(e_60 == "yes")))
e_70 <- apply(all, 1, function(x) check_e(x, perc=0.7)) ; pcg_res <- c(pcg_res, length(which(e_70 == "yes")))
e_80 <- apply(all, 1, function(x) check_e(x, perc=0.8)) ; pcg_res <- c(pcg_res, length(which(e_80 == "yes")))
e_90 <- apply(all, 1, function(x) check_e(x, perc=0.9)) ; pcg_res <- c(pcg_res, length(which(e_90 == "yes")))

#collect results
res <- as.data.frame(matrix(nrow=10, ncol=3)) 
colnames(res) <- c("perc", "num", "type")
res[1:10,1] <- rep(c(">50%", ">60%", ">70%", ">80%", ">90%"),2)
res[1:5,2] <- lncs_res 
res[6:10,2] <- pcg_res 
res[,3] <- c(rep("lnc", 5), rep("pcg", 5))

#3. CHECK HOW MANY ARE EXPRESSED ABOVE MEDIAN OF 1
meds <- apply(lnc, 1, median)
#get distribution of how many lncRNAs are above a certain median value
#ie, 0, 0.5, 1, 2, 3, 4 ...
meds <- as.data.frame(meds)
meds$gene <- rownames(meds)
#12,634 lncRNAs
z <- which(meds$meds < 1) ; 
#11911 have median Expression less than 1
#Don't plot them 
meds2 <- meds[-z,]
#723 remain - make histogram 
#pdf("lncRNAs_medians_greaterThan1.pdf", pointsize=4, width=8.5)
g <- gghistogram(meds2, x="meds", title="723 lncRNAs with Median Expression Greater Than 1FPKM", fill="lightgray", add = "median", rug=FALSE, 
	bins=50, xlab="Median", ylab="Number of lncRNAs")
g <- ggpar(g, xticks.by = 2)
#dev.off()

##greater than 2s
z <- which(meds2$meds < 2) ; 
#11911 have median Expression less than 1
#Don't plot them 
meds3 <- meds2[-z,]
#723 remain - make histogram 
#pdf("lncRNAs_medians_greaterThan2.pdf", pointsize=4, width=8.5)
g2 <- gghistogram(meds3, x="meds", title="317 lncRNAs with Median Expression Greater Than 2FPKM", fill="lightgray", add = "median", rug=FALSE, 
	bins=50, xlab="Median", ylab="Number of lncRNAs")
g2 <- ggpar(g2, xticks.by = 2)
#dev.off()

##greater than 3s
z <- which(meds3$meds < 3) ; 
#11911 have median Expression less than 1
#Don't plot them 
meds4 <- meds3[-z,]
#723 remain - make histogram 
#pdf("lncRNAs_medians_greaterThan2.pdf", pointsize=4, width=8.5)
g3 <- gghistogram(meds4, x="meds", title="169 lncRNAs with Median Expression Greater Than 3FPKM", fill="lightgray", add = "median", rug=FALSE, 
	bins=50, xlab="Median", ylab="Number of lncRNAs")
g3 <- ggpar(g3, xticks.by = 2)
#dev.off()

##greater than 4s
z <- which(meds4$meds < 4) ; 
#11911 have median Expression less than 1
#Don't plot them 
meds5 <- meds4[-z,]
#723 remain - make histogram 
#pdf("lncRNAs_medians_greaterThan2.pdf", pointsize=4, width=8.5)
g4 <- gghistogram(meds5, x="meds", title="105 lncRNAs with Median Expression Greater Than 4FPKM", fill="lightgray", add = "median", rug=FALSE, 
	bins=50, xlab="Median", ylab="Number of lncRNAs")
g4 <- ggpar(g4, xticks.by = 2)
#dev.off()

##greater than 5s
z <- which(meds5$meds < 5) ; 
#11911 have median Expression less than 1
#Don't plot them 
meds6 <- meds5[-z,]
#723 remain - make histogram 
#pdf("lncRNAs_medians_greaterThan2.pdf", pointsize=4, width=8.5)
g5 <- gghistogram(meds6, x="meds", title="70 lncRNAs with Median Expression Greater Than 5FPKM", fill="lightgray", add = "median", rug=FALSE, 
	bins=50, xlab="Median", ylab="Number of lncRNAs")
g5 <- ggpar(g5, xticks.by = 2)
#dev.off()


pdf("arranged_median_summaries.pdf", pointsize=3, width=10, height=12)
grid.arrange(g, g2, g3, g4, g5, ncol = 1, nrow=5)
dev.off()

##Seems like median expression of 4 is pretty good 
#have 105 like that
#make volcano plots for them to get an idea of each lncRNAs expresion 

#subset to those genes
genes <- meds5$gene
z <- which(rownames(lnc) %in% genes)
median4 <- lnc[z,] #105

plot_list <- list()

for(i in 1:nrow(median4)){
geneE <- median4[i,]
gene_name <- rownames(median4)[i]
names <- names(geneE)
gene <- as.data.frame(as.numeric(geneE))
gene$pat <- names
colnames(gene)[1] <- "Expression"
a<-ggviolin(gene,y="Expression",
             xlab="",ylab="FPKM", add = "boxplot")+
             ggtitle(paste(gene_name, "Expression Distribution"))+
             theme(plot.title=element_text(hjust=0.5))

   plot_list[[i]]<-a


}#end for loop

ml<-marrangeGrob(plot_list,nrow=2,ncol=2)
ggsave("lncs_wMedian_greaterThan4_Violins.pdf",ml)


tier4_lncs <- rownames(median4) ; write.table(tier4_lncs, file="meds_greaterthan4.txt", quote=F, row.names=F)




#divide lncRNAs into those with median expression above 1 
#and those with expression less than that 
meds <- apply(lnc, 1, median)
z <- which(meds >=1)
lncs_tier1 <- lnc[z,] #723
lncs_tier2 <- lnc[-z,] #11,911

#PCGS
meds <- apply(all, 1, median)
z <- which(meds >= 1) ; 
#divide lncRNAs into those with median expression above 1 
#and those with expression less than that 
pcgs_tier1 <- all[z,] #9,800
pcgs_tier2 <- all[-z,] #7,831

#4. SAVE LIST OF MEDIANS IN EACH TIER FOR PLOTTING BOXPLOTS 
#tier1 lncRNAs
meds_tier1 <- apply(lncs_tier1, 1, median)
tier1 <- as.data.frame(matrix(nrow=723, ncol=2)) ; colnames(tier1) <- c("Gene", "Expression")
tier1[,1] <- names(meds_tier1) ; tier1[,2] <- meds_tier1
#tier2 pcgs
meds_tier1_pcgs <- apply(pcgs_tier1, 1, median)
tier1_pcgs <- as.data.frame(matrix(nrow=9800, ncol=2)) ; colnames(tier1_pcgs) <- c("Gene", "Expression")
tier1_pcgs[,1] <- names(meds_tier1_pcgs) ; tier1_pcgs[,2] <- meds_tier1_pcgs
#summary(tier1_pcgs), some genes have median E of 13963.279  ? outliers?
#remove all those PCGs with medians of 500 
z <- which(tier1_pcgs$Expression > 500)
tier1_pcgs <- tier1_pcgs[-z,] #9760 are left 

#5. FOR TIER2 LNCRNAS, REMOVE THE ONES THAT ARE EXPRESSED AT 0 IN MORE THAN 20% OF PEOPLE
#tier2 lncRNAs
tier2_rm <- apply(lncs_tier2, 1, function(x) check_e(x, perc=0.2)) ; tier2_rm_res <- which(tier2_rm == "yes")
#8044/11911 expressed above 0 in more than 20% of samples, keep these ones 
lncs_tier2 <- lncs_tier2[tier2_rm_res,]
meds_tier2 <- apply(lncs_tier2, 1, median)
tier2 <- as.data.frame(matrix(nrow=8044, ncol=2)) ; colnames(tier2) <- c("Gene", "Expression")
tier2[,1] <- names(meds_tier2) ; tier2[,2] <- meds_tier2
#tier2 pcgs
tier2_pcgs_rm <- apply(pcgs_tier2, 1, function(x) check_e(x, perc=0.2)) ; tier2_pcgs_rm_res <- which(tier2_pcgs_rm == "yes")
#6458/7831 expressed above 0 in more than 20% of samples, keep these ones 
pcgs_tier2 <- pcgs_tier2[tier2_pcgs_rm_res,]
meds_tier2_pcgs <- apply(pcgs_tier2, 1, median)
tier2_pcgs <- as.data.frame(matrix(nrow=6458, ncol=2)) ; colnames(tier2_pcgs) <- c("Gene", "Expression")
tier2_pcgs[,1] <- names(meds_tier2_pcgs) ; tier2_pcgs[,2] <- meds_tier2_pcgs

#---------------------------------------------------------
#Visualization of Data 
#---------------------------------------------------------
#palette using ggsci package and npg nature style palette 

#2. CHECK HOW MANY ARE EXPRESSED ABOVE 0 IN 50%, 60%, 70%, 80%, 90% of patients 
mypal = pal_npg("nrc", alpha = 0.7)(10)
pdf("arranged_test.pdf", pointsize=5)
p1 = ggplot(res, aes(x=factor(perc), y=num, fill=factor(type))) + 
  geom_bar(stat="identity", colour="black") +  theme_bw() + geom_text(aes(label = num), position = position_stack(vjust = 0.5), size=4) + labs( x="Percent of Patients", y= "Number of Genes with non-zero Expression", title="Data for 68 JP Liver Cancer Patients")+ theme(legend.title=element_blank())
  #scale_fill_manual(values = sample(mypal,2))
p1 + scale_fill_npg() #use nature style colours 
dev.off()

#3A. VISUALIZING LNCRNA TIER1 AND TIER2 EXPRESSION 
#tier1
p1 = ggdensity(tier1, x="Expression", add = "median", rug=FALSE, palette = c("#E64B35B2"), fill="#E64B35B2", title="Gene Expression Distribution Tier 1 lncRNAs", xlab="Median Expression", ylab="Denstiy")+
scale_x_continuous(breaks = c(0,1,2,4,10,20,30,40,50,60,65))
p2 = ggviolin(tier1, y = "Expression", palette = c("#E64B35B2"), fill="#E64B35B2", title="Gene Expression Distribution Tier 1 lncRNAs", add = "boxplot", xlab="", ylab="Median Expression") + 
theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())
#tier2 
p3 = ggdensity(tier2, x="Expression", add = "median", rug=FALSE, palette = c("#91D1C2B2"), fill="#91D1C2B2", title="Gene Expression Distribution Tier 2 lncRNAs", xlab="Median Expression", ylab="Denstiy")+
scale_x_continuous(breaks=seq(0, 1, 0.1), limits=c(0, 1))
p4 = ggboxplot(tier2, y = "Expression", palette = c("#91D1C2B2"), fill="#91D1C2B2", title="Gene Expression Distribution Tier 2 lncRNAs", add = "none",xlab="", ylab="Median Expression") + 
theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())
pdf("arranged_tier1_tier2_summaryE.pdf", pointsize=3, width=10, height=8)
grid.arrange(p1, p3, p2, p4, ncol = 2, nrow=2)
dev.off()

#3B. VISUALIZING LNCRNA TIER1 AND TIER2 EXPRESSION - LOGGED 
tier1$Expression <- log1p(tier1$Expression)
tier2$Expression <- log1p(tier2$Expression)
#tier1
p1 = ggdensity(tier1, x="Expression", add = "median", rug=FALSE, palette = c("#E64B35B2"), fill="#E64B35B2", title="Gene Expression Distribution Tier 1 lncRNAs", xlab="Median Expression", ylab="Denstiy")+
scale_x_continuous(breaks = c(0,1,2,4,10,20,30,40,50,60,65))
p2 = ggviolin(tier1, y = "Expression", palette = c("#E64B35B2"), fill="#E64B35B2", title="Gene Expression Distribution Tier 1 lncRNAs", add = "boxplot", xlab="", ylab="Median Expression") + 
theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())
#tier2 
p3 = ggdensity(tier2, x="Expression", add = "median", rug=FALSE, palette = c("#91D1C2B2"), fill="#91D1C2B2", title="Gene Expression Distribution Tier 2 lncRNAs", xlab="Median Expression", ylab="Denstiy")+
scale_x_continuous(breaks=seq(0, 1, 0.1), limits=c(0, 1))
p4 = ggboxplot(tier2, y = "Expression", palette = c("#91D1C2B2"), fill="#91D1C2B2", title="Gene Expression Distribution Tier 2 lncRNAs", add = "none",xlab="", ylab="Median Expression") + 
theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())
pdf("arranged_tier1_tier2_summaryE_logged.pdf", pointsize=3, width=10, height=8)
grid.arrange(p1, p3, p2, p4, ncol = 2, nrow=2)
dev.off()

#3C. VISUALIZING LNCRNA TIER1 AND TIER2 EXPRESSION - PCGs
#tier1
p1 = ggdensity(tier1_pcgs, x="Expression", add = "median", font.xtick=10, rug=FALSE, palette = c("#E64B35B2"), fill="#E64B35B2", title="Gene Expression Distribution Tier 1 PCGs", xlab="Median Expression", ylab="Denstiy")+
scale_x_continuous(breaks = c(1,20,40,80,100,200,300,400))
p2 = ggviolin(tier1_pcgs, y = "Expression", palette = c("#E64B35B2"), fill="#E64B35B2", title="Gene Expression Distribution Tier 1 PCGs", add = "boxplot", xlab="", ylab="Median Expression") + 
theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())
#tier2 
p3 = ggdensity(tier2_pcgs, x="Expression", add = "median", rug=FALSE, palette = c("#91D1C2B2"), fill="#91D1C2B2", title="Gene Expression Distribution Tier 2 PCGs", xlab="Median Expression", ylab="Denstiy")+
scale_x_continuous(breaks=seq(0, 1, 0.1), limits=c(0, 1))
p4 = ggboxplot(tier2_pcgs, y = "Expression", palette = c("#91D1C2B2"), fill="#91D1C2B2", title="Gene Expression Distribution Tier 2 PCGs", add = "none",xlab="", ylab="Median Expression") + 
theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())
pdf("arranged_tier1_tier2_summaryE_pcgs.pdf", pointsize=3, width=10, height=8)
grid.arrange(p1, p3, p2, p4, ncol = 2, nrow=2)
dev.off()

#3D. VISUALIZING LNCRNA TIER1 AND TIER2 EXPRESSION - PCGs - LOGGED
tier1_pcgs$Expression <- log1p(tier1_pcgs$Expression)
tier2_pcgs$Expression <- log1p(tier2_pcgs$Expression)
#tier1
p1 = ggdensity(tier1_pcgs, x="Expression", add = "median", font.xtick=10, rug=FALSE, palette = c("#E64B35B2"), fill="#E64B35B2", title="Gene Expression Distribution Tier 1 PCGs", xlab="Median Expression", ylab="Denstiy")+
scale_x_continuous(breaks = c(0,1,2,4,5,6))
p2 = ggviolin(tier1_pcgs, y = "Expression", palette = c("#E64B35B2"), fill="#E64B35B2", title="Gene Expression Distribution Tier 1 PCGs", add = "boxplot", xlab="", ylab="Median Expression") + 
theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())
#tier2 
p3 = ggdensity(tier2_pcgs, x="Expression", add = "median", rug=FALSE, palette = c("#91D1C2B2"), fill="#91D1C2B2", title="Gene Expression Distribution Tier 2 PCGs", xlab="Median Expression", ylab="Denstiy")+
scale_x_continuous(breaks=seq(0, 1, 0.1), limits=c(0, 1))
p4 = ggboxplot(tier2_pcgs, y = "Expression", palette = c("#91D1C2B2"), fill="#91D1C2B2", title="Gene Expression Distribution Tier 2 PCGs", add = "none",xlab="", ylab="Median Expression") + 
theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())
pdf("arranged_tier1_tier2_summaryE_logged_pcgs.pdf", pointsize=3, width=10, height=8)
grid.arrange(p1, p3, p2, p4, ncol = 2, nrow=2)
dev.off()

#---------------------------------------------------------
#Save final litsts of genes in each tier to use further
#---------------------------------------------------------

tier1_lncs <- tier1$Gene ; write.table(tier1_lncs, file="tier1_lncs.txt", quote=F, row.names=F)
tier2_lncs <- tier2$Gene ; write.table(tier2_lncs, file="tier2_lncs.txt", quote=F, row.names=F)
tier1_pcgs <- tier1_pcgs$Gene ; write.table(tier1_pcgs, file="tier1_pcgs.txt", quote=F, row.names=F)
tier2_pcgs <- tier2_pcgs$Gene ; write.table(tier2_pcgs, file="tier2_pcgs.txt", quote=F, row.names=F)