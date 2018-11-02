#GTEX_data0.R

###---------------------------------------------------------------
###Load Libraries
###---------------------------------------------------------------

source("source_file.R")
library(ggrastr)

###---------------------------------------------------------------
###Load Data 
###---------------------------------------------------------------

#1. Gene Reads
gene_reads = fread("GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct")

#2. Transcript expected count 
#transcript = fread("GTEx_Analysis_2016-01-15_v7_RSEMv1.2.22_transcript_expected_count.txt")

#3. Sample attributes 
#sample_atts = fread("GTEx_Analysis_v7_Annotations_SampleAttributesDD.xlsx")
atts = fread("GTEx_v7_Annotations_SampleAttributesDS.txt")

atts = filter(atts, SAMPID %in% colnames(gene_reads))
atts = as.data.table(atts)
#for now just need liver
#atts = filter(atts, SMTS == "Liver")
atts = atts[,c(1, 6)]

cols = colnames(gene_reads)[which(colnames(gene_reads) %in% atts$SAMPID)]
cols = c(c("Name", "Description"), cols)

#RNASeq data for all tissue types 
gene_reads = gene_reads[, ..cols]
#1. Gene IDs 
extract3 <- function(row){
	gene <- as.character(row[[1]])
	ens <- gsub("\\..*","",gene)
	return(ens)
}

gene_reads[,1] <- apply(gene_reads[,1:2], 1, extract3) ; 
#4. Subject phenotypes
#sub_pheno = fread("GTEx_Analysis_v7_Annotations_SubjectPhenotypesDD.xlsx")
pheno = fread("GTEx_v7_Annotations_SubjectPhenotypesDS.txt")

#5. UCSC Hg19
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
z <- which(duplicated(ucsc[,6]))
ucsc <- ucsc[-z,]

#Process
gene_reads = as.data.frame(gene_reads)
rownames(gene_reads) = gene_reads$Name
gene_reads = gene_reads[,3:ncol(gene_reads)]
#keep only ZKSCAN3 and TP53 
z = which(rownames(gene_reads) %in% c("ENSG00000189298", "ENSG00000141510"))
exp_data = gene_reads[z,]
exp_data = t(exp_data)
exp_data = as.data.frame(exp_data)
exp_data$SAMPID = rownames(exp_data)
exp_data = merge(exp_data, atts, by = "SAMPID")
exp_data[,2:3] = log1p(exp_data[,2:3])

###-------Correlation plot ZKSCAN3 and TP53 all tissues combined-------------

#Generate better palette
library(viridis)
library(RColorBrewer)
library(colorRamps)
l = length(unique(exp_data$SMTS)) * 2
cc2color = structure(names = unique(exp_data$SMTS),
         sample(colorRampPalette(brewer.pal(12, "Set1"))(l)[1:l %% 2 == 0]))

saveRDS(cc2color, "palette_w30_colours.rds")

# Remove confidence intervals
# Extend the regression lines
exp_data$SMTS = as.factor(exp_data$SMTS)
lm_model0 = lm(exp_data$ENSG00000141510 ~ exp_data$SMTS, data = exp_data)
lm_model1 = lm(exp_data$ENSG00000141510 ~  exp_data$ENSG00000189298 + exp_data$SMTS, data = exp_data)
anova(lm_model1, lm_model1)

cor.test(exp_data$ENSG00000141510 ~  exp_data$ENSG00000189298,method="spearman")$p.value


pdf("zkscan3_tp53_all_gtex_samples_juris_code.pdf")
ggplot(exp_data, aes(ENSG00000189298, ENSG00000141510, fill=SMTS, color=SMTS)) +
        geom_smooth(method='lm', aes(color=SMTS), se=F) +
        #geom_point(alpha=0.5, size=2, shape=21) +
        scale_x_continuous("ZKSCAN3 expression (log)") +
        scale_y_continuous("TP53 expression (log)") +
        ggtitle(paste0("Correlated expression of ZKSCAN3 and TP53\n",
                        "p_lm=", "< 2e-16", "; coef_lm=", signif(tidy(lm_model1)[2,2]))) +
        theme_bw()+
        geom_point_rast(alpha=0.5)
dev.off()

pdf("zkscan3_tp53_all_gtex_samples_point_rast.pdf")
g = ggplot(exp_data, aes(x=ENSG00000189298, y=ENSG00000141510, color=SMTS)) +
  #geom_point(alpha=0.5, size=2, shape=21) + 
  geom_point_rast(alpha=0.5, size=2, shape=21)+
  geom_smooth(method=lm, se=FALSE)+
   scale_color_manual(values = cc2color) + 
   scale_x_continuous("ZKSCAN3 expression (log)") +
        scale_y_continuous("TP53 expression (log)") + theme_bw() +
        ggtitle(paste0("Correlated expression of ZKSCAN3 and TP53\n",
                        "p_lm=", "< 2e-16", "; coef_lm=", signif(tidy(lm_model1)[2,2])))

ggpar(g, legend.title="Tissue")
dev.off()

pdf("tp53_zkscan3_all_gtex_samples_point_rast.pdf")
g = ggplot(exp_data, aes(x=ENSG00000141510, y=ENSG00000189298, color=SMTS)) +
  #geom_point(alpha=0.5, size=2, shape=21) + 
  geom_point_rast(alpha=0.5, size=2, shape=21)+
  geom_smooth(method=lm, se=FALSE)+
   scale_color_manual(values = cc2color) + 
   scale_x_continuous("TP53 expression (log)") +
        scale_y_continuous("ZKSCAN3 expression (log)") + theme_bw() +
        ggtitle(paste0("Correlated expression of TP53 and ZKSCAN3\n",
                        "p_lm=", "< 2e-16", "; coef_lm=", signif(tidy(lm_model1)[2,2])))

ggpar(g, legend.title="Tissue")
dev.off()


###-------Correlation plot NEAT1 and MALAT1 all tissues combined-------------

#keep only NEAT1 and MALAT1 
z = which(rownames(gene_reads) %in% c("ENSG00000245532", "ENSG00000251562"))
exp_data = gene_reads[z,]
exp_data = t(exp_data)
exp_data = as.data.frame(exp_data)
exp_data$SAMPID = rownames(exp_data)
exp_data = merge(exp_data, atts, by = "SAMPID")
exp_data[,2:3] = log1p(exp_data[,2:3])

# Remove confidence intervals
# Extend the regression lines
exp_data$SMTS = as.factor(exp_data$SMTS)
lm_model0 = lm(exp_data$ENSG00000251562 ~ exp_data$SMTS, data = exp_data)
lm_model1 = lm(exp_data$ENSG00000251562 ~  exp_data$ENSG00000245532 + exp_data$SMTS, data = exp_data)
anova(lm_model1, lm_model1)

pdf("neat1_malat1_all_gtex_samples_juris_code.pdf")
ggplot(exp_data, aes(ENSG00000245532, ENSG00000251562, fill=SMTS)) +
        geom_smooth(method='lm', aes(color=SMTS), se=F) +
        geom_point(alpha=0.5, size=2, shape=21) +
        scale_x_continuous("NEAT1 expression (log)") +
        scale_y_continuous("MALAT1 expression (log)") +
        ggtitle(paste0("Correlated expression of NEAT1 and MALAT1\n",
                        "p_lm=", "< 2e-16", "; coef_lm=", signif(tidy(lm_model1)[2,2]))) +
        theme_bw()
dev.off()

pdf("neat1_malat1_all_gtex_samples_point_rast.pdf")
g = ggplot(exp_data, aes(x=ENSG00000245532, y=ENSG00000251562, color=SMTS)) +
  #geom_point(alpha=0.5, size=2, shape=21) + 
  geom_point_rast(alpha=0.5, size=2, shape=21)+
  geom_smooth(method=lm, se=FALSE)+
  scale_color_manual(values = cc2color) + 
   scale_x_continuous("NEAT1 expression (log)") +
        scale_y_continuous("MALAT1 expression (log)") + theme_bw() +
        ggtitle(paste0("Correlated expression of NEAT1 and MALAT1\n",
                        "p_lm=", "< 2e-16", "; coef_lm=", signif(tidy(lm_model1)[2,2])))
ggpar(g, legend.title="Tissue")
dev.off()














