#GTEX_data0.R

###---------------------------------------------------------------
###Load Libraries
###---------------------------------------------------------------

source("source_file.R")

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
#keep only ID1 and EGFR 
z = which(rownames(gene_reads) %in% c("ENSG00000125968", "ENSG00000146648"))
exp_data = gene_reads[z,]
exp_data = t(exp_data)
exp_data = as.data.frame(exp_data)
exp_data$SAMPID = rownames(exp_data)
exp_data = merge(exp_data, atts, by = "SAMPID")
exp_data[,2:3] = log1p(exp_data[,2:3])

###-------Correlation plot ID1 and EGFR all tissues combined-------------

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
lm_model0 = lm(exp_data$ENSG00000146648 ~ exp_data$SMTS, data = exp_data)
lm_model1 = lm(exp_data$ENSG00000146648 ~  exp_data$ENSG00000125968 + exp_data$SMTS, data = exp_data)
anova(lm_model1, lm_model1)

pdf("ID1_EGFR_all_gtex_samples_juris_code.pdf")
ggplot(exp_data, aes(ENSG00000125968, ENSG00000146648, fill=SMTS)) +
        geom_smooth(method='lm', aes(color=SMTS), se=F) +
        geom_point(alpha=0.5, size=2, shape=21) +
        scale_x_continuous("ID1 expression (log)") +
        scale_y_continuous("EGFR expression (log)") +
        ggtitle(paste0("Correlated expression of ID1 and EGFR\n",
                        "p_lm=", "3.86e-14", "; coef_lm=", signif(tidy(lm_model1)[2,2]))) +
        theme_bw()
dev.off()

pdf("ID1_EGFR_all_gtex_samples.pdf")
g = ggplot(exp_data, aes(x=ENSG00000125968, y=ENSG00000146648, color=SMTS)) +
  geom_point(alpha=0.5, size=2, shape=21) + 
  geom_smooth(method=lm, se=FALSE)+
   scale_color_manual(values = cc2color) + 
   scale_x_continuous("ID1 expression (log)") +
        scale_y_continuous("EGFR expression (log)") + theme_bw() +
        ggtitle(paste0("Correlated expression of ID1 and EGFR\n",
                        "p_lm=", "3.86e-14", "; coef_lm=", signif(tidy(lm_model1)[2,2])))
ggpar(g, legend.title="Tissue")
dev.off()


###-------Correlation plot ID1 and EGFR all tissues combined-------------

#keep only ID1 and EGFR 
z = which(rownames(gene_reads) %in% c("ENSG00000125968", "ENSG00000146648"))
exp_data = gene_reads[z,]
exp_data = t(exp_data)
exp_data = as.data.frame(exp_data)
exp_data$SAMPID = rownames(exp_data)
exp_data = merge(exp_data, atts, by = "SAMPID")
exp_data[,2:3] = log1p(exp_data[,2:3])

# Remove confidence intervals
# Extend the regression lines
exp_data$SMTS = as.factor(exp_data$SMTS)
lm_model0 = lm(exp_data$ENSG00000125968 ~ exp_data$SMTS, data = exp_data)
lm_model1 = lm(exp_data$ENSG00000125968 ~  exp_data$ENSG00000146648 + exp_data$SMTS, data = exp_data)
anova(lm_model1, lm_model1)

pdf("id1_egfr_all_gtex_samples_juris_code.pdf")
ggplot(exp_data, aes(ENSG00000146648, ENSG00000125968, fill=SMTS)) +
        geom_smooth(method='lm', aes(color=SMTS), se=F) +
        geom_point(alpha=0.5, size=2, shape=21) +
        scale_x_continuous("EGFR expression (log)") +
        scale_y_continuous("ID1 expression (log)") +
        ggtitle(paste0("Correlated expression of EGFR and ID1\n",
                        "p_lm=", "3.86e- 14", "; coef_lm=", signif(tidy(lm_model1)[2,2]))) +
        theme_bw()
dev.off()

pdf("id1_egfr_all_gtex_samples.pdf")
g = ggplot(exp_data, aes(x=ENSG00000146648, y=ENSG00000125968, color=SMTS)) +
  geom_point(alpha=0.5, size=2, shape=21) + 
  geom_smooth(method=lm, se=FALSE)+
  scale_color_manual(values = cc2color) + 
   scale_x_continuous("EGFR expression (log)") +
        scale_y_continuous("ID1 expression (log)") + theme_bw() +
        ggtitle(paste0("Correlated expression of EGFR and ID1\n",
                        "p_lm=", "3.86e- 14", "; coef_lm=", signif(tidy(lm_model1)[2,2])))
ggpar(g, legend.title="Tissue")
dev.off()














