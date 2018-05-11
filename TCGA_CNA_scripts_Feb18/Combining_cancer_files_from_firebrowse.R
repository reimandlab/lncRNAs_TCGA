###Combining_cancer_files_from_firebrowse.R
library(data.table)
library(plyr)
library(dplyr)
library(stringr)

#1. LIHC
lihc_cna = fread("LIHC.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt")
lihc_cna$canc = "liver"
lihc_cna$Chromosome = paste("chr", lihc_cna$Chromosome, sep="")

#2. KIRC
kirc_cna = fread("KIRC.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt")
kirc_cna$canc = "kidney"
kirc_cna$Chromosome = paste("chr", kirc_cna$Chromosome, sep="")

#3. PAAD
paad_cna = fread("PAAD.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt")
paad_cna$canc = "pancreas"
paad_cna$Chromosome = paste("chr", paad_cna$Chromosome, sep="")

#4. OV 
ov_cna = fread("OV.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt")
ov_cna$canc = "ovary"
ov_cna$Chromosome = paste("chr", ov_cna$Chromosome, sep="")

#5. BRCA
brca_cna = fread("BRCA.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt")
brca_cna$canc = "brca"
brca_cna$Chromosome = paste("chr", brca_cna$Chromosome, sep="")

#6. LUAD
luad_cna = fread("LUAD.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt")
luad_cna$canc = "luad"
luad_cna$Chromosome = paste("chr", luad_cna$Chromosome, sep="")

#7. acc
acc_cna = fread("ACC.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt")
acc_cna$canc = "acc"
acc_cna$Chromosome = paste("chr", acc_cna$Chromosome, sep="")

#8. CESC
cesc_cna = fread("CESC.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt")
cesc_cna$canc = "cesc"
cesc_cna$Chromosome = paste("chr", cesc_cna$Chromosome, sep="")

#9. ESCA
esca_cna = fread("ESCA.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt")
esca_cna$canc = "esca"
esca_cna$Chromosome = paste("chr", esca_cna$Chromosome, sep="")

#10. GBM
gbm_cna = fread("GBM.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt")
gbm_cna$canc = "gbm"
gbm_cna$Chromosome = paste("chr", gbm_cna$Chromosome, sep="")

#11. HNSC
hnsc_cna = fread("HNSC.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt")
hnsc_cna$canc = "hnsc"
hnsc_cna$Chromosome = paste("chr", hnsc_cna$Chromosome, sep="")

#12. KIRP
kirp_cna = fread("KIRP.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt")
kirp_cna$canc = "kirp"
kirp_cna$Chromosome = paste("chr", kirp_cna$Chromosome, sep="")

#13. LGG
lgg_cna = fread("LGG.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt")
lgg_cna$canc = "lgg"
lgg_cna$Chromosome = paste("chr", lgg_cna$Chromosome, sep="")

#14. MESO
meso_cna = fread("MESO.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt")
meso_cna$canc = "meso"
meso_cna$Chromosome = paste("chr", meso_cna$Chromosome, sep="")

#15. READ
read_cna = fread("READ.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt")
read_cna$canc = "read"
read_cna$Chromosome = paste("chr", read_cna$Chromosome, sep="")

#16. SARC
sarc_cna = fread("SARC.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt")
sarc_cna$canc = "sarc"
sarc_cna$Chromosome = paste("chr", sarc_cna$Chromosome, sep="")

#17. STAD
stad_cna = fread("STAD.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt")
stad_cna$canc = "stad"
stad_cna$Chromosome = paste("chr", stad_cna$Chromosome, sep="")

#18. UVM
uvm_cna = fread("UVM.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt")
uvm_cna$canc = "uvm"
uvm_cna$Chromosome = paste("chr", uvm_cna$Chromosome, sep="")

#19. THCA
thca_cna = fread("THCA.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt")
thca_cna$canc = "thca"
thca_cna$Chromosome = paste("chr", thca_cna$Chromosome, sep="")

#20. BLCA
blca_cna = fread("BLCA.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt")
blca_cna$canc = "blca"
blca_cna$Chromosome = paste("chr", blca_cna$Chromosome, sep="")

#21. UCEC
ucec_cna = fread("UCEC.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt")
ucec_cna$canc = "ucec"
ucec_cna$Chromosome = paste("chr", ucec_cna$Chromosome, sep="")

#22. LUSC
lusc_cna = fread("LUSC.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt")
lusc_cna$canc = "lusc"
lusc_cna$Chromosome = paste("chr", lusc_cna$Chromosome, sep="")

#23. COAD
coad_cna = fread("COAD.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt")
coad_cna$canc = "coad"
coad_cna$Chromosome = paste("chr", coad_cna$Chromosome, sep="")

allcnas = rbind(lihc_cna, kirc_cna, paad_cna, ov_cna, brca_cna,
	acc_cna, lusc_cna, luad_cna,
	cesc_cna,
	esca_cna,
	gbm_cna,
	hnsc_cna,
	kirp_cna,
	lgg_cna,
	meso_cna,
	read_cna,
	sarc_cna,
	stad_cna,
	uvm_cna,
	thca_cna,
	blca_cna,
	ucec_cna, 
	coad_cna)

saveRDS(allcnas, file="23cancers_types_ALL_TCGA_candidates_May11th_COPYNUMBER_CNAs.rds")



