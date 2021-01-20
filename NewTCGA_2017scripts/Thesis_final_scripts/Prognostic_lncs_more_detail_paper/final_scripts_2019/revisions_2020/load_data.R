source("/u/kisaev/lncRNAs_TCGA/NewTCGA_2017scripts/Thesis_final_scripts/Prognostic_lncs_more_detail_paper/final_scripts_2019/revisions_2020/load_packages.R")
source("/u/kisaev/lncRNAs_TCGA/NewTCGA_2017scripts/Thesis_final_scripts/Prognostic_lncs_more_detail_paper/final_scripts_2019/revisions_2020/load_functions.R")

#1. lncRNA expression in different cancers
#rna = readRDS("lncRNAs_2019_manuscript/5919_lncRNAs_tcga_all_cancers_March13_wclinical_dataalldat.rds") #<- updated gbm cohort
old = readRDS("lncRNAs_2019_manuscript/5919_lncRNAs_tcga_all_cancers_March13_wclinical_dataalldat.rds") #<- updated gbm cohort
rna = readRDS("lncRNAs_2019_manuscript/5919_lncRNAs_tcga_all_cancers_Jan2021_wclinical_dataalldat.rds") #<- updated gbm cohort

#2. Protein coding genes expression in different cancers
pcg = readRDS("lncRNAs_2019_manuscript/19438_lncRNAs_tcga_all_cancers_Jan2021_wclinical_dataalldat.rds")

#3. Fantom data
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

hg38 = readRDS("lncRNAs_2019_manuscript/grch38_dat_annotables.rds")
hg38$chr = paste("chr", hg38$chr, sep="")

z = which(fantom$gene %in% hg38$ensgene)
fantom = fantom[z,]
hg38_lncs = as.data.table(filter(hg38, ensgene %in% fantom$gene, !(biotype == "protein_coding")))
colnames(hg38_lncs)[1] = "gene"

fantom=merge(fantom, hg38_lncs, by="gene")
colnames(fantom)[2]="hg19_name"
colnames(fantom)[8]="CAT_geneName"

###---------------------------------------------------------------
###Get lncRNAs for each cancer
###---------------------------------------------------------------

#Add cancer type
canc_conversion = readRDS("lncRNAs_2019_manuscript/9753_tcga_id_cancer_type_conversion.txt")
canc_conversion = as.data.frame(canc_conversion)
canc_conversion = canc_conversion[,c(2,4)]
colnames(canc_conversion)[2] = "patient"

rna = merge(rna, canc_conversion, by="patient")
pcg = merge(pcg, canc_conversion, by="patient")

cancers = as.list(unique(rna$Cancer))

#Remove any lncRNAs that are not expressed in any of the patients
z=which(str_detect(colnames(rna), "ENSG"))
sums = apply(rna[,z], 2, sum) #134 with 0 expression in ALL patients
zeroes = names(sums[which(sums ==0)]) #what are they?
z <- which(colnames(rna) %in% zeroes)
rna = rna[,-z]

#only keep lncRNAs that are hg38 in new fantom
z = which(str_detect(colnames(rna), "ENSG"))
lnc_rna = rna[,z]
clin = rna[,-z]
z = which(colnames(lnc_rna) %in% fantom$gene)
lnc_rna = lnc_rna[,z]
rna = cbind(lnc_rna, clin)
rna = as.data.table(rna)

dim(rna)
dim(pcg)

#clean up data
#1. remove discrepancy
z = which(rna$vital_status == "[Discrepancy]")
rna = rna[-z,]
z = which(is.na(as.numeric(rna$age_at_initial_pathologic_diagnosis)))
rna = rna[-z,]
z = which(is.na(as.numeric(rna$OS.time)))
rna = rna[-z,]
z = which(as.numeric(rna$OS.time) == 0)
rna = rna[-z,]
z = which(is.na(as.numeric(rna$PFI.time)) & !(rna$type=="LAML"))
rna = rna[-z,]
z = which((as.numeric(rna$PFI.time) == 0) & (!(rna$type=="LAML")))
rna = rna[-z,]

#fix clinical variables so only one level per variable has NA
z = which(rna$race %in% c("[Not Available]", "[Not Evaluated]", "[Unknown]"))
rna$race[z] = "unknown"

z = which(rna$clinical_stage %in% c("[Not Applicable]", "[Not Available]", "[Discrepancy]"))
rna$clinical_stage[z] = "unknown"

z = which(rna$histological_grade %in% c("[Unknown]", "[Not Available]", "[Discrepancy]"))
rna$histological_grade[z] = "unknown"

print(table(rna$race))

#fix clinical variables so only one level per variable has NA
z = which(pcg$race %in% c("[Not Available]", "[Not Evaluated]", "[Unknown]"))
pcg$race[z] = "unknown"

z = which(pcg$clinical_stage %in% c("[Not Applicable]", "[Not Available]", "[Discrepancy]"))
pcg$clinical_stage[z] = "unknown"

z = which(pcg$histological_grade %in% c("[Unknown]", "[Not Available]", "[Discrepancy]"))
pcg$histological_grade[z] = "unknown"

print(table(pcg$race))

#2. list of cancers to apply function to
cancers = as.list(unique(rna$Cancer))

#3. function that splits data into cancers
get_canc = function(canc){
        canc_data = rna[which(rna$Cancer == canc),]
        return(canc_data)
}

canc_datas = llply(cancers, get_canc)

canc_conv = unique(rna[,c("type", "Cancer")])

z = which(pcg$patient %in% rna$patient)
pcg = pcg[z,]

cols = colnames(rna)[which(colnames(rna) %in% colnames(pcg))]
all = merge(rna, pcg, by = cols)

t = as.data.table(table(all$type))
t = t[order(N)]
print(t)

#------DATA---------------------------------------------------------
#UCSC gene info
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)
#z <- which(ucsc$hg19.ensemblSource.source %in% c("antisense", "lincRNA", "protein_coding"))
#ucsc <- ucsc[z,]
z <- which(duplicated(ucsc[,8]))
ucsc <- ucsc[-z,]

#manual_missing_lncs = as.data.frame(matrix(ncol=9, nrow=3))
#manual_missing_lncs[1,]=c("ENSG00000247844", "na", "CCAT1", "chr8", 127207866, 127219088, "-1",  "lincRNA", "na")
#manual_missing_lncs[2,]=c("ENSG00000233800", "na", "RP11-295G24.4", "chr20", 44448777, 44450092, "1",  "sense_intronic", "na")
#manual_missing_lncs[3,]=c("ENSG00000237907", "na", "LINC01430", "chr8", 127207866, 127219088, "-1",  "antisense_RNA", "na")
#colnames(manual_missing_lncs)=colnames(hg38)

#hg38 = rbind(hg38, manual_missing_lncs)

#remove cancer types with less than 50 patients
pats_num = as.data.table(table(all$type))
pats_num = filter(pats_num, N <50)
canc_rm = pats_num$V1
all = all[-which(all$type %in% canc_rm),]
rna = rna[-which(rna$type %in% canc_rm),]
pcg = pcg[-which(pcg$type %in% canc_rm),]
rna = as.data.table(rna)
pcg = as.data.table(pcg)

print(table(all$type))

#color palette
colours_palette=readRDS("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNAs_2019_manuscript/23_cancers_color_palette.rds")

myColors=colours_palette$color
names(myColors)=colours_palette$type
myColors
colScale <- scale_colour_manual(name = "type",values = myColors)

print("done loading everything YAY")
