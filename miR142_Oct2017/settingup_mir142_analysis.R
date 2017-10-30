###settingup_mir142_analysis.R

###Libraries
source("source_file.R")
library(GenomicRanges)

###Data 

#1. ppl with mutations in miR142
load("mirna142_muts_for_Lina_250117.rsav")
#20 mutations, 17 people in total 
muts_ppl = mirna142_muts$mirna_pre 

#2. Gene expression 
exp = read.table("RNA-seq-new-gencode-only-consolidated-rows-only-mir142-targets.txt")

#3. linker file 
link = read.csv("PCAWGMay2016DataReleasev1_v1.1_v1.2_v1.3_v1.4release_may2016v1.4.csv")

#4. cancer subtype 
canc = fread("pcawg_specimen_histology_August2016_v6.tsv")

#5. all samples that have gone through mutation analysis 
#*you should remove any samples not found in this list from your analysis, so that downstream statistical testing is fair
mutanalysis = fread("mutated_samples.txt")

#6. identified mir142 targets comrpise set of genes tests with annotations 
genes = fread("file_for_network_w_target_type_and_cancer_genes.txt")

#7. copy number data for lymphoma patients (196/208) from JR's file obtained on October 30th 
cna = readRDS("196pcawg_lymphoma_cna_matrix_fromJR_oct30.rds") 

#8. UCSC coordinates 
ucsc <- fread("UCSC_hg19_gene_annotations_downlJuly27byKI.txt", data.table=F)

###Converting RNA-Seq specimen ids to donor ids 
#RNA seq specimen IDs to donor ID 


###CNA figure out how segments match to gene 
segments = rownames(cna) 
segments = str_split_fixed(segments, "::", n=2)

#make sure it is still atttached to patients so that later doesn't get lost, order of cnas...
cna = cbind(segments, cna)

#some rows are segments and some are gene IDs, let's divide those into seperate dataframes

segments1 = cna[which(cna[,1] == "ENCODEmerge"),]
segments2 =  cna[which(!(cna[,1] == "ENCODEmerge")),]

#segments 1
chr = str_split_fixed(segments1[,2], ":", n=4)
coordinates = str_split_fixed(chr[,2], "-", n=2)
chr = chr[,1]
segments1 = cbind(segments1, chr, coordinates)
segments1 = as.data.frame(segments1)
segments1$strand = "*"

#segments 2
chr = str_split_fixed(segments2[,2], ":", n=4)
coordinates = str_split_fixed(chr[,2], "-", n=2)
chr = chr[,1]
segments2 = cbind(segments2, chr, coordinates)
segments2 = as.data.frame(segments2)
segments2$strand = "*"

#seperate gene ID from the rest 
gene = str_split_fixed(segments2[,2], "::", n=3)
gene = cbind(segments2, gene)
gene = gene[which(gene[,204] %in% c("gencode")),]

#basically ignoring segments 1 from now as they are just segments 
#in final CNA file, gene -- 25,687 unique Hugo Gene Symbols
#25765 unique Ensembl gene ids 

##subset to genes we are interested in 
gene <- gene[which(gene[,205] %in% genes$all.targets), ] 


