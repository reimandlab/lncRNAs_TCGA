library(data.table)
library(ggplot2)
library(ggpubr)
library(stringr)
library(dplyr)
library(gridExtra)
library(grid)
library(GenomicRanges)

#Data---------------------------------------------------

#[1] FMREs, 41 in total
#load("dfr_encode_results_signf.rsav")
#head(dfr_encode_results_signf)
#fmres = dfr_encode_results_signf

# 84 unique patients have the FMRE "ENCODEmerge::chr17:57914548-57924685::NA::NA" 

#[2] coding drivers - need to filter for cancer_type==“PANCANCER” and 
#element_type==“gc19_pc.cds”, 47 in total

#old file:
#load("all_results_signf.rsav")
#head(all_results_signf)
#all_results_signf = as.data.table(all_results_signf)
#all_results_signf = filter(all_results_signf, cancer_type == "PANCANCER")
#all_results_signf = filter(all_results_signf, element_type == "gc19_pc.cds")
#coding_drivers = all_results_signf

#cds mutations new june 12 
coding_drivers = fread("july2_cds_drivers.txt")

#[3] mutations in all CRMs, subset of these are FMREs

#old file:

#load("encode_merge_oct2016_mutations__PANCANCER_in_elements.rsav")
#head(mutations_in_elements)
#mutations_in_crms = mutations_in_elements

#new file June 12: 

load("july12_encode_merge__patient_element_snv_list.rsav")
mutations_in_crms = (patient_element_snv_list)

#crm mutations new june 12 
#fmre mutations 
fmres = fread("july12_fmre_drivers.txt")

#[4] mutations in all CDS, subset of these are CDS drivers
load("july12_gc19_pc.cds__patient_element_snv_list.rsav")
cds_mutations = patient_element_snv_list

#[5] all patients in cohort
load("patient2cancertype.rsav")
head(patient2cancer_type) #1844 all together 

patient_table = fread("patient_table.txt")

#[] cds_lnc_drivers

#[] july12_gc19_pc.cds__coords.rsav
load("july12_gc19_pc.cds__coords.rsav")
coords = element_coords
coords = as.data.table(coords)

#[] july12_lncrna.ncrna__coords.rsav
load("july12_lncrna.ncrna__coords.rsav")
lnc_coords = element_coords
lnc_coords = as.data.table(lnc_coords)

#[] july12_encode_merge__coords.rsav
load("july12_encode_merge__coords.rsav")
element_coords = element_coords
element_coords = as.data.table(element_coords)

#[] july12_lncrna.ncrna__patient_element_snv_list.rsav
load("july12_lncrna.ncrna__patient_element_snv_list.rsav")
lncrna_mutations = patient_element_snv_list

#[6] Clinical file 
clin <- fread("pcawg_specimen_histology_August2016_v6.tsv", data.table=F)
conversion <- fread("pcawgConversion.tsv", data.table=F)

#Analysis---------------------------------------------------

#1. Which patients have FMRE mutations  
z = which(names(mutations_in_crms) %in% fmres$id)
mutations_in_crms = mutations_in_crms[z]


#2. Which patients have coding mutations 
z = which(names(cds_mutations) %in% coding_drivers$id)
mutations_in_cds = cds_mutations[z]

#tp53 muts 
tp53_muts = readRDS("tp53_pcawg_muts.rds")

#patients with zkscan FMRE --> 6:27870028
z = which(names(mutations_in_crms) == "ENCODEmerge::chr6:27870028-27871319::NA::NA")
fmre_pats = mutations_in_crms[[z]]

#which TP53 mutations do these patients have?
z = which(tp53_muts$Donor_ID %in% fmre_pats)
tp53_fmre_muts = tp53_muts[z,]
pats = unique(tp53_fmre_muts$Donor_ID)
pats[13]
#[1] "DO28763" <--- supposingly has TP53 mutation but not in CDS file 
tp53_fmre_muts[which(tp53_fmre_muts$Donor_ID == pats[13])]

#add 5 lncRNAs to mutations in cds
lncs = c("lncrna.ncrna::gencode::RN7SK::ENSG00000202198.1",
"lncrna.ncrna::gencode::NEAT1::ENSG00000245532.4",
"lncrna.ncrna::gencode::MALAT1::ENSG00000251562.3",
"lncrna.ncrna::gencode::RPPH1::ENSG00000259001.2",
"lncrna.ncrna::gencode::Z95704.4::ENSG00000248302.2")

z = which(names(lncrna_mutations) %in% lncs)
lncrna_mutations = lncrna_mutations[z]
#bind it to the cds list
mutations_in_cds = c(mutations_in_cds, lncrna_mutations)

#Compare ratios------------------------------------------------

unique_cds = unique(names(mutations_in_cds)) #48 + 5 lncRNAs 
unique_fmre = unique(names(mutations_in_crms)) #30

#Chcekc do we expect overlaps?
all_lnc_cord = as.data.table(filter(lnc_coords, id %in% unique_cds))
colnames(all_lnc_cord)[1:3] = c("Chr", "Start", "End")
all_lnc_cord = makeGRangesFromDataFrame(all_lnc_cord)

all_cds_cord = as.data.table(filter(coords, id %in% unique_cds))
colnames(all_cds_cord)[1:3] = c("Chr", "Start", "End")
all_cds_cord = makeGRangesFromDataFrame(all_cds_cord)

all_fmre_cord = as.data.table(filter(element_coords, id %in% unique_fmre))
colnames(all_fmre_cord)[1:3] = c("Chr", "Start", "End")
all_fmre_cord = makeGRangesFromDataFrame(all_fmre_cord)

findOverlaps(all_fmre_cord, all_cds_cord)
findOverlaps(all_fmre_cord, all_lnc_cord) #so only overlaps between FMREs and lncRNAs 

results_pairs = as.data.frame(matrix(ncol=6)) ; colnames(results_pairs) = c("CDS_mut", "FMRE_mut", "fishers_pval", "fishers_OR", "num_overlap", "canc_overlapping_pats")

#for each cds/fmre combo
for(i in 1:length(unique_cds)){
	for(y in 1:length(unique_fmre)){
		pair = c(unique_cds[i], unique_fmre[y])
		#first make sure that they don't overlap!
		lnc = which(names(lncrna_mutations) %in% pair[1])
		if(length(lnc) == 0){
			cds_cord = as.data.table(filter(coords, id == pair[1]))
			colnames(cds_cord)[1:3] = c("Chr", "Start", "End")
			cds_cord = makeGRangesFromDataFrame(cds_cord)
		}
		if(!(length(lnc) == 0)){
			cds_cord = as.data.table(filter(lnc_coords, id == pair[1]))
			colnames(cds_cord)[1:3] = c("Chr", "Start", "End")
			cds_cord = makeGRangesFromDataFrame(cds_cord)
		}

		fmre_cord = as.data.table(filter(element_coords, id == pair[2]))
		colnames(fmre_cord)[1:3] = c("Chr", "Start", "End")
		fmre_cord = makeGRangesFromDataFrame(fmre_cord)
		check_overlap = findOverlaps(cds_cord, fmre_cord)
		l_overlap = length(check_overlap)
		print(paste("num_overlap =", l_overlap))

		if(l_overlap == 0){

		#get patients that have either of these mutations or both 
		#FMRE
		z1 = which(names(mutations_in_crms) %in% pair)
		pats_crms = as.data.frame(mutations_in_crms[z1])
		pats_crms = as.data.frame(pats_crms[!duplicated(pats_crms), ])
		pats_crms$mut = "FMRE"
		colnames(pats_crms)[1] = "patient"
		#CDS
		z2 = which(names(mutations_in_cds) %in% pair)
		pats_cds = as.data.frame(mutations_in_cds[z2])
		pats_cds = as.data.frame(pats_cds[!duplicated(pats_cds), ])
		pats_cds$mut = "CDS"
		colnames(pats_cds)[1] = "patient"
		#combine 
		patients_wmuts = rbind(pats_crms, pats_cds)
 		patients_wmuts = patients_wmuts[!duplicated(patients_wmuts), ]
		
		#set up contigency table
		#number of pats with both muts, #fmre only, #cds only, #no muts
		#both 
		both = as.data.table(table(patients_wmuts$patient))
		both = filter(both, N ==2)
		both_pats = both$V1
		both = dim(both)[1]
		#overlap
		num_overlap = both

		#test only if have at least 1 patient in common 
		if(!(length(both_pats)==0)){
			canc_both_pats = paste(unique(patient_table$V2[patient_table$V1 %in% both_pats]), collapse="_")

			#fmre only
			fmre = unique(as.character(patients_wmuts$patient[patients_wmuts$mut == "FMRE"]))
			fmre = fmre[-(which(fmre %in% both_pats))]
			FMRE_yes = length(unique(fmre))
			
			#cds only
			cds = unique(patients_wmuts$patient[patients_wmuts$mut == "CDS"])
			cds = cds[-which(cds %in% both_pats)]
			CDS_yes = length(unique(cds))
			
			#none 
			none = length(which(!(patient_table$V1 %in% unique(patients_wmuts$patient))))
			
			#contigency table 
			cont_table = matrix(c(both, CDS_yes,FMRE_yes , none), nrow=2, byrow=T)
			colnames(cont_table) = c("FMRE_yes", "FMRE_no")
			rownames(cont_table) = c("CDS_yes", "CDS__no")
			#p <-tableGrob(cont_table)
			#print(grid.arrange(top=paste(pair[1], pair[2]), p))
			f = fisher.test(cont_table, alt = "greater")
			row = c(pair, f$p.value, f$estimate, num_overlap, canc_both_pats)
			names(row) = colnames(results_pairs)
			results_pairs = rbind(results_pairs, row)
			}
		}
	}
}

results_pairs = results_pairs[-1,]
results_pairs$fishers_pval = as.numeric(results_pairs$fishers_pval)
results_pairs$fdr = p.adjust(results_pairs$fishers_pval, method="fdr")
results_pairs = as.data.table(results_pairs)
results_pairs = results_pairs[order(fishers_pval)]

#colnames(mutations_in_crms)[12] = "reg_id"
#colnames(mutations_in_cds)[12] = "reg_id"
results_pairs$num_overlap = as.numeric(results_pairs$num_overlap)
results_pairs = results_pairs[order(-num_overlap)]
results_pairs = as.data.table(results_pairs)
results_pairs = results_pairs[order(fishers_pval)]
#write.table(results_pairs, file= "819_fmre_cds_pairs_fishers_analysis_with_lncRNAs_July12nd_KI.txt", sep="\t", quote=F, row.names=F)
fdr_sig = filter(results_pairs, fdr <= 0.05)
#write.table(fdr_sig, file= "8_fdr_sig_fmre_cds_pairs_fishers_analysis_with_lncRNAs_July12nd_KI.txt", sep="\t", quote=F, row.names=F)


#-----------------------------------------------------------------------------------------
#DON"T DO FOR NOQ DON"T NEED THIS RIGHT NOW
#-----------------------------------------------------------------------------------------
pdf("FMRE_CDS_pairs_fishers_analysis_results_w_lncRNAs_July12th.pdf", width=9, height=8)

results_pairs = filter(results_pairs, fdr <0.1)

for(i in 1:nrow(results_pairs)){
		pair = c(results_pairs$CDS_mut[i], results_pairs$FMRE_mut[i])
		#get patients that have either of these mutations or both 
		#FMRE
		z1 = which(names(mutations_in_crms) %in% pair)
		pats_crms = as.data.frame(mutations_in_crms[z1])
		pats_crms = as.data.frame(pats_crms[!duplicated(pats_crms), ])
		pats_crms$mut = "FMRE"
		colnames(pats_crms)[1] = "patient"
		#CDS
		z2 = which(names(mutations_in_cds) %in% pair)
		pats_cds = as.data.frame(mutations_in_cds[z2])
		pats_cds = as.data.frame(pats_cds[!duplicated(pats_cds), ])
		pats_cds$mut = "CDS"
		colnames(pats_cds)[1] = "patient"
		#combine 
		patients_wmuts = rbind(pats_crms, pats_cds)
 		patients_wmuts = patients_wmuts[!duplicated(patients_wmuts), ]
		
		#patients_wmuts = patients_wmuts[patients_wmuts$mut=="FMRE",]
		#saveRDS(patients_wmuts, file="TP53_FMRE_patients_that_haveit.rds")
 		#subset to pancreatic patients 
 		#panc = names(patient2cancer_type)[patient2cancer_type == "Lung-AdenoCA"]
		#patients_wmuts = subset(patients_wmuts, patients_wmuts$patient %in% panc)

		#set up contigency table
		#number of pats with both muts, #fmre only, #cds only, #no muts
		#both 
		both = as.data.table(table(patients_wmuts$patient))
		both = filter(both, N ==2)
		both_pats = both$V1
		both = dim(both)[1]
		#overlap
		num_overlap = both
		#test only if have at least 1 patient in common 
			#fmre only
			fmre = unique(patients_wmuts$patient[patients_wmuts$mut == "FMRE"])
			fmre = fmre[-(which(fmre %in% both_pats))]
			FMRE_yes = length(unique(fmre))
			
			#cds only
			cds = unique(patients_wmuts$patient[patients_wmuts$mut == "CDS"])
			cds = cds[-which(cds %in% both_pats)]
			CDS_yes = length(unique(cds))
			
			#none 
			none = length(which(!(patient_table$V1 %in% unique(patients_wmuts$patient))))
	
			#contigency table 
			cont_table = matrix(c(both, CDS_yes,FMRE_yes , none), nrow=2, byrow=T)
			colnames(cont_table) = c("FMRE_yes", "FMRE_no")
			rownames(cont_table) = c("CDS_yes", "CDS__no")
			#pdf("FMRE_CDS_Lungadeno_cancer_pairs_fishers_analysis_results_May10th.pdf", width=10)
			f = fisher.test(cont_table, alt = "greater")
			#summary how many patienst per cancer for the overlapping patietns 
			#z = which(mutations_in_cds$mut_patient %in% both_pats)
			#pp = mutations_in_cds[z,6:7]
			#pp = pp[!duplicated(pp), ]
			pp = patient_table$V2[patient_table$V1 %in% both_pats]
			pp = as.data.table(table(pp))
			pp = pp[order(N)]
			colnames(pp) = c("Cancer", "Freq")
			pp = tableGrob(pp)
			p <-tableGrob(cont_table)
			print(grid.arrange(top=paste(pair[1], pair[2]), p, pp))
			print(grid.text(paste("pvalue =", f$p.value), x = unit(0.5, "npc"), 
          	y = unit(.6, "npc")))
          	#print(grid.text(paste("fdr =", round(results_pairs$fdr[i], digits=4)), x = unit(0.5, "npc"), 
          	#y = unit(.5, "npc")))
          	#dev.off()
}
dev.off()

##check how close together the top pairs are and if they are in the same cancer types 
colnames(mutations_in_crms)[12] = "FMRE_mut"
colnames(mutations_in_cds)[12] = "CDS_mut"

results_pairs$cds_chr = ""
results_pairs$cds_start = ""
results_pairs$cds_end = ""
results_pairs$fmre_chr = ""
results_pairs$fmre_start = ""
results_pairs$fmre_end = ""
results_pairs$same_chr = ""

for(i in 1:nrow(results_pairs)){
	z = which(mutations_in_cds$CDS_mut %in% results_pairs$CDS_mut[i])
	results_pairs$cds_chr[i] = mutations_in_cds$mut_chr[z][1]
	results_pairs$cds_start[i] = mutations_in_cds$reg_starts[z][1]
	results_pairs$cds_end[i] = mutations_in_cds$reg_ends[z][1]

	z = which(mutations_in_crms$FMRE_mut %in% results_pairs$FMRE_mut[i])
	results_pairs$fmre_chr[i] = mutations_in_crms$mut_chr[z][1]
	results_pairs$fmre_start[i] = mutations_in_crms$reg_starts[z][1]
	results_pairs$fmre_end[i] = mutations_in_crms$reg_ends[z][1]

	check = results_pairs$cds_chr[i] == results_pairs$fmre_chr[i]
	results_pairs$same_chr[i] = check
}

#Essentially we need a series of fisher’s exact tests
#for all pairs (CDS_driver_X, FMRE_driver_Y). 

#write.table(results_pairs, file= "1026_fmre_cds_fishers_analysis_May10th_KI_wcancers_coordinates.txt", sep="\t", quote=F, row.names=F)


####Generate heatmap###----------------------------------------------------------------------------------------------------------------

#can you please create a heat map (geom tile) where
#> - genes are Y axis and FMREs are X axis
#> - ordered by their significance in driver analysis (yep)
#> - color is a gradient between darkcyan (weak) and orange (strong)
#> - color is -log10 FDR; FDR>0.1 is set as NA
#> - number of shared patients printed on tile

results_pairs = fread("819_fmre_cds_pairs_fishers_analysis_with_lncRNAs_July12nd_KI.txt")

library(plyr)
library(dplyr)

#keep just gene name for CDS
clean_gene = function(gene){
	#gene = results_pairs$CDS_mut[1]
	r =  unlist(strsplit(gene, "::"))[3]
	return(r)
}

results_pairs$CDS_mut = unlist(llply(results_pairs$CDS_mut, clean_gene))

clean_fmre = function(fmre){
	#fmre = results_pairs$FMRE_mut[1]
	r =  unlist(strsplit(fmre, "::"))[2]
	return(r)
}

results_pairs$FMRE_mut = unlist(llply(results_pairs$FMRE_mut, clean_fmre))
#results_pairs$FMRE_mut = as.character(results_pairs$FMRE_mut)

#1. order by significance in driver analysis 
#using fdr_element column for this 

fmres$id = llply(fmres$id, clean_fmre)
fmres = as.data.table(fmres)
fmres = fmres[order(fdr_element)]

order= as.character(fmres$id)

#relevel order 
results_pairs$FMRE_mut <- factor(results_pairs$FMRE_mut, levels = order)

#relevel pcg order
lnc_coding_drivers = fread("july12_cds_lnc_drivers.txt")
coding_drivers = lnc_coding_drivers
coding_drivers$id = unlist(llply(coding_drivers$id, clean_gene))
coding_drivers = as.data.table(coding_drivers)
coding_drivers = coding_drivers[order(-fdr_element)]
order = as.character(coding_drivers$id)

#relevel order
results_pairs$CDS_mut <- factor(results_pairs$CDS_mut, levels = order)

#results_pairs = results_pairs[match(fmres$id, results_pairs$FMRE_mut),]

#2. color is pvalue (fisher's pvalue) FDR -> -log10 FDR 
#FDR > 0.1 is NA 
results_pairs$fdr_plotting = -log10(results_pairs$fdr)
results_pairs$fdr_plotting[results_pairs$fdr >0.1] = NA

#now just need to plot heatmap 
results_pairs$CDS_mut = as.character(results_pairs$CDS_mut)
results_pairs = as.data.frame(results_pairs)

#shorten fmre name
clean_fmre = function(fmre){
	#fmre = results_pairs$FMRE_mut[1]
	r =  unlist(strsplit(fmre, "-"))[1]
	return(r)
}
results_pairs$FMRE_mut = unlist(llply(as.character(results_pairs$FMRE_mut), clean_fmre))

fmres$id = llply(fmres$id, clean_fmre)
fmres = as.data.table(fmres)
fmres = fmres[order(fdr_element)]

order= as.character(fmres$id)

#relevel order 
results_pairs$FMRE_mut <- factor(results_pairs$FMRE_mut, levels = order)
order = as.character(coding_drivers$id)
results_pairs$CDS_mut <- factor(results_pairs$CDS_mut, levels = order)

#shorten fmre name
clean_fmre = function(fmre){
	#fmre = results_pairs$FMRE_mut[1]
	r =  unlist(strsplit(fmre, "chr"))[2]
	return(r)
}
results_pairs$FMRE_mut = unlist(llply(as.character(results_pairs$FMRE_mut), clean_fmre))
fmres$id = llply(fmres$id, clean_fmre)
fmres = as.data.table(fmres)
fmres = fmres[order(fdr_element)]
order= as.character(fmres$id)
results_pairs$FMRE_mut <- factor(results_pairs$FMRE_mut, levels = order)

pdf("819_pairs_heatmap-log10fdr.pdf", width=9)

g = ggplot(results_pairs, aes(FMRE_mut, CDS_mut)) +
  geom_tile(aes(fill = fdr_plotting)) +
  geom_text(aes(label = round(num_overlap, 1)), size=1.6) +
    scale_fill_gradient(low = "darkcyan", high = "orange", na.value = 'white') +
    xlab("FMRE") + ylab("CDS gene") + theme_bw()
ggpar(g,
 font.tickslab = c(7,"plain", "black"),
 xtickslab.rt = 45)

dev.off()

pdf("819_pairs_heatmap_normal_fdr.pdf", width=9)

g = ggplot(results_pairs, aes(FMRE_mut, CDS_mut)) +
  geom_tile(aes(fill = fdr)) +
  geom_text(aes(label = round(num_overlap, 1)), size=1.6) +
    scale_fill_gradient(low = "darkcyan", high = "orange", na.value = 'white') +
    xlab("FMRE") + ylab("CDS gene") + theme_bw()
ggpar(g,
 font.tickslab = c(7,"plain", "black"),
 xtickslab.rt = 45)

dev.off()

#remove FMREs adn genes wtih no sig cells
results_pairs$fdr_tag = ""
results_pairs$fdr_tag[is.na(results_pairs$fdr_plotting)] = "no"

#num unique CDS #46
unique_cds = length(unique(results_pairs$CDS_mut))

check = as.data.table(table(results_pairs$FMRE_mut, results_pairs$fdr_tag))
check = filter(check, V2 =="no")
colnames(check) = c("fmre", "no", "numNAs")

#check how many pcg prtners each fmre has to compare
check_pcgs = as.data.table(table(results_pairs$FMRE_mut))
colnames(check_pcgs) = c("fmre", "num_pcgs")
check_pcgs = merge(check_pcgs, check, by="fmre")
check_pcgs$match = ""
check_pcgs$match[check_pcgs$num_pcgs == check_pcgs$numNAs] = "remove"
remove = check_pcgs$fmre[check_pcgs$match == "remove"]

results_pairs_rm = results_pairs[-which(results_pairs$FMRE_mut %in% remove),]

#do the same thing but with PCGs
check = as.data.table(table(results_pairs$CDS_mut, results_pairs$fdr_tag))
check = filter(check, V2 =="no")
colnames(check) = c("cds", "no", "numNAs")

#check how many pcg prtners each fmre has to compare
check_pcgs = as.data.table(table(results_pairs$CDS_mut))
colnames(check_pcgs) = c("cds", "num_fmres")
check_pcgs = merge(check_pcgs, check, by="cds")
check_pcgs$match = ""
check_pcgs$match[check_pcgs$num_fmres == check_pcgs$numNAs] = "remove"
remove = check_pcgs$cds[check_pcgs$match == "remove"]

results_pairs_rm = results_pairs_rm[-which(results_pairs_rm$CDS_mut %in% remove),]

pdf("819_pairs_heatmap-log10fdr_removed_unsig.pdf", width=9)

g = ggplot(results_pairs_rm, aes(FMRE_mut, CDS_mut)) +
  geom_tile(aes(fill = fdr_plotting)) +
  geom_text(aes(label = round(num_overlap, 1)), size=3) +
    scale_fill_gradient(low = "darkcyan", high = "orange", na.value = 'transparent') +
    xlab("FMRE") + ylab("CDS gene") + theme_bw()
ggpar(g,
 font.tickslab = c(8,"plain", "black"),
 xtickslab.rt = 45, legend.title="Fisher's FDR")

dev.off()


pdf("819_pairs_heatmap_normal_fdr_removed_unsig.pdf", width=9)

g = ggplot(results_pairs_rm, aes(FMRE_mut, CDS_mut)) +
  geom_tile(aes(fill = fdr)) +
  geom_text(aes(label = round(num_overlap, 1)), size=2.5) +
    scale_fill_gradient(low = "darkcyan", high = "orange", na.value = 'transparent') +
    xlab("FMRE") + ylab("CDS gene") + theme_bw()
ggpar(g,
 font.tickslab = c(7,"plain", "black"),
 xtickslab.rt = 45, legend.title="Fisher's FDR")

dev.off()


fdr_targets = fread("FMRE_table.txt")









