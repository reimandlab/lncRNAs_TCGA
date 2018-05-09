library(data.table)
library(ggplot2)
library(ggpubr)
library(stringr)
library(dplyr)
library(gridExtra)
library(grid)

#Data---------------------------------------------------

#[1] FMREs, 41 in total
load("dfr_encode_results_signf.rsav")
head(dfr_encode_results_signf)
fmres = dfr_encode_results_signf

#[2] coding drivers - need to filter for cancer_type==“PANCANCER” and 
#element_type==“gc19_pc.cds”, 47 in total
load("all_results_signf.rsav")
head(all_results_signf)
all_results_signf = as.data.table(all_results_signf)
all_results_signf = filter(all_results_signf, cancer_type == "PANCANCER")
all_results_signf = filter(all_results_signf, element_type == "gc19_pc.cds")
coding_drivers = all_results_signf

#[3] mutations in all CRMs, subset of these are FMREs
load("encode_merge_oct2016_mutations_sanger_neutral__PANCANCER_in_elements.rsav")
head(mutations_in_elements)
mutations_in_crms = mutations_in_elements

#[4] mutations in all CDS, subset of these are CDS drivers
load("gc19_pc.cds_oct2016_mutations__PANCANCER_in_elements.rsav")
cds_mutations = mutations_in_elements

#[5] all patients in cohort
load("patient2cancertype.rsav")
head(patient2cancer_type) #1844 all together 

#[6] Clinical file 
clin <- fread("pcawg_specimen_histology_August2016_v6.tsv", data.table=F)
conversion <- fread("pcawgConversion.tsv", data.table=F)

#Analysis---------------------------------------------------

#1. Which patients have FMRE mutations  
z = which(mutations_in_crms$reg_id %in% fmres$id)
mutations_in_crms = mutations_in_crms[z,]
unique(mutations_in_crms$mut_patient)

cols = c()
for(y in 1:ncol(conversion)){
	z = which(conversion[,y] %in% mutations_in_crms$mut_patient)
	if(!(length(z)==0)){
		cols = c(cols, y)
	}
}

conversion = conversion[,c(5, cols)]
for(i in 1:nrow(mutations_in_crms)){
	pat = mutations_in_crms$mut_patient[i]
	z = which(conversion[,2] == pat)
	if(length(z)==0){
		z = which(conversion[,3] == pat)
	}
	if(length(z)==0){
		z = which(conversion[,4] == pat)
	}
	if(length(z)==0){
		z = which(conversion[,5] == pat)
	}
	if(length(z)==0){
		z = which(conversion[,6] == pat)
	}
	if(length(z)==0){
		z = which(conversion[,7] == pat)
	}
	if(length(z)==0){
		z = which(conversion[,8] == pat)
	}
	if(!(length(z)==0)){
		mutations_in_crms$mut_patient[i] = conversion[z,1]
	}
}

#keep patients that have ICGC id 
z = which(str_detect(mutations_in_crms$mut_patient, "DO")) #377 have ICGC donor ID - not sure why others don't 
mutations_in_crms = mutations_in_crms[z,]

#2. Which patients have coding mutations 
z = which(cds_mutations$reg_id %in% coding_drivers$id)
mutations_in_cds = cds_mutations[z,]

#Compare ratios------------------------------------------------

unique_cds = unique(mutations_in_cds$reg_id) #47
unique_fmre = unique(mutations_in_crms$reg_id) #40

results_pairs = as.data.frame(matrix(ncol=6)) ; colnames(results_pairs) = c("CDS_mut", "FMRE_mut", "fishers_pval", "fishers_OR", "num_overlap", "canc_overlapping_pats")

#for each cds/fmre combo
for(i in 1:length(unique_cds)){
	for(y in 1:length(unique_fmre)){
		pair = c(unique_cds[i], unique_fmre[y])
		#get patients that have either of these mutations or both 
		#FMRE
		z1 = which(mutations_in_crms$reg_id %in% pair)
		pats_crms = as.data.frame(mutations_in_crms$mut_patient[z1])
		pats_crms$mut = "FMRE"
		colnames(pats_crms)[1] = "patient"
		#CDS
		z2 = which(mutations_in_cds$reg_id %in% pair)
		pats_cds = as.data.frame(mutations_in_cds$mut_patient[z2])
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
			canc_both_pats = paste(unique(mutations_in_cds$mut_cc[which(mutations_in_cds$mut_patient == both_pats)]), collapse="_")

			#fmre only
			fmre = unique(patients_wmuts$patient[patients_wmuts$mut == "FMRE"])
			fmre = fmre[-(which(fmre %in% both_pats))]
			FMRE_yes = length(unique(fmre))
			
			#cds only
			cds = unique(patients_wmuts$patient[patients_wmuts$mut == "CDS"])
			cds = cds[-which(cds %in% both_pats)]
			CDS_yes = length(unique(cds))
			
			#none 
			none = length(which(!(names(patient2cancer_type) %in% unique(patients_wmuts$patient))))
			
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

results_pairs = results_pairs[-1,]
results_pairs$fishers_pval = as.numeric(results_pairs$fishers_pval)
results_pairs$fdr = p.adjust(results_pairs$fishers_pval, method="fdr")
results_pairs = as.data.table(results_pairs)
results_pairs = results_pairs[order(fishers_pval)]
write.table(results_pairs, file= "449_fmre_cds_fishers_analysis_May8th_KI.txt", sep="\t", quote=F, row.names=F)

colnames(mutations_in_crms)[12] = "reg_id"
colnames(mutations_in_cds)[12] = "reg_id"
results_pairs$num_overlap = as.numeric(results_pairs$num_overlap)
results_pairs = results_pairs[order(-num_overlap)]

pdf("FMRE_CDS_pairs_fishers_analysis_results_May9th.pdf", width=10)
for(i in 1:nrow(results_pairs)){
	pair = c(results_pairs$CDS_mut[i], results_pairs$FMRE_mut[i])
		z1 = which(mutations_in_crms$reg_id %in% pair)
		pats_crms = as.data.frame(mutations_in_crms$mut_patient[z1])
		pats_crms$mut = "FMRE"
		colnames(pats_crms)[1] = "patient"
		#CDS
		z2 = which(mutations_in_cds$reg_id %in% pair)
		pats_cds = as.data.frame(mutations_in_cds$mut_patient[z2])
		pats_cds$mut = "CDS"
		colnames(pats_cds)[1] = "patient"
		#combine 
		patients_wmuts = rbind(pats_crms, pats_cds)
 		patients_wmuts = patients_wmuts[!duplicated(patients_wmuts), ]
		patients_wmuts = patients_wmuts[patients_wmuts$mut=="FMRE",]
		saveRDS(patients_wmuts, file="TP53_FMRE_patients_that_haveit.rds")
 		#subset to pancreatic patients 
 		#panc = names(patient2cancer_type)[patient2cancer_type == "Ovary-AdenoCA"]
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
			none = length(which(!(names(patient2cancer_type) %in% unique(patients_wmuts$patient))))
			
			#contigency table 
			cont_table = matrix(c(both, CDS_yes,FMRE_yes , none), nrow=2, byrow=T)
			colnames(cont_table) = c("FMRE_yes", "FMRE_no")
			rownames(cont_table) = c("CDS_yes", "CDS__no")
			#pdf("FMRE_CDS_Ovarian_cancer_pairs_fishers_analysis_results_May8th.pdf", width=10)
			f = fisher.test(cont_table, alt = "greater")
			#summary how many patienst per cancer for the overlapping patietns 
			z = which(mutations_in_cds$mut_patient %in% both_pats)
			pp = mutations_in_cds[z,6:7]
			pp = pp[!duplicated(pp), ]
			pp = as.data.table(table(pp$mut_cc))
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

write.table(results_pairs, file= "449_fmre_cds_fishers_analysis_May8th_KI_wcancers_coordinates.txt", sep="\t", quote=F, row.names=F)










