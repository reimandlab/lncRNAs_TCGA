source("/u/kisaev/lncRNAs_TCGA/NewTCGA_2017scripts/Thesis_final_scripts/Prognostic_lncs_more_detail_paper/final_scripts_2019/revisions_2020/load_data.R")

#------FEATURES-----------------------------------------------------
setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNAs_2019_manuscript")

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = subset(allCands, data == "TCGA") #173 unique lncRNA-cancer combos, #166 unique lncRNAs 
allCands$combo = unique(paste(allCands$gene, allCands$cancer, sep="_"))

#library(TCGAbiolinks)

clean_cols=c(
"Age"         ,                   "ARID1A mutation",
"ATRX status"   ,                 "CDKN2A Epigenetically Silenced",
"Chr.19.20.co.gain",              "Chr.7.gain.Chr.10.loss",
"ER Status"      ,                "Ethnicity",
"Expression subtype",            "FAT1 mutation",
"Grade"              ,            "Height",
"HER2 Status"         ,           "Histology",
"IDH Status"           ,          "Immune infiltration (ESTIMATE)",
"Karnofsky Performance Score",    "MGMT Promoter Status",
"Molecular Subtype"           ,   "MSI status",
"Mutation Count"               ,  "PAM50",
"PBRM1 Mutation"    ,             "Percent aneuploidy",
"PIK3CA Mutation"    ,            "Ploidy (ABSOLUTE)",
"PR Status"   ,                   "Purity (ESTIMATE)",
"SETD2 Mutation" ,                "Stage",
"STROMA (ESTIMATE)" ,             "Telomere Maintenance",
"TERT Promoter Status" ,          "Total Mutation Rate",
"Tumour Purity  (ABSOLUTE)",      "WHO Class",
"X1p.19q.codeletion")

#--------------------------------------------------------------------
#Clinical files - use TCGAbiolinks via previous script 
#--------------------------------------------------------------------

#--------FDR & Summarize Results-------------------------------------

all_clin = readRDS("13_data_sets_biolinks_results.rds")

fdr_sum = function(dtt){

  #first add fdr for p-values
  dtt$fdr = p.adjust(dtt$pval, method="fdr") #<- this fdr is the for the correlation test p-value
  #look at just chisq tests
  dtt$chisq = as.numeric(dtt$chisq)
  dtt$chisq_fdr = p.adjust(dtt$chisq, method="fdr")
  dtt$clin_pval_fdr = as.numeric(dtt$clin_pval)
  dtt$clin_pval_fdr = p.adjust(dtt$clin_pval, method="fdr")

  dtt$clin_vs_combo_anova_fdr = as.numeric(dtt$clin_vs_combo_anova)
  dtt$clin_vs_combo_anova_fdr = p.adjust(dtt$clin_vs_combo_anova, method="fdr")

  dtt = as.data.table(dtt)
  #dtt = filter(dtt, fdr < 0.05)

  #remove OS.time, OS, lncRNA tag... 
  z = which(dtt$colname %in% c("OS", "OS.time", "lncRNA_tag", "vital_status", 
    "Tissue.source.site", "Whole.genome", "SNP6", "HM450", "HM27", "Vital.status..1.dead.", 
    "COC", "Status", "Survival..months.", "OS.event", "Days.to.date.of.Death", "OS.event", "BCR", 
    "Tumor", "X2009stagegroup", "time_of_follow.up", "CDE_ID.3045435", "batch", "Survival", "Exome.data", "CDE_ID.3104937","OS.Time",
    "Country", "os_days", "CDE_ID.2006657", "icd_o_3_site", "WGS"))
  if(!(length(z)==0)){
    dtt = dtt[-z,]
  }

 #remove if it's just the same lncRNA correalted with itself 
 z = which(dtt$lnc == dtt$colname)
 if(!(length(z)==0)){
 dtt = dtt[-z,]}
 dtt = as.data.table(dtt)
 dtt = dtt[order(fdr)]
 print(unique(dtt$colname))
 return(dtt)

}#end fdr_sum

clean_up = llply(all_clin, fdr_sum)
clean_up = ldply(clean_up, data.table)
clean_up = as.data.table(clean_up)
clean_up = clean_up[order(fdr)]

#keep going with those associations where either spearman fdr or chisq fdr is sig 

z = which((clean_up$chisq_fdr < 0.05) | (clean_up$fdr < 0.05))
clean_up$sig_tests[z] = "*"
clean_up = as.data.table(filter(clean_up, sig_tests == "*"))

canc_conv = rna[,c("type", "Cancer")]
canc_conv = unique(canc_conv)
colnames(canc_conv)[2] = "canc"
clean_up = merge(clean_up, canc_conv, by="canc")
clean_up = clean_up[order(fdr)]
clean_up$type = factor(clean_up$type, levels=unique(clean_up$type))

clean_up$colname[which(str_detect(clean_up$colname, "age_at"))] = "Age"
clean_up$colname[clean_up$colname == "age"] = "Age"

#fdr = as.data.table(filter(clean_up, fdr < 0.05))
fdr = clean_up
fdr$cor = as.numeric(fdr$cor)
fdr$kw_pval = as.numeric(fdr$kw_pval)
#fdr = filter(fdr, (is.na(cor) | (abs(cor) >= 0.3)), kw_pval < 0.05)

write.csv(fdr, file="cleaned_clinical_variables_associations_data_sept19_precleanup.csv", quote=F, row.names=F)

fdr = clean_up
saveRDS(clean_up, file="correlation_results_clinical_lncRNA_exp_July19_using_biolinks.rds")
write.table(clean_up, file="correlation_results_clinical_lncRNA_exp_July19_using_biolnks.txt", row.names=F, quote=F)

#-------PLOT summary results-------------------------------------------

############## POST MANUAL CLEANUP ###################################################################################

#post manual cleanup of variables 
clin_results = read.csv("cleaned_clinical_variables_associations_data_sept28_post_cleanup.csv")
#clin_results$clin_pval = as.numeric(clin_results$clin_pval)
#clin_results$clin_pval_fdr = p.adjust(clin_results$clin_pval, method="fdr")

clin_results = as.data.table(clin_results)
#clin_results$clin_vs_combo_anova_fdr = p.adjust(clin_results$clin_vs_combo_anova, method="fdr")

#keep only those with significant chisq associations 
clin_results = as.data.table(filter(clin_results, chisq_fdr < 0.05 | (!(is.na(cor) & fdr < 0.05)))) #245 left

clin_results$canc_lnc_clin = paste(clin_results$combo, clin_results$colname)

dups = unique(clin_results$canc_lnc_clin[which(duplicated(clin_results$canc_lnc_clin))])

#new_dat = as.data.table(filter(clin_results, !(canc_lnc_clin %in% dups)))
new_dat = clin_results

#for(i in 1:length(dups)){
#  dup = dups[i]
#  dupdat = as.data.table(filter(clin_results, canc_lnc_clin == dup))
#  dupdat = dupdat[order(-lnc_concordance, -concordance_combo_model, clin_vs_combo_anova)]
#  keep = dupdat[1,]
#  new_dat = rbind(new_dat, keep)
#}

new_dat$combo = paste(new_dat$lnc, new_dat$canc, sep = "_")
new_dat = as.data.table(filter(new_dat, combo %in% allCands$combo))

#new_dat = 113 unique canc-lncRNA-clinical associations that are significant 
#17 unique colnames 
#unique(new_dat$colname)
# [1] "X1p.19q.codeletion"             "tumor_grade"                   
# [3] "treatment_outcome_first_course" "TERT.promoter.status"          
# [5] "TERT.expression.status"         "Spread to Lymph nodes"         
# [7] "Ethnicity "                     "PR.Status"                     
# [9] "PAM50.mRNA"                     "Original.Subtype"              
#[11] "MGMT.promoter.status"           "IDH.status"                    
#[13] "HER2.Final.Status"              "Sex"                           
#[15] "ER.Status"                      "clinical_stage"                
#[17] "Chr.7.gain.Chr.10.loss"         "Chr.19.20.co.gain"       

#10 cancer types 
length(which(new_dat$clin_pval_fdr < 0.05)) #194/245 also significnatly associated with survival 
clin_results = new_dat

#113 unique associations between a lncRNA and a clinical variable 

#which combos are better once lncRNA is used
#look at only those where clinical variable also associated with survival
clin_results$combo = paste(clin_results$name, clin_results$type)
#clinsig = as.data.table(filter(clin_results, sig_tag == "V")) #clin also sig
#clinsig_cause_lnc = as.data.table(filter(clinsig, clin_vs_combo_anova < 0.05)) #97/101
#clinsig_notcause_lnc = as.data.table(filter(clinsig, clin_vs_combo_anova > 0.05)) #4/101

z = which(str_detect(clin_results$colname, "PAM5"))
clin_results$colname[z] = "PAM50"
clin_results$name[clin_results$name == "HOXA-AS4"] = "HOXA10-AS"

#get order - 113 unique combos
t = as.data.table(table(clin_results$colname))
t = as.data.table(filter(t, N > 0))
t = t[order(-N)]

clin_results$colname = factor(clin_results$colname, levels = t$V1)

#get order 
t = as.data.table(table(clin_results$name))
t = as.data.table(filter(t, N > 0))
t = t[order(-N)]
clin_results$name = factor(clin_results$name, levels = t$V1)

#get order 
t = as.data.table(table(clin_results$type))
t = as.data.table(filter(t, N > 0))
t = t[order(-N)]

clin_results$type = factor(clin_results$type, levels = t$V1)

clin_results$concordance_combo_model = as.numeric(clin_results$concordance_combo_model)
clin_results$clin_concordance = as.numeric(clin_results$clin_concordance)

z = which((clin_results$concordance_combo_model > clin_results$clin_concordance) & (clin_results$clin_vs_combo_anova_fdr < 0.05))
clin_results$better[z] = "V"

#mark the clinical variables that are also associated with survival 
z = which(clin_results$clin_pval_fdr < 0.05)
clin_results$clin_sig[z] = "V"
clin_results$fdr_fig = ""
clin_results$fdr_fig[!(is.na(clin_results$chisq_fdr))] = clin_results$chisq_fdr[!(is.na(clin_results$chisq_fdr))]
clin_results$fdr_fig[(is.na(clin_results$chisq_fdr))] = clin_results$fdr[(is.na(clin_results$chisq_fdr))]
clin_results$fdr_fig = as.numeric(clin_results$fdr_fig)
clin_results$fdr_fig[clin_results$fdr_fig == 0] = 0.00000000001

#keep only unique combos
z = which(duplicated(clin_results$canc_lnc_clin))
dups = clin_results$canc_lnc_clin[z]
unique = as.data.table(filter(clin_results, !(canc_lnc_clin %in% dups)))
get_unique = function(combo){
 z = which(clin_results$canc_lnc_clin == combo)
 dat = clin_results[z,]  
 dat = dat[order(-concordance_combo_model)]
 return(dat[1,])
}

dups_dat = ldply(llply(dups, get_unique))
clin_results = rbind(dups_dat, unique)
length(unique(clin_results$canc_lnc_clin))

clin_results$better[clin_results$better == "V"] = "*"

pdf("summary_biolinks_subtypes_lncRNA_exp_April18.pdf", height=6, width=10)
#make geom_tile plot
ggplot(clin_results, aes(name, colname)) +
  geom_tile(aes(fill = -log10(fdr_fig), color=clin_sig, width=0.7, height=0.7), size=0.55) +
  theme_bw() + #geom_text(aes(label = better), size=2.5) + 
  theme(legend.title=element_blank(), legend.position="bottom", axis.title.x=element_blank(), 
    axis.text.x = element_text(angle = 45, hjust = 1, size=8)) +
    scale_fill_gradient(low = "tan1", high = "darkred")+
    facet_grid(cols = vars(type), scales = "free", space = "free")+
     theme(strip.background = element_rect(colour="black", fill="white", 
                                       size=1.5, linetype="solid"))+
     scale_color_manual(values=c("black", "white"))
    #+ scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

dev.off()

write.csv(clin_results, file="cleaned_clinical_variables_associations_data_sept28_post_cleanup_final_figure_data.csv", quote=F, row.names=F)

clin_results$anova_sig_combo_clin = ""
clin_results$anova_sig_combo_clin[clin_results$clin_vs_combo_anova_fdr < 0.05] = "Sig"

pdf("summary_clinical_concordances_vs_lnc_scatterplot_april18_wide.pdf", width=10, height=10)
g = ggplot(clin_results, aes(clin_concordance, concordance_combo_model, label=canc_lnc_clin)) +
  geom_point(aes(colour=type, 
       shape=anova_sig_combo_clin), fill="white", size=1.75) +
 #scale_size(range = c(0, 3))+
    #scale_colour_manual(values = mypal[c(2:5, 9,8)]) +
    #scale_fill_manual(values = sample(mypal5,9)) +  
    scale_colour_brewer(palette="Set1")+
    xlab("Clinical Concordance") + ylab("lncRNA & Clinical Combined Concordance") + theme_classic() +
    theme(legend.position = "top", axis.text = element_text(size=12), 
      legend.text=element_text(size=10), legend.title=element_text(size=10)) +
     xlim(0.5,1) + ylim(0.5,1) + geom_abline(intercept=0) + 
     geom_text_repel(data = subset(clin_results, 
      canc_lnc_clin %in% c("RP11-279F6.3 KIRP Stage", "RP5-1086K13.1 LGG X1p.19q.codeletion")),min.segment.length = unit(0, 'lines'), 
                     nudge_y = .2)
g
dev.off()

z = which((clin_results$concordance_combo_model > clin_results$clin_concordance) & (clin_results$clin_vs_combo_anova_fdr < 0.05))
clin_results$better[z] = "V"



