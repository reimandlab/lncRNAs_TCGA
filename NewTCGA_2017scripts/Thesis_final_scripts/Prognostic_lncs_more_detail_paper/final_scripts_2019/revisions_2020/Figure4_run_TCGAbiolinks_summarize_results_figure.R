source("/u/kisaev/lncRNAs_TCGA/NewTCGA_2017scripts/Thesis_final_scripts/Prognostic_lncs_more_detail_paper/final_scripts_2019/revisions_2020/load_data.R")

#------FEATURES-----------------------------------------------------
setwd("/.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNAs_2019_manuscript")

allCands = readRDS("final_candidates_TCGA_PCAWG_results_100CVsofElasticNet_June15.rds")
allCands = subset(allCands, data == "TCGA") #173 unique lncRNA-cancer combos, #166 unique lncRNAs 
allCands$combo = unique(paste(allCands$gene, allCands$cancer, sep="_"))

#--------------------------------------------------------------------
#Clinical files - use TCGAbiolinks via previous script 
#--------------------------------------------------------------------

#--------FDR & Summarize Results-------------------------------------

#1. adjust pvalues and remove unimportant/uninformative columns 
all_clin = readRDS("12_data_sets_biolinks_results.rds")
all_clin_dat = as.data.table(ldply(all_clin)) #combine all dataframes into one 

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
    "Country", "os_days", "CDE_ID.2006657", "icd_o_3_site", "WGS", "Pan.Glioma.RNA.Expression.Cluster", 
    "Pan.Glioma.DNA.Methylation.Cluster", 
    "Supervised.DNA.Methylation.Cluster",
    "IDH.specific.RNA.Expression.Cluster",
    "Random.Forest.Sturm.Cluster",
    "Telomere.length.estimate.in.blood.normal..Kb.",
    "Telomere.length.estimate.in.tumor..Kb.",
    "Oncogene.Negative.or.Positive.Groups",    
    "Transversion.High.Low",
    "iCluster.Group",
    "IDH.codel.subtype",    
    "Age..years.at.diagnosis.",
    "IDH.codel.subtype",
    "IDH.specific.DNA.Methylation.Cluster",
    "mRNA_K4",
    "MethyLevel",
    "miRNA.cluster",
    "SCNA.cluster",
    "protein.cluster",
    "OncoSign",
    "batch",
    "tumor_type.KIRP.path.",
    "CDE_ID.3203106",
    "CDE_ID.3203222",
    "Largest.Dimension",
    "Middle.Dimension",
    "Smallest.Dimension",
    "DNA.copy.cluster..Murray.",
    "Meth.Cluster..Laird.group.",
    "miRNA.clusters..4.group.NMF..Robertson.group.",
    "mRNA.clusters..3.group.NMF..Rathmell.group.",
    "cluster.of.clusters..CoCA...k.4..Katie.Hoadley.",
    "NF2.mutation.1",
    "Country","Days.to.Last.Follow.up","Days.to.Death",
    "Recurrence",
    "Days.to.Recurrence",
    "sequenced",
    "microRNA", "SNP6",
    "WGS", "RNAseq", "Agilent.expression", "methylation",
    "normal_meth", "complete", "Methylation.Cluster", 
    "Copy.Number.Cluster", 
    "MicroRNA.Expression.Cluster", 
    "Gene.Expression.Cluster", 
    "days_to_last_followup",
    "days_to_last_known_alive",
    "number_of_first_degree_relatives_with_cancer_diagnosis",     
    "number_of_lymphnodes_examined",
    "person_neoplasm_cancer_status",
  "vital_status",     "RNA",
  "Methylation",
  "miRNA",
  "Copy.Number", 
  "PARADIGM",
  "Survival.Data.Form",
  "Days.to.Date.of.Last.Contact",
  "Days.to.date.of.Death",
  "OS.event",
  "OS.Time",
  "SigClust.Unsupervised.mRNA",
  "SigClust.Intrinsic.mRNA",
  "miRNA.Clusters",
  "methylation.Clusters",
  "CN.Clusters",
  "Integrated.Clusters..with.PAM50.",
  "Integrated.Clusters..no.exp.",
  "Integrated.Clusters..unsup.exp.",
  "HM27","mRNA_cluster", 
  "HM450"))

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

#2. keep going with those associations where either spearman fdr or chisq fdr is sig 

z = which((clean_up$chisq_fdr < 0.05) | (clean_up$fdr < 0.05))
clean_up$sig_tests[z] = "*"
clean_up = as.data.table(filter(clean_up, sig_tests == "*")) 

#3. add short cancer name 
canc_conv = rna[,c("type", "Cancer")]
canc_conv = unique(canc_conv)
colnames(canc_conv)[2] = "canc"
clean_up = merge(clean_up, canc_conv, by="canc")
clean_up = clean_up[order(fdr)]
clean_up$type = factor(clean_up$type, levels=unique(clean_up$type))
clean_up$cor = as.numeric(clean_up$cor)
clean_up$kw_pval = as.numeric(clean_up$kw_pval)

#keep only those with significant chisq and spearman test associations 
clean_up = as.data.table(filter(clean_up, chisq_fdr < 0.05 | (!(is.na(cor)) & fdr < 0.05 & abs(cor) > 0.25))) #297 left
clean_up$combo = paste(clean_up$lnc, clean_up$canc, sep = "_")
clean_up$canc_lnc_clin = paste(clean_up$combo, clean_up$colname)

#save data set to plot results for individual clinical variables just the significant ones 
saveRDS(clean_up, file="process_tcga_biolinks_results_for_plotting.rds")

#-------PLOT summary results-------------------------------------------

#8 cancer types 
length(which(clean_up$clin_pval_fdr < 0.05)) #191/234 also significnatly associated with survival 

clean_up$colname[which(str_detect(clean_up$colname, "age_at"))] = "Age"
#clean_up$colname[clean_up$colname == "age"] = "Age"
#clean_up$colname[clean_up$colname == "Age..years.at.diagnosis."] = "Age"
clean_up$colname[which(str_detect(clean_up$colname, "Tumor.Stage..Clinical"))] = "Cliincal Stage"
clean_up$colname[which(str_detect(clean_up$colname, "clinical_stage"))] = "Cliincal Stage"
clean_up$colname[which(str_detect(clean_up$colname, "histological_grade"))] = "Grade"

clean_up$colname[which(clean_up$colname == "expression_subtype")] = "Expression.Subtype"

z = which(clean_up$colname == "Vital.Status")
clean_up = clean_up[-z,]

z = which(clean_up$colname == "PBRM1")
clean_up = clean_up[-z,]

#which combos are better once lncRNA is used
#look at only those where clinical variable also associated with survival
clean_up$name = sapply(clean_up$lnc, get_name)
clean_up$combo = paste(clean_up$name, clean_up$type)

#get order 
clean_up$name[clean_up$name == "HOXA-AS4"] = "HOXA10-AS"
t = as.data.table(table(clean_up$name))
t = as.data.table(filter(t, N > 0))
t = t[order(-N)]
clean_up$name = factor(clean_up$name, levels = t$V1)

#get order 
t = as.data.table(table(clean_up$type))
t = as.data.table(filter(t, N > 0))
t = t[order(-N)]

z = which(str_detect(clean_up$colname, "PAM5"))
#clean_up$colname[z] = "PAM50"

clean_up$type = factor(clean_up$type, levels = t$V1)

clean_up$concordance_combo_model = as.numeric(clean_up$concordance_combo_model)
clean_up$clin_concordance = as.numeric(clean_up$clin_concordance)

z = which((clean_up$concordance_combo_model > clean_up$clin_concordance) & (clean_up$clin_vs_combo_anova_fdr < 0.05))
clean_up$better = ""
clean_up$better[z] = "V"

#mark the clinical variables that are also associated with survival 
z = which(clean_up$clin_pval_fdr < 0.05)
clean_up$clin_sig[z] = "V"
clean_up$fdr_fig = ""
clean_up$fdr_fig[!(is.na(clean_up$chisq_fdr))] = clean_up$chisq_fdr[!(is.na(clean_up$chisq_fdr))]
clean_up$fdr_fig[(is.na(clean_up$chisq_fdr))] = clean_up$fdr[(is.na(clean_up$chisq_fdr))]
clean_up$fdr_fig = as.numeric(clean_up$fdr_fig)
clean_up$fdr_fig[clean_up$fdr_fig == 0] = 0.00000000001

#keep only unique combos
clean_up$canc_lnc_clin = paste(clean_up$combo, clean_up$colname)
z = which(duplicated(clean_up$canc_lnc_clin))
dups = clean_up$canc_lnc_clin[z]
unique = as.data.table(filter(clean_up, !(canc_lnc_clin %in% dups)))
get_unique = function(combo){
 z = which(clean_up$canc_lnc_clin == combo)
 dat = clean_up[z,]  
 dat = dat[order(-concordance_combo_model)]
 return(dat[1,])
}

dups_dat = ldply(llply(dups, get_unique))
clean_up = rbind(dups_dat, unique)
length(unique(clean_up$canc_lnc_clin))

clean_up$better[clean_up$better == "V"] = "*"

#get order - 113 unique combos
t = as.data.table(table(clean_up$colname))
t = as.data.table(filter(t, N > 0))
t = t[order(N)]

clean_up$colname = factor(clean_up$colname, levels = t$V1)

pdf("/u/kisaev/summary_biolinks_subtypes_lncRNA_exp_April18.pdf", height=7, width=10)
#make geom_tile plot
ggplot(clean_up, aes(name, colname)) +
  geom_tile(aes(fill = -log10(fdr_fig), color=clin_sig, width=0.7, height=0.7), size=0.55) +
  theme_bw() + #geom_text(aes(label = better), size=2.5) + 
  theme(legend.title=element_blank(), legend.position="bottom", axis.title.x=element_blank(), 
    axis.text.x = element_text(angle = 90, hjust = 1, size=6)) +
    scale_fill_gradient(low = "tan1", high = "darkred")+
    facet_grid(cols = vars(type), scales = "free", space = "free")+
     theme(strip.background = element_rect(colour="black", fill="white", 
                                       size=1.5, linetype="solid"))+
     scale_color_manual(values=c("black", "white"))
    #+ scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

dev.off()

write.csv(clean_up, file="/u/kisaev/cleaned_clinical_variables_associations_data_sept28_post_cleanup_final_figure_data.csv", quote=F, row.names=F)

clean_up$anova_sig_combo_clin = ""
clean_up$anova_sig_combo_clin[clean_up$clin_vs_combo_anova_fdr < 0.05] = "Sig"

pdf("/u/kisaev/summary_clinical_concordances_vs_lnc_scatterplot_april18_wide.pdf", width=6, height=6)
g = ggplot(clean_up, aes(clin_concordance, concordance_combo_model, label=canc_lnc_clin)) +
  geom_point(aes(colour=type, 
       shape=anova_sig_combo_clin), size=1.75) +
 #scale_size(range = c(0, 3))+
    #scale_colour_manual(values = mypal[c(2:5, 9,8)]) +
    #scale_fill_manual(values = sample(mypal5,9)) +  
    scale_colour_brewer(palette="Set1")+
    xlab("Clinical Concordance") + ylab("lncRNA & Clinical Combined Concordance") + theme_classic() +
    theme(legend.position = "top", axis.text = element_text(size=12), 
      legend.text=element_text(size=10), legend.title=element_text(size=10)) +
     xlim(0.5,1) + ylim(0.5,1) + geom_abline(intercept=0) + 
     geom_text_repel(data = subset(clean_up, 
      canc_lnc_clin %in% c("RP11-279F6.3 KIRP Cliincal Stage", "RP5-1086K13.1 LGG X1p.19q.codeletion")),min.segment.length = unit(0, 'lines'), 
                     nudge_y = .2)
g
dev.off()

z = which((clean_up$concordance_combo_model > clean_up$clin_concordance) & (clean_up$clin_vs_combo_anova_fdr < 0.05))
clean_up$better[z] = "V"



