#this file contains the order of all scripts that was run in this analysis

main=/u/kisaev/lncRNAs_TCGA/NewTCGA_2017scripts/Thesis_final_scripts/Prognostic_lncs_more_detail_paper/final_scripts_2019/revisions_2020
cd /.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ
module load rstats

#1. figure 1 related analysis
Rscript $main/Figure_1B_code.R
Rscript $main/001_clustering_RNAs_tSNE.R
Rscript $main/001_clustering_RNAs_umap_lgg_idh_only.R
Rscript $main/001_clustering_RNAs_umap_kirc_only_mrna_clusters.R

#2. elastic net analysis
sbatch $main/actual_elastic_net_script.sh
Rscript $main/007B_analyzing_tests_REAL_DATA_1000_en_runs.R

#3. figure 2 related analysis
#comparing OS vs PFI as endpoints for survival
Rscript $main/003_survival_OS_analysis_candidates.R
Rscript $main/003_survival_PFI_analysis_candidates.R
Rscript $main/003_survival_KM_OS_vs_PFI.R

#using lncRNA candidates and clinical variables available run 1000 cross-validations using
#lncRNA only, lncRNA+clinical or clinical only models to get c-index distribution
Rscript $main/002_univariate_cross_validations_candidates.R
#Rscript $main/Figure2c_figure_code.R
Rscript $main/002_multivariate_cross_validations_candidates.R
Rscript $main/Figure2c_figure_code.R #currently figure 2B <-use palette defined here for all cancers
Rscript $main/002_multivariate_prepare_for_figure.R
Rscript $main/002_multivariate_vs_univariate_cross_validations_candidates_figure.R
Rscript $main/002_multivariate_vs_univariate_cross_validations_candidates_figure_get_mini_plots.R
Rscript $main/002_multivariate_vs_univaraite_cross_validations_candidates_figure_local_plotting.R

#comparing univariate and multivariate models using lncRNAs and clinical variables
#Rscript $main/002_risk_scores_multivariate_models_forrestplot.R
#Rscript $main/002_risk_scores_univariate_models_forrestplot.R
#Rscript $main/002_risk_scores_all_models_summary_figure.R

#lncRNA versus pcgs correlation and survival model comparisons
Rscript $main/004_lncRNA_pcg_correlations_analysis.R
Rscript $main/004_lncRNA_pcg_correlations_medians_analysis.R
Rscript $main/004_lncRNA_pcg_correlation.R

#4. Figure 3 related analysis

#pcawg validation
Rscript $main/003_pcawg_validation_simple.R #just make KM plots for all PCAWG candidates as we did for TCGA
Rscript $main/003_pcawg_validation_risk_scores_TCGA_model.R #use TCGA model fit on lncRNA and apply to PCAWG data

#5. Figure 4 related code

#lncRNAs versus TCGAbiolinkd additional clinical variables
Rscript $main/Figure4_run_TCGAbiolinks.R

#lncRNAs versus TCGAbiolinks main analysis
Rscript $main/Figure4_run_TCGAbiolinks_main_analysis.R

#lncRNAs versus TCGAbiolinks summarize results and figure
Rscript $main/Figure4_run_TCGAbiolinks_summarize_results_figure.R

#lncRNAs versus TCGAbiolinks individual plots
Rscript $main/Figure4_run_TCGAbiolinks_individual_plots.R

#lncRNAs versus Marker genes in some specific cancers
Rscript $main/Figure4_run_TCGAbiolinks_gene_corrlations.R

#lncRNA vs XCELL immune cells
Rscript $main/Cibersort_integration_association_covariates.R

#6. Figure 6 related code

#Run LIMMA
Rscript $main/figure6_de_analysis.R
Rscript $main/figure6abcd.R
