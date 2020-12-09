#this file contains the order of all scripts that was run in this analysis 

main=/u/kisaev/lncRNAs_TCGA/NewTCGA_2017scripts/Thesis_final_scripts/Prognostic_lncs_more_detail_paper/final_scripts_2019/revisions_2020
cd /.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ
module load rstats

#1. figure 1 related analysis

#2. elastic net analysis

#3. figure 2 related analysis

#using lncRNA candidates and clinical variables available run 1000 cross-validations using
#lncRNA only, lncRNA+clinical or clinical only models to get c-index distribution
Rscript $main/002_univariate_cross_validations_candidates.R
Rscript $main/Figure2c_figure_code.R
Rscript $main/002_multivariate_cross_validations_candidates.R
Rscript $main/002_multivariate_prepare_for_figure.R

#comparing univariate and multivariate models using lncRNAs and clinical variables
Rscript $main/002_risk_scores_multivariate_models_forrestplot.R
Rscript $main/002_risk_scores_univariate_models_forrestplot.R
Rscript $main/002_risk_scores_all_models_summary_figure.R

#comparing OS vs PFI as endpoints for survival
Rscript $main/003_survival_OS_analysis_candidates.R
Rscript $main/003_survival_PFI_analysis_candidates.R
Rscript $main/003_survival_KM_OS_vs_PFI.R

#lncRNA versus pcgs correlation and survival model comparisons
Rscript $main/004_lncRNA_pcg_correlation_analysis.R #just make KM plots for all PCAWG candidates as we did for TCGA
Rscript $main/004_lncRNA_pcg_correlation.R #just make KM plots for all PCAWG candidates as we did for TCGA

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
