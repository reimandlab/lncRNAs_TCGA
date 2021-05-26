#!/bin/bash

cd /.mounts/labs/reimandlab/private/users/kisaev/Thesis/TCGA_FALL2017_PROCESSED_RNASEQ/lncRNAs_2019_manuscript

qsub -cwd -P reimandlab -b y -N tcgabiolinks -l h_vmem=55g -l h_rt=86400 "module load rstats;Rscript /u/kisaev/lncRNAs_TCGA/NewTCGA_2017scripts/Thesis_final_scripts/Prognostic_lncs_more_detail_paper/final_scripts_2019/revisions_2020/Figure4_run_TCGAbiolinks_main_analysis.R"
