#!/bin/bash

#for j in `seq 1 28`; do qsub -cwd -b y -N real$j -l h_vmem=35g "module load R/3.4.0;Rscript V1_march_2019_REAL_elsatic_net_lncRNAs_main_script.R $j"; done
#for j in `seq 1 28`; do qsub -cwd -b y -N real$j -l h_vmem=45g "module load R/3.4.0;Rscript V1_1000_runsmarch_2019_REAL_elsatic_net_lncRNAs_main_script.R $j"; done
#for j in `seq 9 9`; do qsub -cwd -b y -N real$j -l h_vmem=33g "module load R/3.4.0;Rscript V1_march_2019_REAL_elsatic_net_lncRNAs_main_script.R $j"; done
#for j in `seq 22 22`; do qsub -cwd -b y -N real$j -l h_vmem=33g "module load R/3.4.0;Rscript V1_march_2019_REAL_elsatic_net_lncRNAs_main_script.R $j"; done
#for j in `seq 6 6`; do qsub -cwd -b y -N real$j -l h_vmem=33g "module load R/3.4.0;Rscript V1_march_2019_REAL_elsatic_net_lncRNAs_main_script.R $j"; done
#for j in `seq 20 20`; do qsub -cwd -b y -N real$j -l h_vmem=30g "module load R/3.4.0;Rscript V1_march_2019_REAL_elsatic_net_lncRNAs_main_script.R $j"; done
#for j in `seq 17 17`; do qsub -cwd -b y -N real$j -l h_vmem=30g "module load R/3.4.0;Rscript V1_march_2019_REAL_elsatic_net_lncRNAs_main_script.R $j"; done
#for j in `seq 1 28`; do qsub -cwd -b y -N real1000$j -l h_vmem=45g "module load R/3.4.0;Rscript V1_1000_runsmarch_2019_REAL_elsatic_net_lncRNAs_main_script.R $j"; done
#for j in `seq 23 23`; do qsub -cwd -b y -N real1000$j -l h_vmem=45g "module load R/3.4.0;Rscript V1_1000_runsmarch_2019_REAL_elsatic_net_lncRNAs_main_script.R $j"; done

#module load .python/3.6
for j in `seq 1 30`; do qsub -cwd -P reimandlab -b y -N real1000$j -l h_vmem=55g -l h_rt=86400 "module load rstats;Rscript /u/kisaev/lncRNAs_TCGA/NewTCGA_2017scripts/Thesis_final_scripts/Prognostic_lncs_more_detail_paper/final_scripts_2019/revisions_2020/elastic_net_os_2021/V1_1000_runsmarch_2019_REAL_elsatic_net_lncRNAs_main_script.R $j"; done
