#!/bin/bash

for i in `seq 1 100`; do for j in `seq 1 30`; do qsub -cwd -P reimandlab -b y -N EN_perms$i$j -l h_vmem=70g -l h_rt=500:00:00 "module load rstats;Rscript /u/kisaev/lncRNAs_TCGA/NewTCGA_2017scripts/Thesis_final_scripts/Prognostic_lncs_more_detail_paper/final_scripts_2019/revisions_2020/elastic_net_os_2021/V1_1000_runsmarch_2019_random_permutations_elsatic_net_lncRNAs_main_script.R $i $j"; done ; done
