#!/bin/bash

for j in `seq 1 28`; do qsub -cwd -b y -N real$j -l h_vmem=35g "module load R/3.4.0;Rscript V1_march_2019_REAL_elsatic_net_lncRNAs_main_script.R $j"; done



for j in `seq 9 9`; do qsub -cwd -b y -N real$j -l h_vmem=33g "module load R/3.4.0;Rscript V1_march_2019_REAL_elsatic_net_lncRNAs_main_script.R $j"; done

for j in `seq 22 22`; do qsub -cwd -b y -N real$j -l h_vmem=33g "module load R/3.4.0;Rscript V1_march_2019_REAL_elsatic_net_lncRNAs_main_script.R $j"; done

for j in `seq 6 6`; do qsub -cwd -b y -N real$j -l h_vmem=33g "module load R/3.4.0;Rscript V1_march_2019_REAL_elsatic_net_lncRNAs_main_script.R $j"; done
