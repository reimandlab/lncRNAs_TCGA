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

for j in `seq 1 29`; do qsub -cwd -b y -N real1000$j -l h_vmem=55g "module load rstats;Rscript V1_1000_runsmarch_2019_REAL_elsatic_net_lncRNAs_main_script.R $j"; done
