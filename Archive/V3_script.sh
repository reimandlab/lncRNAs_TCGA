#!/bin/bash

for i in `seq 1 100`; do for j in `seq 1 28`; do qsub -cwd -b y -N EN_geneexp$i$j -l h_vmem=35g "module load R/3.4.0;Rscript V3_march_2019_permutations_elsatic_net_lncRNAs_main_script.R $i $j"; done ; done



