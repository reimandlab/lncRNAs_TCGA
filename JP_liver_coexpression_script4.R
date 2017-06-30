#---------------------------------------------------------
#JP_liver_coexpression_script2.R
#---------------------------------------------------------

#Author: Karina_Isaev
#Date_started: June 28th 2017
#dir: Thesis/pcawg_liver_JP

#------------
#Description:
#------------

#Have two RNA-Seq files, one for lncRNAs (all from UCSC including low confidence) 
#and one for PCGs. Conduct here LM and NGB regression to identify 
#signficantly co-expressed lncRNA-PCGs in liver cancer patients 
#using an array of confounders 
#------------
#confounders:
#------------
#[1]. PCG CNA
#[2]. lncRNA CNA
#[3]. Clinical features 
#[4]. is methylation available?

#script4 - co-expression analysis 



