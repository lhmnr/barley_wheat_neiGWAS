library(dplyr)
library(readr)

################################# settings
#you need to adjust this before running

#day number for easier identification of outputs
daynumber = "20240513"
#do not use slashes or points, since this will be used for the output folder name

#pool settings:
species = "barley"
#barley or wheat
poolName = "NFNB_pooled"
poolName_2 = "Blr_pooled"
#barley pools:
#NFNB_pooled, PM_2015, Blr_pooled, Scald_pooled, SFNB_pooled
#wheat pools:
#YR_pooled, LR_2017, SR_2017, YLS_2017, STB_2019, YS_2019

#################################

#load and prepare
GWAS_results = read.csv(paste0("results/",species,"/neiGWAS/data_collection/neiGWAS_results_",poolName,".csv"))
GWAS_results = arrange(GWAS_results, pos)
GWAS_results = arrange(GWAS_results, chr)
GWAS_results = GWAS_results[,c(1,((1:(ncol(GWAS_results)/2))*2))]
nSNPs = nrow(GWAS_results)
distances = ncol(GWAS_results)-2
sig_SNPs = read.csv(paste0("results/",species,"/sigSNP/data_collection/",poolName_2,"_sig_SNPs.csv"))

#average -log(p) (average over nei dist)
sum(GWAS_results[,c(3:ncol(GWAS_results))])/nSNPs/distances

#-log(p) of sig SNPs (average over nei dist)
for(i in 1:nrow(sig_SNPs)){
  SNP_ID = which(GWAS_results$pos==sig_SNPs[i,3]&GWAS_results$chr==sig_SNPs[i,2])
  sum(GWAS_results[SNP_ID,c(3:ncol(GWAS_results))])/distances
}

#max -log(p) of sig SNPs
for(i in 1:nrow(sig_SNPs)){
  SNP_ID = which(GWAS_results$pos==sig_SNPs[i,3]&GWAS_results$chr==sig_SNPs[i,2])
  max(GWAS_results[SNP_ID,c(3:ncol(GWAS_results))])
}

