library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)

#this script is to evaluate whether there are significant SNPs and get an overwiev
#and to check how p inflatedness changes over distance
#it helps to choose, which "phenos"/experiments should be used for further analysis

#this script is an adapted version of 3_SNP_analysis
#for asymmetric neiGWAS results

################################# settings

#day number for easier identification of outputs
daynumber = "20240513"

#pool settings:
species = "barley"
#barley or wheat

################################# prepare

#you need to manually copy the neiGWAS_results csv files into data_collection
#note: all neiGWAS output datasets need to have the same distances
res_col = system(paste0("ls ./results/",species,"/neiGWAS_asym/data_collection/*.csv"), intern=TRUE)

#importing data
GWAS_results_all = data.frame(NULL)
for(i in 1:length(res_col)){
  res_i = read.csv(res_col[i])
  filename = unlist(strsplit(res_col[i],split="/",fixed=TRUE))[6]
  filename = unlist(strsplit(filename,split=".",fixed=TRUE))[1]
  res_i = mutate(res_i, phenotype = paste0(unlist(strsplit(filename,split="_",fixed=TRUE))[4],"_",unlist(strsplit(filename,split="_",fixed=TRUE))[5]))
  res_i = select(res_i, c(ncol(res_i),1,2,3:(ncol(res_i)-1)))
  res_i = mutate(res_i, nSNPs = nrow(res_i))
  if(i==1){
    GWAS_results_all = res_i
  }
  else{
    GWAS_results_all = rows_append(GWAS_results_all, res_i)
  }
}

#filtering for significant SNPs
sig_nei = {}
top_p_nei = {}
top_dist_nei = {}
for(j in 1:nrow(GWAS_results_all)){
  sig_j = FALSE
  top_p_j = 0
  top_dist_j = "NA"
  for(i in (1:((ncol(GWAS_results_all)-4)/9))*9-1){
    if(GWAS_results_all[j,i]>-log10(0.05/GWAS_results_all[j,ncol(GWAS_results_all)])){
      sig_j = TRUE
      if(GWAS_results_all[j,i]>top_p_j){
        top_p_j = GWAS_results_all[j,i]
        top_dist_j = unlist(strsplit(colnames(GWAS_results_all)[i],split="_",fixed=TRUE))[4]
      }
    }
  }
  top_p_nei = c(top_p_nei, top_p_j)
  top_dist_nei = c(top_dist_nei, top_dist_j)
  if(sig_j){
    sig_nei = c(sig_nei,j)
  }
}
sig_sxn = {}
top_p_sxn = {}
top_dist_sxn = {}
for(j in 1:nrow(GWAS_results_all)){
  sig_j = FALSE
  top_p_j = 0
  top_dist_j = "NA"
  for(i in (1:((ncol(GWAS_results_all)-4)/9))*9+2){
    if(GWAS_results_all[j,i]>-log10(0.05/GWAS_results_all[j,ncol(GWAS_results_all)])){
      sig_j = TRUE
      if(GWAS_results_all[j,i]>top_p_j){
        top_p_j = GWAS_results_all[j,i]
        top_dist_j = unlist(strsplit(colnames(GWAS_results_all)[i],split="_",fixed=TRUE))[4]
      }
    }
  }
  top_p_sxn = c(top_p_sxn, top_p_j)
  top_dist_sxn = c(top_dist_sxn, top_dist_j)
  if(sig_j){
    sig_sxn = c(sig_sxn,j)
  }
}
GWAS_results_all = mutate(GWAS_results_all, top_dist_nei = top_dist_nei)
GWAS_results_all = mutate(GWAS_results_all, top_p_nei = top_p_nei)
GWAS_results_all = mutate(GWAS_results_all, top_dist_sxn = top_dist_sxn)
GWAS_results_all = mutate(GWAS_results_all, top_p_sxn = top_p_sxn)
GWAS_results_sig_nei = slice(GWAS_results_all, sig_nei)
GWAS_results_sig_sxn = slice(GWAS_results_all, sig_sxn)

#create output folder if necessary
if(!file.exists(paste0("results/",species,"/SNP_analysis_asym"))){
  dir.create(path = paste0("results/",species,"/SNP_analysis_asym"))
}
if(!file.exists(paste0("results/",species,"/SNP_analysis_asym/sig_SNPs_vs_dist_",daynumber))){
  dir.create(path = paste0("results/",species,"/SNP_analysis_asym/sig_SNPs_vs_dist_",daynumber))
}

#preparing
sig_SNPs_nei = mutate(GWAS_results_sig_nei, phenochrpos = paste0(phenotype,"_",chr,"_",pos))
sig_SNPs_sxn = mutate(GWAS_results_sig_sxn, phenochrpos = paste0(phenotype,"_",chr,"_",pos))
cols_to_keep_nei = c(ncol(sig_SNPs_nei),ncol(sig_SNPs_nei)-5,9*(1:((ncol(sig_SNPs_nei)-8)/9))-1)
cols_to_keep_sxn = c(ncol(sig_SNPs_sxn),ncol(sig_SNPs_sxn)-5,9*(1:((ncol(sig_SNPs_sxn)-8)/9))+2)
sig_SNPs_nei = select(sig_SNPs_nei, cols_to_keep_nei)
sig_SNPs_sxn = select(sig_SNPs_sxn, cols_to_keep_sxn)
for(i in 1:nrow(sig_SNPs_nei)){
  for(j in 3:ncol(sig_SNPs_nei)){
    sig_SNPs_nei[i,j] = sig_SNPs_nei[i,j]-log10(sig_SNPs_nei[i,2])
  }
}
for(i in 1:nrow(sig_SNPs_sxn)){
  for(j in 3:ncol(sig_SNPs_sxn)){
    sig_SNPs_sxn[i,j] = sig_SNPs_sxn[i,j]-log10(sig_SNPs_sxn[i,2])
  }
}
sig_SNPs_nei = sig_SNPs_nei[,c(1,3:ncol(sig_SNPs_nei))]
sig_SNPs_sxn = sig_SNPs_sxn[,c(1,3:ncol(sig_SNPs_sxn))]
sig_SNPs_nei = gather(sig_SNPs_nei, key = "distance", value = "p_log10", -phenochrpos)
sig_SNPs_sxn = gather(sig_SNPs_sxn, key = "distance", value = "p_log10", -phenochrpos)
sig_SNPs_nei = mutate(sig_SNPs_nei, distance = sapply(strsplit(distance,split="_",fixed=TRUE), function(x) (x[4])))
sig_SNPs_sxn = mutate(sig_SNPs_sxn, distance = sapply(strsplit(distance,split="_",fixed=TRUE), function(x) (x[4])))

#discrete plot
x_axis_order = as.factor(unique(as.numeric(sig_SNPs_nei$distance)))
sig_SNPs_plot_nei = ggplot(sig_SNPs_nei, aes(x=distance,y=p_log10,group=phenochrpos,colour=phenochrpos)) + geom_line() + geom_hline(yintercept = -log10(0.05)) + ggtitle("neighbour parity (nei) p-values of significant SNPs of each distance of all phenotypes") + scale_x_discrete(limits = x_axis_order) + xlab("distance of included neighbours") + ylab("bonferroni-corrected -log10 p-value")
ggsave(filename = paste0("results/",species,"/SNP_analysis_asym/sig_SNPs_vs_dist_",daynumber,"/nei_significant_SNPs_vs_distance.png"),plot=sig_SNPs_plot_nei,width=4000,height=3000,units="px",device="png")
sig_SNPs_plot_sxn = ggplot(sig_SNPs_sxn, aes(x=distance,y=p_log10,group=phenochrpos,colour=phenochrpos)) + geom_line() + geom_hline(yintercept = -log10(0.05)) + ggtitle("neighbour genotype (sxn) p-values of significant SNPs of each distance of all phenotypes") + scale_x_discrete(limits = x_axis_order) + xlab("distance of included neighbours") + ylab("bonferroni-corrected -log10 p-value")
ggsave(filename = paste0("results/",species,"/SNP_analysis_asym/sig_SNPs_vs_dist_",daynumber,"/sxn_significant_SNPs_vs_distance.png"),plot=sig_SNPs_plot_sxn,width=4000,height=3000,units="px",device="png")
#additional threshold lines can be added with line below
sig_SNPs_plot_nei = sig_SNPs_plot_nei + geom_hline(yintercept = -log10(0.01)) + geom_hline(yintercept = -log10(0.001))
ggsave(filename = paste0("results/",species,"/SNP_analysis_asym/sig_SNPs_vs_dist_",daynumber,"/nei_significant_SNPs_vs_distance_3t.png"),plot=sig_SNPs_plot_nei,width=4000,height=3000,units="px",device="png")
sig_SNPs_plot_sxn = sig_SNPs_plot_sxn + geom_hline(yintercept = -log10(0.01)) + geom_hline(yintercept = -log10(0.001))
ggsave(filename = paste0("results/",species,"/SNP_analysis_asym/sig_SNPs_vs_dist_",daynumber,"/sxn_significant_SNPs_vs_distance_3t.png"),plot=sig_SNPs_plot_sxn,width=4000,height=3000,units="px",device="png")
rm(sig_SNPs_plot_nei)
rm(sig_SNPs_plot_sxn)

#clean
rm(top_p_nei)
rm(top_p_sxn)

#save data
GWAS_results_sig_nei = select(GWAS_results_sig_nei, c(1,2,3,ncol(GWAS_results_sig_nei)-3,ncol(GWAS_results_sig_nei)-2,ncol(GWAS_results_sig_nei)-4,4:(ncol(GWAS_results_sig_nei)-5)))
GWAS_results_sig_sxn = select(GWAS_results_sig_sxn, c(1,2,3,ncol(GWAS_results_sig_sxn)-1,ncol(GWAS_results_sig_sxn),ncol(GWAS_results_sig_sxn)-4,4:(ncol(GWAS_results_sig_sxn)-5)))
GWAS_results_sig_nei = arrange(GWAS_results_sig_nei, -top_p_nei)
GWAS_results_sig_sxn = arrange(GWAS_results_sig_sxn, -top_p_sxn)
write_csv(GWAS_results_sig_nei, file = paste0("results/",species,"/SNP_analysis_asym/sig_SNPs_vs_dist_",daynumber,"/nei_significant_SNPs_all_dist.csv"))
write_csv(GWAS_results_sig_sxn, file = paste0("results/",species,"/SNP_analysis_asym/sig_SNPs_vs_dist_",daynumber,"/sxn_significant_SNPs_all_dist.csv"))

#clean
rm(top_p_j)
rm(top_dist_nei)
rm(top_dist_sxn)
rm(top_dist_j)
rm(sig_nei)
rm(sig_sxn)
rm(res_i)
rm(cols_to_keep_nei)
rm(cols_to_keep_sxn)
rm(sig_SNPs_nei)
rm(sig_SNPs_sxn)
rm(filename)
rm(sig_j)
rm(x_axis_order)
rm(i)
rm(j)

################################# qq analysis

inflatedness_nei = data.frame(NULL)
inflatedness_sxn = data.frame(NULL)
pheno = unique(GWAS_results_all$phenotype)
for(i in 1:length(pheno)){
  GWAS_results_i = filter(GWAS_results_all, phenotype == pheno[i])
  inflatedness_nei_i = {}
  inflatedness_sxn_i = {}
  dist_i = {}
  for(j in 1:((ncol(GWAS_results_i)-8)/9)){
    GWAS_results_ij_nei = GWAS_results_i[,(9*j)-2]
    GWAS_results_ij_sxn = GWAS_results_i[,(9*j)+1]
    dist_ij = unlist(strsplit(colnames(GWAS_results_i)[(9*j)-2],split="_",fixed=TRUE))[3]
    #get expected p
    nSNPs = length(GWAS_results_ij_nei)
    expected = -log10((nSNPs:1)/(nSNPs + 1))
    observed_nei = sort(-log10(GWAS_results_ij_nei))
    observed_sxn = sort(-log10(GWAS_results_ij_sxn))
    #model slope
    slope_nei = observed_nei/expected
    slope_sxn = observed_sxn/expected
    inflatedness_nei_ij = sum(slope_nei)/nSNPs
    inflatedness_sxn_ij = sum(slope_sxn)/nSNPs
    inflatedness_nei_i = c(inflatedness_nei_i, inflatedness_nei_ij)
    inflatedness_sxn_i = c(inflatedness_sxn_i, inflatedness_sxn_ij)
    dist_i = c(dist_i, dist_ij)
  }
  inflatedness_nei_i = data.frame("dist"=dist_i,"x"=inflatedness_nei_i)
  inflatedness_sxn_i = data.frame("dist"=dist_i,"x"=inflatedness_sxn_i)
  colnames(inflatedness_nei_i)<- c("dist",pheno[i])
  colnames(inflatedness_sxn_i)<- c("dist",pheno[i])
  if(i==1){
    inflatedness_nei = inflatedness_nei_i
    inflatedness_sxn = inflatedness_sxn_i
  }
  else{
    inflatedness_nei = merge(inflatedness_nei, inflatedness_nei_i)
    inflatedness_sxn = merge(inflatedness_sxn, inflatedness_sxn_i)
  }
}

#save data
inflatedness_nei = arrange(inflatedness_nei, as.numeric(dist))
inflatedness_sxn = arrange(inflatedness_sxn, as.numeric(dist))
write_csv(inflatedness_nei, file = paste0("results/",species,"/SNP_analysis_asym/sig_SNPs_vs_dist_",daynumber,"/nei_inflatedness_all_dist.csv"))
write_csv(inflatedness_sxn, file = paste0("results/",species,"/SNP_analysis_asym/sig_SNPs_vs_dist_",daynumber,"/sxn_inflatedness_all_dist.csv"))

#preparing
inflatedness_nei = gather(inflatedness_nei, key = "pheno", value = "inflation", -dist)
inflatedness_sxn = gather(inflatedness_sxn, key = "pheno", value = "inflation", -dist)

#discrete plot
x_axis_order = as.factor(unique(as.numeric(inflatedness_nei$dist)))
inflatedness_nei_plot = ggplot(inflatedness_nei, aes(x=dist,y=inflation,group=pheno,colour=pheno)) + geom_line() + ggtitle("neighbour parity (nei) -log10(p-value)-inflatedness vs. distance of all phenotypes") + scale_x_discrete(limits = x_axis_order) + xlab("distance of included neighbours") + ylab("-log10(p) observed/expected average")
inflatedness_sxn_plot = ggplot(inflatedness_sxn, aes(x=dist,y=inflation,group=pheno,colour=pheno)) + geom_line() + ggtitle("neighbour genotype (sxn) -log10(p-value)-inflatedness vs. distance of all phenotypes") + scale_x_discrete(limits = x_axis_order) + xlab("distance of included neighbours") + ylab("-log10(p) observed/expected average")
ggsave(filename = paste0("results/",species,"/SNP_analysis_asym/sig_SNPs_vs_dist_",daynumber,"/nei_inflatedness_vs_distance_disc.png"),plot=inflatedness_nei_plot,width=4000,height=3000,units="px",device="png")
ggsave(filename = paste0("results/",species,"/SNP_analysis_asym/sig_SNPs_vs_dist_",daynumber,"/sxn_inflatedness_vs_distance_disc.png"),plot=inflatedness_sxn_plot,width=4000,height=3000,units="px",device="png")
rm(inflatedness_nei_plot)
rm(inflatedness_sxn_plot)

#continuous plot
inflatedness_nei_plot = ggplot(inflatedness_nei, aes(x=as.numeric(dist),y=inflation,group=pheno,colour=pheno)) + geom_line() + ggtitle("neighbour parity (nei) -log10(p-value)-inflatedness vs. distance of all phenotypes") + xlab("distance of included neighbours") + ylab("-log10(p) observed/expected average")
inflatedness_sxn_plot = ggplot(inflatedness_sxn, aes(x=as.numeric(dist),y=inflation,group=pheno,colour=pheno)) + geom_line() + ggtitle("neighbour genotype (sxn) -log10(p-value)-inflatedness vs. distance of all phenotypes") + xlab("distance of included neighbours") + ylab("-log10(p) observed/expected average")
ggsave(filename = paste0("results/",species,"/SNP_analysis_asym/sig_SNPs_vs_dist_",daynumber,"/nei_inflatedness_vs_distance.png"),plot=inflatedness_nei_plot,width=4000,height=3000,units="px",device="png")
ggsave(filename = paste0("results/",species,"/SNP_analysis_asym/sig_SNPs_vs_dist_",daynumber,"/sxn_inflatedness_vs_distance.png"),plot=inflatedness_sxn_plot,width=4000,height=3000,units="px",device="png")
rm(inflatedness_nei_plot)
rm(inflatedness_sxn_plot)

#clean
rm(res_col)
rm(GWAS_results_i)
rm(GWAS_results_ij_nei)
rm(GWAS_results_ij_sxn)
rm(inflatedness_nei_i)
rm(inflatedness_sxn_i)
rm(inflatedness_nei_ij)
rm(inflatedness_sxn_ij)
rm(expected)
rm(observed_nei)
rm(observed_sxn)
rm(slope_nei)
rm(slope_sxn)
rm(dist_i)
rm(dist_ij)
rm(nSNPs)
rm(pheno)
rm(i)
rm(j)
rm(x_axis_order)
detach(package:ggplot2)
detach(package:tidyr)

