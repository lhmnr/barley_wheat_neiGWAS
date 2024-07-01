library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)

#this script is to evaluate whether there are significant SNPs and get an overwiev
#and to check how p inflatedness changes over distance
#it helps to choose, which "phenos"/experiments should be used for further analysis

################################# settings

#day number for easier identification of outputs
daynumber = "20240425"

#pool settings:
species = "wheat"
#barley or wheat

################################# prepare

#you need to manually copy the neiGWAS_results csv files into data_collection
#note: all neiGWAS output datasets need to have the same distances
res_col = system(paste0("ls ./results/",species,"/neiGWAS/data_collection/*.csv"), intern=TRUE)

#importing data
GWAS_results_all = data.frame(NULL)
for(i in 1:length(res_col)){
  res_i = read.csv(res_col[i])
  filename = unlist(strsplit(res_col[i],split="/",fixed=TRUE))[6]
  filename = unlist(strsplit(filename,split=".",fixed=TRUE))[1]
  res_i = mutate(res_i, phenotype = paste0(unlist(strsplit(filename,split="_",fixed=TRUE))[3],"_",unlist(strsplit(filename,split="_",fixed=TRUE))[4]))
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
sig = {}
top_p = {}
top_dist = {}
for(j in 1:nrow(GWAS_results_all)){
  sig_j = FALSE
  top_p_j = 0
  top_dist_j = "NA"
  for(i in (1:((ncol(GWAS_results_all)-4)/2))*2+3){
    if(GWAS_results_all[j,i]>-log10(0.05/GWAS_results_all[j,ncol(GWAS_results_all)])){
      sig_j = TRUE
      if(GWAS_results_all[j,i]>top_p_j){
        top_p_j = GWAS_results_all[j,i]
        top_dist_j = unlist(strsplit(colnames(GWAS_results_all)[i],split="_",fixed=TRUE))[3]
      }
    }
  }
  top_p = c(top_p, top_p_j)
  top_dist = c(top_dist, top_dist_j)
  if(sig_j){
    sig = c(sig,j)
  }
}
GWAS_results_all = mutate(GWAS_results_all, top_dist = top_dist)
GWAS_results_all = mutate(GWAS_results_all, top_p = top_p)
GWAS_results_sig = slice(GWAS_results_all, sig)

#create output folder if necessary
if(!file.exists(paste0("results/",species,"/SNP_analysis"))){
  dir.create(path = paste0("results/",species,"/SNP_analysis"))
}
if(!file.exists(paste0("results/",species,"/SNP_analysis/sig_SNPs_vs_dist_",daynumber))){
  dir.create(path = paste0("results/",species,"/SNP_analysis/sig_SNPs_vs_dist_",daynumber))
}

#preparing
sig_SNPs = mutate(GWAS_results_sig, phenochrpos = paste0(phenotype,"_",chr,"_",pos))
cols_to_keep = c(ncol(sig_SNPs),ncol(sig_SNPs)-3,2*(1:((ncol(sig_SNPs)-6)/2))+3)
sig_SNPs = select(sig_SNPs, cols_to_keep)
for(i in 1:nrow(sig_SNPs)){
  for(j in 3:ncol(sig_SNPs)){
    sig_SNPs[i,j] = sig_SNPs[i,j]-log10(sig_SNPs[i,2])
  }
}
sig_SNPs = sig_SNPs[,c(1,3:ncol(sig_SNPs))]
sig_SNPs = gather(sig_SNPs, key = "distance", value = "p_log10", -phenochrpos)
sig_SNPs = mutate(sig_SNPs, distance = sapply(strsplit(distance,split="_",fixed=TRUE), function(x) (x[3])))

#discrete plot
x_axis_order = as.factor(unique(as.numeric(sig_SNPs$distance)))
sig_SNPs_plot = ggplot(sig_SNPs, aes(x=distance,y=p_log10,group=phenochrpos,colour=phenochrpos)) + geom_line() + geom_hline(yintercept = -log10(0.05)) + ggtitle("p-values of significant SNPs of each distance of all phenotypes") + scale_x_discrete(limits = x_axis_order) + xlab("distance of included neighbours") + ylab("bonferroni-corrected -log10 p-value")
ggsave(filename = paste0("results/",species,"/SNP_analysis/sig_SNPs_vs_dist_",daynumber,"/significant_SNPs_vs_distance.png"),plot=sig_SNPs_plot,width=4000,height=3000,units="px",device="png")
#additional threshold lines can be added with line below
sig_SNPs_plot = sig_SNPs_plot + geom_hline(yintercept = -log10(0.01)) + geom_hline(yintercept = -log10(0.001))
ggsave(filename = paste0("results/",species,"/SNP_analysis/sig_SNPs_vs_dist_",daynumber,"/significant_SNPs_vs_distance_3t.png"),plot=sig_SNPs_plot,width=4000,height=3000,units="px",device="png")
rm(sig_SNPs_plot)

#save data
GWAS_results_sig = select(GWAS_results_sig, c(1,2,3,ncol(GWAS_results_sig)-2,ncol(GWAS_results_sig)-1,ncol(GWAS_results_sig),4:ncol(GWAS_results_sig)-3))
GWAS_results_sig = arrange(GWAS_results_sig, -top_p)
write_csv(GWAS_results_sig, file = paste0("results/",species,"/SNP_analysis/sig_SNPs_vs_dist_",daynumber,"/significant_SNPs_all_dist.csv"))

#clean
rm(top_p)
rm(top_p_j)
rm(top_dist)
rm(top_dist_j)
rm(sig)
rm(res_i)
rm(cols_to_keep)
rm(sig_SNPs)
rm(filename)
rm(sig_j)
rm(x_axis_order)
rm(i)
rm(j)

################################# qq analysis

inflatedness = data.frame(NULL)
pheno = unique(GWAS_results_all$phenotype)
for(i in 1:length(pheno)){
  GWAS_results_i = filter(GWAS_results_all, phenotype == pheno[i])
  inflatedness_i = {}
  dist_i = {}
  for(j in 1:((ncol(GWAS_results_i)-6)/2)){
    GWAS_results_ij = GWAS_results_i[,(2*j)+2]
    dist_ij = unlist(strsplit(colnames(GWAS_results_i)[(2*j)+2],split="_",fixed=TRUE))[2]
    #get expected p
    nSNPs = length(GWAS_results_ij)
    expected = -log10((nSNPs:1)/(nSNPs + 1))
    observed = sort(-log10(GWAS_results_ij))
    #model slope
    slope = observed/expected
    inflatedness_ij = sum(slope)/nSNPs
    inflatedness_i = c(inflatedness_i, inflatedness_ij)
    dist_i = c(dist_i, dist_ij)
  }
  inflatedness_i = data.frame("dist"=dist_i,"x"=inflatedness_i)
  colnames(inflatedness_i)<- c("dist",pheno[i])
  if(i==1){
    inflatedness = inflatedness_i
  }
  else{
    inflatedness = merge(inflatedness, inflatedness_i)
  }
}

#save data
inflatedness = arrange(inflatedness, as.numeric(dist))
write_csv(inflatedness, file = paste0("results/",species,"/SNP_analysis/sig_SNPs_vs_dist_",daynumber,"/inflatedness_all_dist.csv"))

#preparing
inflatedness = gather(inflatedness, key = "pheno", value = "inflation", -dist)

#discrete plot
x_axis_order = as.factor(unique(as.numeric(inflatedness$dist)))
inflatedness_plot = ggplot(inflatedness, aes(x=dist,y=inflation,group=pheno,colour=pheno)) + geom_line() + ggtitle("-log10(p-value)-inflatedness vs. distance of all phenotypes") + scale_x_discrete(limits = x_axis_order) + xlab("distance of included neighbours") + ylab("-log10(p) observed/expected average")
ggsave(filename = paste0("results/",species,"/SNP_analysis/sig_SNPs_vs_dist_",daynumber,"/inflatedness_vs_distance_disc.png"),plot=inflatedness_plot,width=4000,height=3000,units="px",device="png")
rm(inflatedness_plot)

#continuous plot
inflatedness_plot = ggplot(inflatedness, aes(x=as.numeric(dist),y=inflation,group=pheno,colour=pheno)) + geom_line() + ggtitle("-log10(p-value)-inflatedness vs. distance of all phenotypes") + xlab("distance of included neighbours") + ylab("-log10(p) observed/expected average")
ggsave(filename = paste0("results/",species,"/SNP_analysis/sig_SNPs_vs_dist_",daynumber,"/inflatedness_vs_distance.png"),plot=inflatedness_plot,width=4000,height=3000,units="px",device="png")
rm(inflatedness_plot)

#clean
rm(res_col)
rm(GWAS_results_i)
rm(GWAS_results_ij)
rm(inflatedness_i)
rm(inflatedness_ij)
rm(expected)
rm(observed)
rm(slope)
rm(dist_i)
rm(dist_ij)
rm(nSNPs)
rm(pheno)
rm(i)
rm(j)
rm(x_axis_order)
detach(package:ggplot2)
detach(package:tidyr)

