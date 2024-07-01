library(dplyr)
library(readr)

#this script is to further analyse the transcript lists
#and fuse it with the neiGWAS results

################################# settings
#you need to adjust this before running

#day number for easier identification of outputs
daynumber = "20240530"
#do not use slashes or points, since this will be used for the output folder name

#pool settings:
species = "barley"
#barley or wheat
poolName = "PM_2015"
#barley pools:
#NFNB_pooled, PM_2015, Blr_pooled, Scald_pooled, SFNB_pooled
#wheat pools:
#YR_pooled, LR_2017, SR_2017, YLS_2017, STB_2019, YS_2019

################################# loading and preparing data

locus_info = read.csv(paste0("genelists/collection/locus_info.csv"))
sig_SNPs = read.csv(paste0("results/",species,"/sigSNP/data_collection/",poolName,"_sig_SNPs.csv"))
gwas_results = read.csv(paste0("results/",species,"/neiGWAS/data_collection/neiGWAS_results_",poolName,".csv"))

#assign loci to sig. SNPs
sig_SNP_loci = data.frame("phenotype"=character(),"chr"=character(),"locus"=character(),"SNP_pos"=numeric(),"SNP_top_dist"=numeric(),"SNP_top_p"=numeric(),"SNP_ID"=numeric())
for(i in 1:nrow(sig_SNPs)){
  chr_i = sig_SNPs[i,2]
  if(species=="barley"){
    chr_i = unlist(strsplit(chr_i,split="r",fixed=TRUE))[2]
  }
  locus_info_i = filter(locus_info, chr == chr_i)
  locus_info_i = filter(locus_info, taxa == species)
  correct_loci = {}
  for(j in 1:nrow(locus_info_i)){
    sublocus_info_j = read.csv(paste0("genelists/collection/subloci_",locus_info_i[j,1],".csv"))
    for(k in 1:nrow(sublocus_info_j)){
      if(sig_SNPs[i,3]>sublocus_info_j[k,2]&sig_SNPs[i,3]<sublocus_info_j[k,3]){
        correct_loci = c(correct_loci, locus_info_i[j,1])
      }
    }
  }
  correct_loci = unique(correct_loci)
  for(j in 1:length(correct_loci)){
    row_j = data.frame("phenotype"=sig_SNPs[i,1],"chr"=sig_SNPs[i,2],"locus"=correct_loci[j],"SNP_pos"=sig_SNPs[i,3],"SNP_top_dist"=sig_SNPs[i,4],"SNP_top_p"=sig_SNPs[i,5],"SNP_ID"=sig_SNPs[i,6])
    sig_SNP_loci = rows_append(sig_SNP_loci, row_j)
  }
}

loci = unique(sig_SNP_loci$locus)
locus_info = filter(locus_info, locus %in% loci)
loci = loci[sort.list(loci)]
gwas_results = arrange(gwas_results, pos)
gwas_results = arrange(gwas_results, chr)
gwas_results = select(gwas_results, c(1,2,((1:((ncol(gwas_results)-2)/2))*2)+1))
nei_dist = sapply(strsplit(colnames(gwas_results[,3:ncol(gwas_results)]),split="_",fixed=TRUE), function(x) (x[2]))

#clean
rm(locus_info_i)
rm(row_j)
rm(sublocus_info_j)
rm(correct_loci)
rm(chr_i)
rm(i)
rm(j)
rm(k)
rm(sig_SNPs)

#save
if(!file.exists(paste0("results/",species,"/gene_analysis"))){
  dir.create(path = paste0("results/",species,"/gene_analysis"))
}
if(!file.exists(paste0("results/",species,"/gene_analysis/",daynumber,"_",poolName))){
  dir.create(path = paste0("results/",species,"/gene_analysis/",daynumber,"_",poolName))
}
write_csv(locus_info, file = paste0("results/",species,"/gene_analysis/",daynumber,"_",poolName,"/locus_info.csv"))
write_csv(sig_SNP_loci, file = paste0("results/",species,"/gene_analysis/",daynumber,"_",poolName,"/sig_SNP_loci.csv"))

################################# locus analysis

library(ggplot2)

for(i in 1:length(loci)){
  genelist = read.csv(paste0("genelists/collection/genelist_",loci[i],".csv"))
  sig_SNPs_i = filter(sig_SNP_loci,locus==locus_info[i,1])
  gene_SNP_dist = {}
  nearest_SNP = {}
  sig_SNPs_inside = {}
  for(j in 1:nrow(genelist)){
    gene_pos_j = (genelist[j,6] + genelist[j,7])/2
    gene_SNP_dist_j = abs(sig_SNPs_i[1,4]-gene_pos_j)
    nearest_SNP_j = sig_SNPs_i[1,4]
    sig_SNPs_inside_j = 0
    for(k in 1:nrow(sig_SNPs_i)){
      if(abs(sig_SNPs_i[k,4]-gene_pos_j)<gene_SNP_dist_j){
        gene_SNP_dist_j = abs(sig_SNPs_i[k,4]-gene_pos_j)
        nearest_SNP_j = sig_SNPs_i[k,4]
      }
      if((sig_SNPs_i[k,4]<=genelist[j,7])&(sig_SNPs_i[k,4]>=genelist[j,6])){
        sig_SNPs_inside_j = sig_SNPs_inside_j + 1
      }
    }
    gene_SNP_dist = c(gene_SNP_dist, gene_SNP_dist_j)
    nearest_SNP = c(nearest_SNP, nearest_SNP_j)
    sig_SNPs_inside = c(sig_SNPs_inside, sig_SNPs_inside_j)
  }
  genelist = mutate(genelist, "nearest_sigSNP" = nearest_SNP)
  genelist = mutate(genelist, "gene_SNP_dist" = gene_SNP_dist)
  genelist = mutate(genelist, "sig_SNPs_inside" = sig_SNPs_inside)
  genelist = arrange(genelist, gene_SNP_dist)
  best_dist = sig_SNPs_i[sort.list(sig_SNPs_i[,6])[1],5]
  gwas_results_i = gwas_results[,c(1,2,which(as.numeric(nei_dist)==best_dist)+2)]
  gwas_results_i = filter(gwas_results_i, chr == sig_SNPs_i[1,2])
  p_int = {}
  for(j in 1:nrow(genelist)){
    gene_pos_j = (genelist[j,6] + genelist[j,7])/2
    gwas_results_ij_l = filter(gwas_results_i, pos <= gene_pos_j)
    gwas_results_ij_r = filter(gwas_results_i, pos >= gene_pos_j)
    SNP_j_l = gwas_results_ij_l[nrow(gwas_results_ij_l),]
    SNP_j_r = gwas_results_ij_r[1,]
    weight_j = (gene_pos_j - SNP_j_l[1,2])/(SNP_j_r[1,2] - SNP_j_l[1,2])
    p_int_j = (SNP_j_l[1,3]*(1-weight_j))+(SNP_j_r[1,3]*weight_j)
    p_int = c(p_int, p_int_j)
  }
  genelist = mutate(genelist, "interpolated_p" = p_int)
  genelist = mutate(genelist, "interpolated_p_log10" = -log10(p_int))
  weight_i = sum(1/genelist$interpolated_p)
  genelist = mutate(genelist, "score" = 1/(-(interpolated_p_log10 - log10(weight_i))))
  genelist = select(genelist, c(13,3,4,1,5,6,7,17,18,19,22,14,15,16,20,21,23,24,25,8,9,10,11,12))
  genelist = arrange(genelist, interpolated_p)
  #save
  write_csv(genelist, file = paste0("results/",species,"/gene_analysis/",daynumber,"_",poolName,"/",loci[i],"_transcript_list.csv"))
  #plot
  genelist = genelist[which(!duplicated(genelist$gene_ID)),]
  genelist = genelist[,c(6,7)]
  genelist = mutate(genelist, "position" = (transcript_pos_l+transcript_pos_r)/2)
  genelist = mutate(genelist, "type" = "gene")
  genelist = genelist[,c(3,4)]
  sig_SNPs_i = mutate(sig_SNPs_i, "type" = "sig. SNP")
  sig_SNPs_i = sig_SNPs_i[,c(4,8)]
  colnames(sig_SNPs_i)<- c("position","type")
  genelist = rows_append(genelist, sig_SNPs_i)
  gene_hist = ggplot(genelist) + geom_histogram(aes(x=position,fill=type)) + labs(title=paste0("gene density in locus ",loci[i]))
  ggsave(filename = paste0("results/",species,"/gene_analysis/",daynumber,"_",poolName,"/",loci[i],"_gene_distribution.png"),plot=gene_hist,width=4000,height=3000,units="px",device="png")
  rm(gene_hist)
}

rm(genelist)
rm(gwas_results)
rm(gwas_results_i)
rm(gwas_results_ij_l)
rm(gwas_results_ij_r)
rm(sig_SNPs_i)
rm(SNP_j_l)
rm(SNP_j_r)
rm(best_dist)
rm(gene_pos_j)
rm(gene_SNP_dist)
rm(gene_SNP_dist_j)
rm(i)
rm(j)
rm(k)
rm(loci)
rm(nearest_SNP)
rm(nearest_SNP_j)
rm(nei_dist)
rm(p_int)
rm(p_int_j)
rm(sig_SNPs_inside)
rm(sig_SNPs_inside_j)
rm(weight_i)
rm(weight_j)
detach(package:ggplot2)

