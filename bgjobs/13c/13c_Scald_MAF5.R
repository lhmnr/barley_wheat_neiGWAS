library(rNeighborGWAS)
library(gaston)
library(dplyr)
library(readr)

#this script is to plot effect size

#for asymmetric neiGWAS results

################################# settings
#you need to adjust this before running

#day number for easier identification of outputs
daynumber = "20240513"

#pool settings:
species = "barley"
#barley or wheat
poolName = "Scald_pooled-MAF5"
#barley pools:
#NFNB_pooled, PM_2015, Blr_pooled, Scald_pooled, SFNB_pooled
#wheat pools:
#YR_pooled, LR_2017, SR_2017, YLS_2017, STB_2019, YS_2019

################################# loading and preparing data
#before running check manually in folders if all required files are present

gwas_results = read.csv(paste0("./results/",species,"/neiGWAS_asym/data_collection/neiGWAS_results_asym_",poolName,".csv"))
nSNPs = nrow(gwas_results)

gwas_results = arrange(gwas_results, pos)
gwas_results = arrange(gwas_results, chr)

#filtering for significant SNPs
sig_self = {}
top_p_self = {}
top_dist_self = {}
top_b_self = {}
for(j in 1:nrow(gwas_results)){
  sig_j = FALSE
  top_p_j = 0
  top_dist_j = "NA"
  top_b_j = "NA"
  for(i in (1:((ncol(gwas_results)-2)/9))*9-5){
    if(gwas_results[j,i]>-log10(0.05/nSNPs)){
      sig_j = TRUE
      if(gwas_results[j,i]>top_p_j){
        top_p_j = gwas_results[j,i]
        top_dist_j = unlist(strsplit(colnames(gwas_results)[i],split="_",fixed=TRUE))[4]
        top_b_j = gwas_results[j,i+1]
      }
    }
  }
  top_p_self = c(top_p_self, top_p_j)
  top_dist_self = c(top_dist_self, top_dist_j)
  top_b_self = c(top_b_self, top_b_j)
  if(sig_j){
    sig_self = c(sig_self,j)
  }
}
#
sig_nei = {}
top_p_nei = {}
top_dist_nei = {}
top_b_nei = {}
for(j in 1:nrow(gwas_results)){
  sig_j = FALSE
  top_p_j = 0
  top_dist_j = "NA"
  top_b_j = "NA"
  for(i in (1:((ncol(gwas_results)-2)/9))*9-2){
    if(gwas_results[j,i]>-log10(0.05/nSNPs)){
      sig_j = TRUE
      if(gwas_results[j,i]>top_p_j){
        top_p_j = gwas_results[j,i]
        top_dist_j = unlist(strsplit(colnames(gwas_results)[i],split="_",fixed=TRUE))[4]
        top_b_j = gwas_results[j,i+1]
      }
    }
  }
  top_p_nei = c(top_p_nei, top_p_j)
  top_dist_nei = c(top_dist_nei, top_dist_j)
  top_b_nei = c(top_b_nei, top_b_j)
  if(sig_j){
    sig_nei = c(sig_nei,j)
  }
}
sig_sxn = {}
top_p_sxn = {}
top_dist_sxn = {}
top_b_sxn = {}
for(j in 1:nrow(gwas_results)){
  sig_j = FALSE
  top_p_j = 0
  top_dist_j = "NA"
  top_b_j = "NA"
  for(i in (1:((ncol(gwas_results)-2)/9))*9+1){
    if(gwas_results[j,i]>-log10(0.05/nSNPs)){
      sig_j = TRUE
      if(gwas_results[j,i]>top_p_j){
        top_p_j = gwas_results[j,i]
        top_dist_j = unlist(strsplit(colnames(gwas_results)[i],split="_",fixed=TRUE))[4]
        top_b_j = gwas_results[j,i+1]
      }
    }
  }
  top_p_sxn = c(top_p_sxn, top_p_j)
  top_dist_sxn = c(top_dist_sxn, top_dist_j)
  top_b_sxn = c(top_b_sxn, top_b_j)
  if(sig_j){
    sig_sxn = c(sig_sxn,j)
  }
}
gwas_results = mutate(gwas_results, top_dist_self = top_dist_self)
gwas_results = mutate(gwas_results, top_p_self = top_p_self)
gwas_results = mutate(gwas_results, top_b_self = top_b_self)
gwas_results = mutate(gwas_results, top_dist_nei = top_dist_nei)
gwas_results = mutate(gwas_results, top_p_nei = top_p_nei)
gwas_results = mutate(gwas_results, top_b_nei = top_b_nei)
gwas_results = mutate(gwas_results, top_dist_sxn = top_dist_sxn)
gwas_results = mutate(gwas_results, top_p_sxn = top_p_sxn)
gwas_results = mutate(gwas_results, top_b_sxn = top_b_sxn)
#add base info
if(species=="barley"){
  genotypes <- read.delim("./geno/barley/caigebarley__53355variants__807individuals_imp.vcf.gz", header=FALSE, comment.char="#")
}
if(species=="wheat"){
  genotypes <- read.delim("./geno/wheat/caigewheat__10511variants__801individuals_imp.vcf.gz", header=FALSE, comment.char="#")
}
genotypes = select(genotypes, c(1,2,4,5))
gwas_results = mutate(gwas_results, chr_pos = paste0(chr,"_",pos))
genotypes = mutate(genotypes, chr_pos = paste0(V1,"_",V2))
genotypes = genotypes[genotypes$chr_pos %in% gwas_results$chr_pos,3:5]
colnames(genotypes)<-c("base_ref","base_alt","chr_pos")
gwas_results = merge(gwas_results, genotypes)
gwas_results = arrange(gwas_results, pos)
gwas_results = arrange(gwas_results, chr)
#extract sig. SNPs
gwas_results_sig_self = slice(gwas_results, sig_self)
gwas_results_sig_nei = slice(gwas_results, sig_nei)
gwas_results_sig_sxn = slice(gwas_results, sig_sxn)
gwas_results_sig_self = select(gwas_results_sig_self, c(2,3,ncol(gwas_results_sig_self)-1,ncol(gwas_results_sig_self),ncol(gwas_results_sig_self)-10,ncol(gwas_results_sig_self)-9,ncol(gwas_results_sig_self)-8,4:(ncol(gwas_results_sig_self)-11)))
gwas_results_sig_nei = select(gwas_results_sig_nei, c(2,3,ncol(gwas_results_sig_nei)-1,ncol(gwas_results_sig_nei),ncol(gwas_results_sig_nei)-7,ncol(gwas_results_sig_nei)-6,ncol(gwas_results_sig_nei)-5,4:(ncol(gwas_results_sig_nei)-11)))
gwas_results_sig_sxn = select(gwas_results_sig_sxn, c(2,3,ncol(gwas_results_sig_sxn)-1,ncol(gwas_results_sig_sxn),ncol(gwas_results_sig_sxn)-4,ncol(gwas_results_sig_sxn)-3,ncol(gwas_results_sig_sxn)-2,4:(ncol(gwas_results_sig_sxn)-11)))

#create output folder if necessary
if(!file.exists(paste0("results/",species,"/sigSNP_asym"))){
  dir.create(path = paste0("results/",species,"/sigSNP_asym"))
}
if(!file.exists(paste0("results/",species,"/sigSNP_asym"))){
  dir.create(path = paste0("results/",species,"/sigSNP_asym"))
}
if(!file.exists(paste0("results/",species,"/sigSNP_asym/",daynumber,"_",poolName))){
  dir.create(path = paste0("results/",species,"/sigSNP_asym/",daynumber,"_",poolName))
}
if(!file.exists(paste0("results/",species,"/sigSNP_asym/",daynumber,"_",poolName,"/effect_vs_dist"))){
  dir.create(path = paste0("results/",species,"/sigSNP_asym/",daynumber,"_",poolName,"/effect_vs_dist"))
}
if(!file.exists(paste0("results/",species,"/sigSNP_asym/data_collection"))){
  dir.create(path = paste0("results/",species,"/sigSNP_asym/data_collection"))
}

write_csv(gwas_results_sig_self, file = paste0("results/",species,"/sigSNP_asym/",daynumber,"_",poolName,"/",poolName,"_self_sig_SNPs_p_and_effect_all_dist_with_bases.csv"))
write_csv(gwas_results_sig_nei, file = paste0("results/",species,"/sigSNP_asym/",daynumber,"_",poolName,"/",poolName,"_nei_sig_SNPs_p_and_effect_all_dist_with_bases.csv"))
write_csv(gwas_results_sig_sxn, file = paste0("results/",species,"/sigSNP_asym/",daynumber,"_",poolName,"/",poolName,"_sxn_sig_SNPs_p_and_effect_all_dist_with_bases.csv"))

sig_SNPs_nei = gwas_results_sig_nei[,1:7]
sig_SNPs_sxn = gwas_results_sig_sxn[,1:7]
sig_SNPs_nei = mutate(sig_SNPs_nei, "phenotype" = poolName)
sig_SNPs_sxn = mutate(sig_SNPs_sxn, "phenotype" = poolName)
sig_SNPs_nei = select(sig_SNPs_nei,c(8,1,2,5,6,7,3,4))
sig_SNPs_sxn = select(sig_SNPs_sxn,c(8,1,2,5,6,7,3,4))
colnames(sig_SNPs_nei)<-c("phenotype","chr","pos","top_dist","p_nei","effect_nei","base_ref","base_alt")
colnames(sig_SNPs_sxn)<-c("phenotype","chr","pos","top_dist","p_sxn","effect_sxn","base_ref","base_alt")

write_csv(sig_SNPs_nei, file = paste0("results/",species,"/sigSNP_asym/",daynumber,"_",poolName,"/",poolName,"_nei_sig_SNPs_with_bases.csv"))
write_csv(sig_SNPs_sxn, file = paste0("results/",species,"/sigSNP_asym/",daynumber,"_",poolName,"/",poolName,"_sxn_sig_SNPs_with_bases.csv"))

gwas_results_sig_nei = mutate(gwas_results_sig_nei, SNP_info = paste0(chr,"_",pos,"_(",base_ref,"_or_",base_alt,")"))
gwas_results_sig_sxn = mutate(gwas_results_sig_sxn, SNP_info = paste0(chr,"_",pos,"_(ref_",base_ref,"_alt_",base_alt,")"))
nei_sig_SNPs_p = gwas_results_sig_nei[,c(ncol(gwas_results_sig_nei),((1:((ncol(gwas_results_sig_nei)-8)/9))*9)+3)]
sxn_sig_SNPs_p = gwas_results_sig_sxn[,c(ncol(gwas_results_sig_sxn),((1:((ncol(gwas_results_sig_sxn)-8)/9))*9)+6)]
nei_sig_SNPs_effect = gwas_results_sig_nei[,c(ncol(gwas_results_sig_nei),((1:((ncol(gwas_results_sig_nei)-8)/9))*9)+4)]
sxn_sig_SNPs_effect = gwas_results_sig_sxn[,c(ncol(gwas_results_sig_sxn),((1:((ncol(gwas_results_sig_sxn)-8)/9))*9)+7)]

library(tidyr)
nei_sig_SNPs_effect = gather(nei_sig_SNPs_effect, key = "distance", value = "b", -SNP_info)
sxn_sig_SNPs_effect = gather(sxn_sig_SNPs_effect, key = "distance", value = "b", -SNP_info)
nei_sig_SNPs_p = gather(nei_sig_SNPs_p, key = "distance", value = "p_log10", -SNP_info)
sxn_sig_SNPs_p = gather(sxn_sig_SNPs_p, key = "distance", value = "p_log10", -SNP_info)
nei_sig_SNPs_effect = mutate(nei_sig_SNPs_effect, distance = sapply(strsplit(distance,split="_",fixed=TRUE), function(x) (x[3])))
sxn_sig_SNPs_effect = mutate(sxn_sig_SNPs_effect, distance = sapply(strsplit(distance,split="_",fixed=TRUE), function(x) (x[3])))
nei_sig_SNPs_p = mutate(nei_sig_SNPs_p, distance = sapply(strsplit(distance,split="_",fixed=TRUE), function(x) (x[4])))
sxn_sig_SNPs_p = mutate(sxn_sig_SNPs_p, distance = sapply(strsplit(distance,split="_",fixed=TRUE), function(x) (x[4])))
nei_sig_SNPs_effect = mutate(nei_sig_SNPs_effect, b = 2*b)
sxn_sig_SNPs_effect = mutate(sxn_sig_SNPs_effect, b = 2*b)
detach(package:tidyr)

library(ggplot2)

#discrete plot
x_axis_order = as.factor(sort(unique(as.numeric(nei_sig_SNPs_p$distance))))
if(length(x_axis_order)==0){
  x_axis_order = as.factor(sort(unique(as.numeric(sxn_sig_SNPs_p$distance))))
}
nei_sigSNPs_effect_disc = ggplot(nei_sig_SNPs_effect, aes(x=distance,y=b,group=SNP_info,colour=SNP_info)) + geom_line() + geom_hline(yintercept = 0) + ggtitle("effects of base difference (nei) in sig. SNPs on suseptibility for each distance") + scale_x_discrete(limits = x_axis_order) + xlab("distance of included neighbours") + ylab("effect on suseptibilty (1x for different base, 0,5x for heterozygous)")
nei_sigSNPs_p_disc = ggplot(nei_sig_SNPs_p, aes(x=distance,y=p_log10,group=SNP_info,colour=SNP_info)) + geom_line() + geom_hline(yintercept = -log10(0.05/nrow(gwas_results))) + ggtitle("p-values of effects of base difference (nei) in sig. SNPs on suseptibility for each distance") + scale_x_discrete(limits = x_axis_order) + xlab("distance of included neighbours") + ylab("-log10 p-value")
ggsave(filename = paste0("results/",species,"/sigSNP_asym/",daynumber,"_",poolName,"/effect_vs_dist/",poolName,"_nei_sigSNPs_effect_vs_distance_discrete.png"),plot=nei_sigSNPs_effect_disc,width=4000,height=3000,units="px",device="png")
ggsave(filename = paste0("results/",species,"/sigSNP_asym/",daynumber,"_",poolName,"/effect_vs_dist/",poolName,"_nei_sigSNPs_p_vs_distance_discrete.png"),plot=nei_sigSNPs_p_disc,width=4000,height=3000,units="px",device="png")
rm(nei_sigSNPs_effect_disc)
rm(nei_sigSNPs_p_disc)
sxn_sigSNPs_effect_disc = ggplot(sxn_sig_SNPs_effect, aes(x=distance,y=b,group=SNP_info,colour=SNP_info)) + geom_line() + geom_hline(yintercept = 0) + ggtitle("effects of alternative base (sxn) in sig. SNPs on suseptibility for each distance") + scale_x_discrete(limits = x_axis_order) + xlab("distance of included neighbours") + ylab("effect on suseptibilty (1x for alternatvie base, 0,5x for heterozygous)")
sxn_sigSNPs_p_disc = ggplot(sxn_sig_SNPs_p, aes(x=distance,y=p_log10,group=SNP_info,colour=SNP_info)) + geom_line() + geom_hline(yintercept = -log10(0.05/nrow(gwas_results))) + ggtitle("p-values of effects of alternative base (sxn) in sig. SNPs on suseptibility for each distance") + scale_x_discrete(limits = x_axis_order) + xlab("distance of included neighbours") + ylab("-log10 p-value")
ggsave(filename = paste0("results/",species,"/sigSNP_asym/",daynumber,"_",poolName,"/effect_vs_dist/",poolName,"_sxn_sigSNPs_effect_vs_distance_discrete.png"),plot=sxn_sigSNPs_effect_disc,width=4000,height=3000,units="px",device="png")
ggsave(filename = paste0("results/",species,"/sigSNP_asym/",daynumber,"_",poolName,"/effect_vs_dist/",poolName,"_sxn_sigSNPs_p_vs_distance_discrete.png"),plot=sxn_sigSNPs_p_disc,width=4000,height=3000,units="px",device="png")
rm(sxn_sigSNPs_effect_disc)
rm(sxn_sigSNPs_p_disc)

#continuous plot
nei_sig_SNPs_effect = mutate(nei_sig_SNPs_effect, distance = as.numeric(distance))
sxn_sig_SNPs_effect = mutate(sxn_sig_SNPs_effect, distance = as.numeric(distance))
nei_sig_SNPs_p = mutate(nei_sig_SNPs_p, distance = as.numeric(distance))
sxn_sig_SNPs_p = mutate(sxn_sig_SNPs_p, distance = as.numeric(distance))
nei_sigSNPs_effect_cont = ggplot(nei_sig_SNPs_effect, aes(x=distance,y=b,group=SNP_info,colour=SNP_info)) + geom_line() + geom_hline(yintercept = 0) + ggtitle("effects of base difference (nei) in sig. SNPs on suseptibility for each distance") + xlab("distance of included neighbours") + ylab("effect on suseptibilty (1x for different base, 0,5x for heterozygous)")
nei_sigSNPs_p_cont = ggplot(nei_sig_SNPs_p, aes(x=distance,y=p_log10,group=SNP_info,colour=SNP_info)) + geom_line() + geom_hline(yintercept = -log10(0.05/nrow(gwas_results))) + ggtitle("p-values of effects of base difference (nei) in sig. SNPs on suseptibility for each distance") + xlab("distance of included neighbours") + ylab("-log10 p-value")
ggsave(filename = paste0("results/",species,"/sigSNP_asym/",daynumber,"_",poolName,"/effect_vs_dist/",poolName,"_nei_sigSNPs_effect_vs_distance_continuous.png"),plot=nei_sigSNPs_effect_cont,width=4000,height=3000,units="px",device="png")
ggsave(filename = paste0("results/",species,"/sigSNP_asym/",daynumber,"_",poolName,"/effect_vs_dist/",poolName,"_nei_sigSNPs_p_vs_distance_continuous.png"),plot=nei_sigSNPs_p_cont,width=4000,height=3000,units="px",device="png")
rm(nei_sigSNPs_effect_cont)
rm(nei_sigSNPs_p_cont)
sxn_sigSNPs_effect_cont = ggplot(sxn_sig_SNPs_effect, aes(x=distance,y=b,group=SNP_info,colour=SNP_info)) + geom_line() + geom_hline(yintercept = 0) + ggtitle("effects of alternative base (sxn) in sig. SNPs on suseptibility for each distance") + xlab("distance of included neighbours") + ylab("effect on suseptibilty (1x for alternative base, 0,5x for heterozygous)")
sxn_sigSNPs_p_cont = ggplot(sxn_sig_SNPs_p, aes(x=distance,y=p_log10,group=SNP_info,colour=SNP_info)) + geom_line() + geom_hline(yintercept = -log10(0.05/nrow(gwas_results))) + ggtitle("p-values of effects of alternative base (sxn) in sig. SNPs on suseptibility for each distance") + xlab("distance of included neighbours") + ylab("-log10 p-value")
ggsave(filename = paste0("results/",species,"/sigSNP_asym/",daynumber,"_",poolName,"/effect_vs_dist/",poolName,"_sxn_sigSNPs_effect_vs_distance_continuous.png"),plot=sxn_sigSNPs_effect_cont,width=4000,height=3000,units="px",device="png")
ggsave(filename = paste0("results/",species,"/sigSNP_asym/",daynumber,"_",poolName,"/effect_vs_dist/",poolName,"_sxn_sigSNPs_p_vs_distance_continuous.png"),plot=sxn_sigSNPs_p_cont,width=4000,height=3000,units="px",device="png")
rm(sxn_sigSNPs_effect_cont)
rm(sxn_sigSNPs_p_cont)

#cleaning up
remove(sig_j)
remove(sig_self)
remove(sig_nei)
remove(sig_sxn)
remove(nei_sig_SNPs_effect)
remove(nei_sig_SNPs_p)
remove(sxn_sig_SNPs_effect)
remove(sxn_sig_SNPs_p)
remove(gwas_results)
remove(gwas_results_sig_self)
remove(gwas_results_sig_nei)
remove(gwas_results_sig_sxn)
remove(genotypes)
remove(nSNPs)
remove(i)
remove(j)
remove(top_dist_j)
remove(top_dist_self)
remove(top_dist_nei)
remove(top_dist_sxn)
remove(top_p_j)
remove(top_p_self)
remove(top_p_nei)
remove(top_p_sxn)
remove(top_b_j)
remove(top_b_self)
remove(top_b_nei)
remove(top_b_sxn)
remove(x_axis_order)
detach(package:ggplot2)

