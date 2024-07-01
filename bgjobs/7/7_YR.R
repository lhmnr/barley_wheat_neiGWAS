library(rNeighborGWAS)
library(gaston)
library(dplyr)
library(readr)

#this script is to calculate LD for each chromosome,
#using haploid genotypes for better LD values,
#to be used for gene mapping and as an overview

#note that the LD values differ between phenos,
#as different genotypes/accessions are present in different experiments,
#which causes different allele frequencies

################################# settings
#you need to adjust this before running

#day number for easier identification of outputs
daynumber = "20240607"
#do not use slashes or points, since this will be used for the output folder name

#pool settings:
species = "wheat"
#barley or wheat
poolName = "YR_pooled"
#barley pools:
#NFNB_pooled, PM_2015, Blr_pooled, Scald_pooled, SFNB_pooled
#wheat pools:
#YR_pooled, LR_2017, SR_2017, YLS_2017, STB_2019, YS_2019

#cores to use
cores = 2L

interesting_distances = c(22,10,23)
#index of dist in (1,1.5,2,2.3,2.9,3,3.2,3.7,4,4.2,4.3,4.5,5,5.2,5.4,5.7,5.9,6,6.1,6.4,6.5,6.8,7,7.1,7.25,7.3,7.7,7.9,8,8.1,8.3,8.5,8.6,8.7,8.98,9,9.1,9.3,9.45,9.5,9.9,9.95,10)
#for YR: (22,10,23)
#for YS: (7)
#for NFNB_pooled: (9,43,3)
#for PM_2015: (25,16)
#for Blr_pooled: (19)
#for Scald_pooled: (6,4)
dist_chr_link = data.frame(dist = c("6.8","4.2","7.0"), chr = c("4A","3B","4A"))
#for YR: ("6.8","4.2","7.0") and ("4A","3B","4A")
#for YS: ("3.2") and ("2A")
#for NFNB_pooled: ("4.0","10.0","2.0") and ("chr7H","chr6H","chr5H")
#for PM_2015: ("7.25","5.7") and ("chr1H","chr6H")
#for Blr_pooled: ("6.1") and ("chr6H")
#for Scald_pooled: ("3.0","2.3") and ("chr7H","chr3H")

################################# loading and preparing data
#before running check manually in folders if all required files are present

haplo = read.csv(paste0("results/",species,"/filtered_inputs/data_collection/haplo_",poolName,".csv"))

################################# get log-p-values

#importing data
GWAS_results = read.csv(paste0("results/",species,"/neiGWAS/data_collection/neiGWAS_results_",poolName,".csv"))
GWAS_results = mutate(GWAS_results, phenotype = poolName)
GWAS_results = select(GWAS_results, c(ncol(GWAS_results),1,2,3:(ncol(GWAS_results)-1)))

interesting_distances = (interesting_distances*2.0)+3

GWAS_results_for_LD_p_hm = GWAS_results[,c(2,3)]
for(sig_dist in interesting_distances){
  GWAS_results_for_LD_p_hm = mutate(GWAS_results_for_LD_p_hm, new_column = GWAS_results[,sig_dist])
  colnames(GWAS_results_for_LD_p_hm)<- c("chr", "pos", 3:ncol(GWAS_results_for_LD_p_hm))
}
colnames(GWAS_results_for_LD_p_hm)<- c("chr", "pos", colnames(GWAS_results)[interesting_distances])
GWAS_results_for_LD_p_hm = arrange(GWAS_results_for_LD_p_hm, pos)
GWAS_results_for_LD_p_hm = arrange(GWAS_results_for_LD_p_hm, chr)
rm(interesting_distances)
rm(GWAS_results)

GWAS_results_self = read.csv(file = paste0("results/",species,"/selfGWAS_compare/data_collection/selfGWAS_results_",poolName,".csv"))
GWAS_results_self = arrange(GWAS_results_self[,c(1,2,4)], pos)
GWAS_results_self = arrange(GWAS_results_self, chr)

################################# calculating LD
#this takes long time...

if(!file.exists(paste0("results/",species,"/LD/chromosome"))){
  dir.create(path = paste0("results/",species,"/LD/chromosome"))
}
if(!file.exists(paste0("results/",species,"/LD/chromosome/",daynumber,"_",poolName))){
  dir.create(path = paste0("results/",species,"/LD/chromosome/",daynumber,"_",poolName))
}
if(!file.exists(paste0("results/",species,"/LD/chromosome/",daynumber,"_",poolName,"/heatmap"))){
  dir.create(path = paste0("results/",species,"/LD/chromosome/",daynumber,"_",poolName,"/heatmap"))
}

snps_all = mutate(haplo[,c(1,2)], SNP_ID = c(1:nrow(haplo)))
colnames(snps_all)<-c("chr","pos","SNP_ID")
nhaplo = ncol(haplo)-2

chromosomes = unique(snps_all$chr)

getStackedHeatmapColourmapRB <- function(red_res, blue_res) {
  colours = {}
  for(i in 0:(red_res-1)){
    for(j in 0:(blue_res-1)){
      colours = c(colours, rgb(i/(red_res-1), 0.0, j/(blue_res-1)))
    }
  }
  return(colours)
}
colourscheme_RB = getStackedHeatmapColourmapRB(32,32)

getStackedHeatmapColourmapRGB <- function(red_res, green_res, blue_res) {
  colours = {}
  for(i in 0:(red_res-1)){
    for(j in 0:(green_res-1)){
      for(k in 0:(blue_res-1)){
        colours = c(colours, rgb(i/(red_res-1), j/(green_res-1), k/(blue_res-1)))
      }
    }
  }
  return(colours)
}
colourscheme_RGB = getStackedHeatmapColourmapRGB(32,32,32)

library(foreach)
library(doParallel)
registerDoParallel(cores=cores) #this may not work on Windows

foreach(s = 1:length(chromosomes), .packages="dplyr") %dopar% {
  snps_s = filter(snps_all, chr==chromosomes[s])
  LD_s = data.frame(NULL)
  for(i in 1:nrow(snps_s)){
    snp_i = snps_s[i,]
    LD_i = {}
    for(j in 1:nrow(snps_s)){
      if(j<i){
        r2 = LD_s[j,i]
      }
      else{
        snp_j = snps_s[j,]
        counts = c(0,0,0)#I&J,I,J
        for(k in 3:ncol(haplo)){
          counts[1] = counts[1] + (haplo[snp_i[1,3],k]*haplo[snp_j[1,3],k])
          counts[2] = counts[2] + (haplo[snp_i[1,3],k])
          counts[3] = counts[3] + (haplo[snp_j[1,3],k])
        }
        r2 = (((nhaplo*counts[1])-(counts[2]*counts[3]))*((nhaplo*counts[1])-(counts[2]*counts[3])))/(((nhaplo*counts[2])-(counts[2]*counts[2]))*((nhaplo*counts[3])-(counts[3]*counts[3])))
      }
      LD_i = c(LD_i, r2)
    }
    if(i==1){
      LD_s = as.data.frame(t(LD_i))
    }
    else{
      LD_s = rows_append(LD_s, as.data.frame(t(LD_i)))
    }
  }
  snps_s = mutate(snps_s, chrpos = paste0(chr,"_",pos))
  colnames(LD_s)<-snps_s$chrpos
  rownames(LD_s)<-snps_s$chrpos
  write_csv(LD_s, file = paste0("results/",species,"/LD/chromosome/",daynumber,"_",poolName,"/SNPwise_LD_of_chromosome_",chromosomes[s],".csv"))
  hm_s = as.matrix(LD_s)
  png(filename = paste0("results/",species,"/LD/chromosome/",daynumber,"_",poolName,"/heatmap/LD_heatmap_of_chromosome_",chromosomes[s],".png"),width=3000,height=3000,units="px")
  par(mar = c(9.1,1.1,1.1,9.1))
  image(hm_s, axes = FALSE, x = c(1:nrow(hm_s)), y = c(1:ncol(hm_s)), xlab = "", ylab = "") + axis(side = 1, at = c(1:nrow(hm_s)),labels = rownames(hm_s), las = 2) + axis(side = 4, at = c(1:ncol(hm_s)),labels = colnames(hm_s), las = 2)
  dev.off()
  png(filename = paste0("results/",species,"/LD/chromosome/",daynumber,"_",poolName,"/heatmap/LD_heatmap_of_chromosome_",chromosomes[s],"_alt.png"),width=3000,height=3000,units="px")
  par(mar = c(9.1,1.1,1.1,9.1))
  image(hm_s, col = rgb(0.01*(0:100), 0.0, 0.0), zlim = c(0.0,1.0), axes = FALSE, x = c(1:nrow(hm_s)), y = c(1:ncol(hm_s)), xlab = "", ylab = "") + axis(side = 1, at = c(1:nrow(hm_s)),labels = rownames(hm_s), las = 2) + axis(side = 4, at = c(1:ncol(hm_s)),labels = colnames(hm_s), las = 2)
  dev.off()
  dist_chr_link_s = mutate(dist_chr_link, ID = c(1:nrow(dist_chr_link)))
  dist_chr_link_s = filter(dist_chr_link_s, chr == chromosomes[s])
  for(d in dist_chr_link_s$ID){
    pvalues_s = GWAS_results_for_LD_p_hm[,c(1,2,d+2)] #these are actually -log p values
    pvalues_s = filter(pvalues_s, chr==chromosomes[s])
    pvalues_s = pvalues_s[order(pvalues_s$pos),]
    p_s = matrix(0.0, ncol = nrow(pvalues_s), nrow = nrow(pvalues_s))
    for(i in 1:nrow(pvalues_s)){
      for(j in 1:nrow(pvalues_s)){
        p_s[i,j] = sqrt(pvalues_s[i,3]+pvalues_s[j,3]) #using sqrt for better visuality
      }
    }
    p_s = p_s/max(p_s)
    hm_s_r = round((hm_s*32)-0.5)
    p_s_r = round((p_s*32)-0.5)
    phm_s = matrix(0.0, ncol = ncol(hm_s), nrow = nrow(hm_s))
    for(i in 1:nrow(hm_s)){
      for(j in 1:ncol(hm_s)){
        if(hm_s_r[i,j]==32){hm_s_r[i,j]=31}
        if(p_s_r[i,j]==32){p_s_r[i,j]=31}
        phm_s[i,j] = (32*hm_s_r[i,j]) + p_s_r[i,j]
      }
    }
    png(filename = paste0("results/",species,"/LD/chromosome/",daynumber,"_",poolName,"/heatmap/LD_and_p_heatmap_of_chromosome_",chromosomes[s],"_and_dist_",dist_chr_link$dist[d],".png"),width=3000,height=3000,units="px")
    par(mar = c(9.1,1.1,1.1,9.1))
    image(phm_s, col = colourscheme_RB, zlim = c(0,1023), axes = FALSE, x = c(1:nrow(phm_s)), y = c(1:ncol(phm_s)), xlab = "", ylab = "") + axis(side = 1, at = c(1:nrow(phm_s)),labels = rownames(hm_s), las = 2) + axis(side = 4, at = c(1:ncol(phm_s)),labels = colnames(hm_s), las = 2)
    dev.off()
    p_nei_s = matrix(0.0, ncol = nrow(pvalues_s), nrow = nrow(pvalues_s))
    for(i in 1:nrow(pvalues_s)){
      for(j in 1:nrow(pvalues_s)){
        p_nei_s[i,j] = sqrt(pvalues_s[i,3]+pvalues_s[j,3]) #using sqrt for better visuality
      }
    }
    pvalues_self_s = filter(GWAS_results_self, chr==chromosomes[s]) #these are actually -log p values
    pvalues_self_s = pvalues_self_s[order(pvalues_self_s$pos),]
    p_self_s = matrix(0.0, ncol = nrow(pvalues_self_s), nrow = nrow(pvalues_self_s))
    for(i in 1:nrow(pvalues_self_s)){
      for(j in 1:nrow(pvalues_self_s)){
        p_self_s[i,j] = sqrt(pvalues_self_s[i,3]+pvalues_self_s[j,3]) #using sqrt for better visuality
      }
    }
    m = max(c(max(p_nei_s),max(p_self_s)))
    p_nei_s = p_nei_s/m
    p_self_s = p_self_s/m
    hm_s_r = round((hm_s*32)-0.5)
    p_nei_s_r = round((p_nei_s*32)-0.5)
    p_self_s_r = round((p_self_s*32)-0.5)
    phm_s = matrix(0.0, ncol = ncol(hm_s), nrow = nrow(hm_s))
    for(i in 1:nrow(hm_s)){
      for(j in 1:ncol(hm_s)){
        if(hm_s_r[i,j]==32){hm_s_r[i,j]=31}
        if(p_nei_s_r[i,j]==32){p_nei_s_r[i,j]=31}
        if(p_self_s_r[i,j]==32){p_self_s_r[i,j]=31}
        phm_s[i,j] = (1024*hm_s_r[i,j]) + (32*p_self_s_r[i,j]) + p_nei_s_r[i,j]
      }
    }
    png(filename = paste0("results/",species,"/LD/chromosome/",daynumber,"_",poolName,"/heatmap/LD,_self-p_and_nei-p_heatmap_of_chromosome_",chromosomes[s],"_and_dist_",dist_chr_link$dist[d],".png"),width=3000,height=3000,units="px")
    par(mar = c(9.1,1.1,1.1,9.1))
    image(phm_s, col = colourscheme_RGB, zlim = c(0,32767), axes = FALSE, x = c(1:nrow(phm_s)), y = c(1:ncol(phm_s)), xlab = "", ylab = "") + axis(side = 1, at = c(1:nrow(phm_s)),labels = rownames(hm_s), las = 2) + axis(side = 4, at = c(1:ncol(phm_s)),labels = colnames(hm_s), las = 2)
    dev.off()
  }
  paste0(chromosomes[s]," is done")
}

