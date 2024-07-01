library(rNeighborGWAS)
library(gaston)
library(dplyr)
library(readr)

#this script is to make heatmaps of the region around significant SNPs,
#using haploid genotypes for better LD values

#note: the region consists of up to 100 closest SNPs
#this maximum could be changed in line 116 and 117

################################# settings
#you need to adjust this before running

#day number for easier identification of outputs
daynumber = "20240607"
#do not use slashes or points, since this will be used for the output folder name

#pool settings:
species = "wheat"
#barley or wheat
poolName = "YS_2019"
#barley pools:
#NFNB_pooled, PM_2015, Blr_pooled, Scald_pooled, SFNB_pooled
#wheat pools:
#YR_pooled, LR_2017, SR_2017, YLS_2017, STB_2019, YS_2019

################################# loading and preparing data
#before running check manually in folders if all required files are present

haplo = read.csv(paste0("results/",species,"/filtered_inputs/data_collection/haplo_",poolName,".csv"))
#geno = read.csv(paste0("results/",species,"/filtered_inputs/data_collection/geno_",poolName,".csv"))
#pheno = read.csv(paste0("results/",species,"/filtered_inputs/data_collection/pheno_",poolName,".csv"))
#gmap = read.csv(paste0("results/",species,"/filtered_inputs/data_collection/gmap_",poolName,".csv"))
#smap = read.csv(paste0("results/",species,"/filtered_inputs/data_collection/smap_",poolName,".csv"))
#af = read.csv(paste0("results/",species,"/filtered_inputs/data_collection/af_",poolName,".csv"))

#geno = as.matrix(geno)
#colnames(geno)<-c(1:ncol(geno))
#smap = as.matrix(smap)
#af = af[,1]

################################# extract significant SNPs
#this part is largely but not completely identical with sigSNP_LD

#importing data
GWAS_results = read.csv(paste0("results/",species,"/neiGWAS/data_collection/neiGWAS_results_",poolName,".csv"))
GWAS_results = mutate(GWAS_results, phenotype = poolName)
GWAS_results = select(GWAS_results, c(ncol(GWAS_results),1,2,3:(ncol(GWAS_results)-1)))
nSNPs = nrow(GWAS_results)

#filtering for significant SNPs
sig = {}
top_p = {}
top_dist = {}
top_dist_ID = {}
for(j in 1:nrow(GWAS_results)){
  sig_j = FALSE
  top_p_j = 0
  top_dist_j = "NA"
  top_dist_ID_j = -1
  for(i in (1:((ncol(GWAS_results)-3)/2))*2+3){
    if(GWAS_results[j,i]>-log10(0.05/nSNPs)){
      sig_j = TRUE
      if(GWAS_results[j,i]>top_p_j){
        top_p_j = GWAS_results[j,i]
        top_dist_j = unlist(strsplit(colnames(GWAS_results)[i],split="_",fixed=TRUE))[3]
        top_dist_ID_j = i
      }
    }
  }
  top_p = c(top_p, top_p_j)
  top_dist = c(top_dist, top_dist_j)
  top_dist_ID = c(top_dist_ID, top_dist_ID_j)
  if(sig_j){
    sig = c(sig,j)
  }
}
GWAS_results = mutate(GWAS_results, top_dist = top_dist)
GWAS_results = mutate(GWAS_results, top_p = top_p)
GWAS_results = mutate(GWAS_results, top_dist_ID = top_dist_ID)
GWAS_results_sig = slice(GWAS_results, sig)

GWAS_results_sig = select(GWAS_results_sig, c(1,2,3,ncol(GWAS_results_sig)-2,ncol(GWAS_results_sig)-1,ncol(GWAS_results_sig),4:ncol(GWAS_results_sig)-3))
GWAS_results_sig = arrange(GWAS_results_sig, -top_p)
sig_SNPs = select(GWAS_results_sig, c(1:5))

GWAS_results_for_LD_p_hm = GWAS_results[,c(2,3)]
for(sig_dist in GWAS_results_sig$top_dist_ID){
  GWAS_results_for_LD_p_hm = mutate(GWAS_results_for_LD_p_hm, new_column = GWAS_results[,sig_dist])
  colnames(GWAS_results_for_LD_p_hm)<- c("chr", "pos", 3:ncol(GWAS_results_for_LD_p_hm))
}

rm(top_p)
rm(top_p_j)
rm(top_dist)
rm(top_dist_j)
rm(top_dist_ID)
rm(top_dist_ID_j)
rm(sig)
rm(nSNPs)
rm(GWAS_results)
rm(GWAS_results_sig)
rm(sig_j)

################################# determine range for each sig. SNP

#choose range to include 100 closest SNPs
pos_max = {}
pos_min = {}
n_snps = {}
for(i in 1:nrow(sig_SNPs)){
  snps_i = filter(haplo[,c(1,2)], chr == sig_SNPs$chr[i])
  snps_i = mutate(snps_i, dist = abs(pos - sig_SNPs$pos[i]))
  snps_i = arrange(snps_i, dist)
  if(nrow(snps_i)>100){
    snps_i = snps_i[1:100,]
  }
  pos_max = c(pos_max, max(snps_i$pos))
  pos_min = c(pos_min, min(snps_i$pos))
  n_snps = c(n_snps, nrow(snps_i))
}
sig_SNPs_range = sig_SNPs[,1:3]
sig_SNPs_range = mutate(sig_SNPs_range, pos_min = pos_min)
sig_SNPs_range = mutate(sig_SNPs_range, pos_max = pos_max)
sig_SNPs_range = mutate(sig_SNPs_range, snps = n_snps)

rm(i)
rm(snps_i)
rm(pos_max)
rm(pos_min)
rm(n_snps)

################################# calculating LD
#this can take a while...

if(!file.exists(paste0("results/",species,"/LD"))){
  dir.create(path = paste0("results/",species,"/LD"))
}
if(!file.exists(paste0("results/",species,"/LD/region"))){
  dir.create(path = paste0("results/",species,"/LD/region"))
}
if(!file.exists(paste0("results/",species,"/LD/region/",daynumber,"_",poolName))){
  dir.create(path = paste0("results/",species,"/LD/region/",daynumber,"_",poolName))
}
if(!file.exists(paste0("results/",species,"/LD/region/",daynumber,"_",poolName,"/heatmaps"))){
  dir.create(path = paste0("results/",species,"/LD/region/",daynumber,"_",poolName,"/heatmaps"))
}
write_csv(sig_SNPs_range, file = paste0("results/",species,"/LD/region/",daynumber,"_",poolName,"/",poolName,"_sig_SNPs_region.csv"))

snps_all = mutate(haplo[,c(1,2)], SNP_ID = c(1:nrow(haplo)))
colnames(snps_all)<-c("chr","pos","SNP_ID")
nhaplo = ncol(haplo)-2

getStackedHeatmapColourmapRB <- function(red_res, blue_res) {
  colours = {}
  for(i in 0:(red_res-1)){
    for(j in 0:(blue_res-1)){
      colours = c(colours, rgb(i/(red_res-1), 0.0, j/(blue_res-1)))
    }
  }
  return(colours)
}

heatmap_cs_def = getStackedHeatmapColourmapRB(32,32)

for(s in 1:nrow(sig_SNPs_range)){
  snps_s = filter(snps_all, chr==sig_SNPs_range$chr[s]&pos>=sig_SNPs_range$pos_min[s]&pos<=sig_SNPs_range$pos_max[s])
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
  write_csv(LD_s, file = paste0("results/",species,"/LD/region/",daynumber,"_",poolName,"/",poolName,"_SNPwise_LD_around_",sig_SNPs_range$chr[s],"_",sig_SNPs_range$pos[s],".csv"))
  #
  hm_s = as.matrix(LD_s)
  colnames(hm_s)<-paste0(snps_s$chr,"_",snps_s$pos)
  rownames(hm_s)<-paste0(snps_s$chr,"_",snps_s$pos)
  png(filename = paste0("results/",species,"/LD/region/",daynumber,"_",poolName,"/heatmaps/",poolName,"_LD_heatmap_around_",sig_SNPs_range$chr[s],"_",sig_SNPs_range$pos[s],".png"),width=1600,height=800,units="px")
  par(mar = c(9.1,1.1,1.1,9.1))
  image(hm_s, axes = FALSE, x = c(1:nrow(hm_s)), y = c(1:ncol(hm_s)), xlab = "", ylab = "") + axis(side = 1, at = c(1:nrow(hm_s)),labels = rownames(hm_s), las = 2) + axis(side = 4, at = c(1:ncol(hm_s)),labels = colnames(hm_s), las = 2)
  dev.off()
  png(filename = paste0("results/",species,"/LD/region/",daynumber,"_",poolName,"/heatmaps/",poolName,"_LD_heatmap_around_",sig_SNPs_range$chr[s],"_",sig_SNPs_range$pos[s],"_alt.png"),width=1600,height=800,units="px")
  par(mar = c(9.1,1.1,1.1,9.1))
  image(hm_s, col = rgb(0.01*(0:100), 0.0, 0.0), zlim = c(0.0,1.0), axes = FALSE, x = c(1:nrow(hm_s)), y = c(1:ncol(hm_s)), xlab = "", ylab = "") + axis(side = 1, at = c(1:nrow(hm_s)),labels = rownames(hm_s), las = 2) + axis(side = 4, at = c(1:ncol(hm_s)),labels = colnames(hm_s), las = 2)
  dev.off()
  #
  pvalues_s = GWAS_results_for_LD_p_hm[,c(1,2,s+2)] #these are actually -log p values
  pvalues_s = filter(pvalues_s, chr==sig_SNPs_range$chr[s]&pos>=sig_SNPs_range$pos_min[s]&pos<=sig_SNPs_range$pos_max[s])
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
  png(filename = paste0("results/",species,"/LD/region/",daynumber,"_",poolName,"/heatmaps/",poolName,"_LD_and_p_heatmap_around_",sig_SNPs_range$chr[s],"_",sig_SNPs_range$pos[s],".png"),width=1600,height=800,units="px")
  par(mar = c(9.1,1.1,1.1,9.1))
  image(phm_s, col = heatmap_cs_def, zlim = c(0,1023), axes = FALSE, x = c(1:nrow(phm_s)), y = c(1:ncol(phm_s)), xlab = "", ylab = "") + axis(side = 1, at = c(1:nrow(phm_s)),labels = rownames(hm_s), las = 2) + axis(side = 4, at = c(1:ncol(phm_s)),labels = colnames(hm_s), las = 2)
  dev.off()
}

rm(hm_s)
rm(hm_s_r)
rm(LD_s)
rm(p_s)
rm(p_s_r)
rm(phm_s)
rm(pvalues_s)
rm(snp_i)
rm(snp_j)
rm(snps_s)
rm(counts)
rm(heatmap_cs_def)
rm(i)
rm(j)
rm(k)
rm(LD_i)
rm(nhaplo)
rm(r2)
rm(s)
rm(sig_dist)

################################# plot with selfGWAS p-values

#import selfGWAS p-values
#note: selfGWAS_results files have to be manually copied from folders into selfGWAS_compare/data_collection
GWAS_results_self = read.csv(file = paste0("results/",species,"/selfGWAS_compare/data_collection/selfGWAS_results_",poolName,".csv"))

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

heatmap_cs_def = getStackedHeatmapColourmapRGB(32,32,32)

for(s in 1:nrow(sig_SNPs_range)){
  snps_s = filter(snps_all, chr==sig_SNPs_range$chr[s]&pos>=sig_SNPs_range$pos_min[s]&pos<=sig_SNPs_range$pos_max[s])
  LD_s = read.csv(file = paste0("results/",species,"/LD/region/",daynumber,"_",poolName,"/",poolName,"_SNPwise_LD_around_",sig_SNPs_range$chr[s],"_",sig_SNPs_range$pos[s],".csv"))
  hm_s = as.matrix(LD_s)
  colnames(hm_s)<-paste0(snps_s$chr,"_",snps_s$pos)
  rownames(hm_s)<-paste0(snps_s$chr,"_",snps_s$pos)
  pvalues_nei_s = GWAS_results_for_LD_p_hm[,c(1,2,s+2)] #these are actually -log p values
  pvalues_nei_s = filter(pvalues_nei_s, chr==sig_SNPs_range$chr[s]&pos>=sig_SNPs_range$pos_min[s]&pos<=sig_SNPs_range$pos_max[s])
  pvalues_nei_s = pvalues_nei_s[order(pvalues_nei_s$pos),]
  p_nei_s = matrix(0.0, ncol = nrow(pvalues_nei_s), nrow = nrow(pvalues_nei_s))
  for(i in 1:nrow(pvalues_nei_s)){
    for(j in 1:nrow(pvalues_nei_s)){
      p_nei_s[i,j] = sqrt(pvalues_nei_s[i,3]+pvalues_nei_s[j,3]) #using sqrt for better visuality
    }
  }
  pvalues_self_s = GWAS_results_self[,c(1,2,4)] #these are actually -log p values
  pvalues_self_s = filter(pvalues_self_s, chr==sig_SNPs_range$chr[s]&pos>=sig_SNPs_range$pos_min[s]&pos<=sig_SNPs_range$pos_max[s])
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
  png(filename = paste0("results/",species,"/LD/region/",daynumber,"_",poolName,"/heatmaps/",poolName,"_LD,_self-p_and_nei-p_heatmap_around_",sig_SNPs_range$chr[s],"_",sig_SNPs_range$pos[s],".png"),width=1600,height=800,units="px")
  par(mar = c(9.1,1.1,1.1,9.1))
  image(phm_s, col = heatmap_cs_def, zlim = c(0,32767), axes = FALSE, x = c(1:nrow(phm_s)), y = c(1:ncol(phm_s)), xlab = "", ylab = "") + axis(side = 1, at = c(1:nrow(phm_s)),labels = rownames(hm_s), las = 2) + axis(side = 4, at = c(1:ncol(phm_s)),labels = colnames(hm_s), las = 2)
  dev.off()
}

rm(haplo)
rm(GWAS_results_for_LD_p_hm)
rm(GWAS_results_self)
rm(hm_s)
rm(hm_s_r)
rm(LD_s)
rm(p_nei_s)
rm(p_nei_s_r)
rm(p_self_s)
rm(p_self_s_r)
rm(phm_s)
rm(pvalues_nei_s)
rm(pvalues_self_s)
rm(snps_s)
rm(i)
rm(j)
rm(m)
rm(s)

