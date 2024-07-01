library(rNeighborGWAS)
library(gaston)
library(dplyr)
library(readr)

#this script is to make heatmaps of the region around significant SNPs,
#using haploid genotypes for better LD values

#for asymmetric neiGWAS results

#note: the region consists of up to 100 closest SNPs
#this maximum could be changed in line 92 and 93

################################# settings
#you need to adjust this before running

#day number for easier identification of outputs
daynumber = "20240607"
#do not use slashes or points, since this will be used for the output folder name

#pool settings:
species = "barley"
#barley or wheat
poolName = "NFNB_pooled-MAF5"
#barley pools:
#NFNB_pooled, PM_2015, Blr_pooled, Scald_pooled, SFNB_pooled
#wheat pools:
#YR_pooled, LR_2017, SR_2017, YLS_2017, STB_2019, YS_2019

################################# loading and preparing data
#before running check manually in folders if all required files are present

#importing data
haplo = read.csv(paste0("results/",species,"/filtered_inputs/data_collection/haplo_",poolName,".csv"))

gwas_results = read.csv(paste0("results/",species,"/neiGWAS_asym/data_collection/neiGWAS_results_asym_",poolName,".csv"))
gwas_results = mutate(gwas_results, phenotype = poolName)
gwas_results = select(gwas_results, c(ncol(gwas_results),1,2,3:(ncol(gwas_results)-1)))
nSNPs = nrow(gwas_results)

sig_SNPs_nei = read.csv(paste0("results/",species,"/sigSNP_asym/data_collection/",poolName,"_nei_sig_SNPs_with_bases.csv"))
sig_SNPs_sxn = read.csv(paste0("results/",species,"/sigSNP_asym/data_collection/",poolName,"_sxn_sig_SNPs_with_bases.csv"))
sig_SNPs_nei = sig_SNPs_nei[,1:5]
sig_SNPs_sxn = sig_SNPs_sxn[,1:5]
colnames(sig_SNPs_nei)<-c("phenotype","chr","pos","top_dist","p")
colnames(sig_SNPs_sxn)<-c("phenotype","chr","pos","top_dist","p")
sig_SNPs_nei$phenotype = as.character(sig_SNPs_nei$phenotype)
sig_SNPs_sxn$phenotype = as.character(sig_SNPs_sxn$phenotype)
sig_SNPs_nei$chr = as.character(sig_SNPs_nei$chr)
sig_SNPs_sxn$chr = as.character(sig_SNPs_sxn$chr)
sig_SNPs_nei$pos = as.integer(sig_SNPs_nei$pos)
sig_SNPs_sxn$pos = as.integer(sig_SNPs_sxn$pos)
sig_SNPs_nei$top_dist = as.numeric(sig_SNPs_nei$top_dist)
sig_SNPs_sxn$top_dist = as.numeric(sig_SNPs_sxn$top_dist)
sig_SNPs_nei$p = as.numeric(sig_SNPs_nei$p)
sig_SNPs_sxn$p = as.numeric(sig_SNPs_sxn$p)
sig_SNPs_nei = mutate(sig_SNPs_nei, "type" = "nei")
sig_SNPs_sxn = mutate(sig_SNPs_sxn, "type" = "sxn")
sig_SNPs = rows_append(sig_SNPs_nei, sig_SNPs_sxn)
sig_SNPs = arrange(select(sig_SNPs, c(1,2,3,6,4,5)), -p)

gwas_results_for_LD_p_hm = gwas_results[,c(2,3,5)]
for(i in 1:nrow(sig_SNPs)){
  sig_dist = sig_SNPs$top_dist[i]
  if(sig_SNPs$type[i]=="nei"){
    gwas_results_for_LD_p_hm = mutate(gwas_results_for_LD_p_hm, new_column = gwas_results[,paste0("p_nei_log10_",sig_dist)])
    gwas_results_for_LD_p_hm = mutate(gwas_results_for_LD_p_hm, new_column_2 = gwas_results[,paste0("p_sxn_log10_",sig_dist)])
  }else{
    gwas_results_for_LD_p_hm = mutate(gwas_results_for_LD_p_hm, new_column = gwas_results[,paste0("p_sxn_log10_",sig_dist)])
    gwas_results_for_LD_p_hm = mutate(gwas_results_for_LD_p_hm, new_column_2 = gwas_results[,paste0("p_nei_log10_",sig_dist)])
  }
  colnames(gwas_results_for_LD_p_hm)<- c("chr", "pos", 3:ncol(gwas_results_for_LD_p_hm))
}

rm(sig_SNPs_nei)
rm(sig_SNPs_sxn)
rm(i)
rm(sig_dist)
rm(nSNPs)
rm(gwas_results)

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
sig_SNPs_range = sig_SNPs[,1:4]
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

if(!file.exists(paste0("results/",species,"/LD_asym"))){
  dir.create(path = paste0("results/",species,"/LD_asym"))
}
if(!file.exists(paste0("results/",species,"/LD_asym/region"))){
  dir.create(path = paste0("results/",species,"/LD_asym/region"))
}
if(!file.exists(paste0("results/",species,"/LD_asym/region/",daynumber,"_",poolName))){
  dir.create(path = paste0("results/",species,"/LD_asym/region/",daynumber,"_",poolName))
}
if(!file.exists(paste0("results/",species,"/LD_asym/region/",daynumber,"_",poolName,"/heatmaps"))){
  dir.create(path = paste0("results/",species,"/LD_asym/region/",daynumber,"_",poolName,"/heatmaps"))
}
write_csv(sig_SNPs_range, file = paste0("results/",species,"/LD_asym/region/",daynumber,"_",poolName,"/",poolName,"_sig_SNPs_region.csv"))

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
  write_csv(LD_s, file = paste0("results/",species,"/LD_asym/region/",daynumber,"_",poolName,"/",poolName,"_SNPwise_LD_around_",sig_SNPs_range$type[s],"_",sig_SNPs_range$chr[s],"_",sig_SNPs_range$pos[s],".csv"))
  #
  hm_s = as.matrix(LD_s)
  colnames(hm_s)<-paste0(snps_s$chr,"_",snps_s$pos)
  rownames(hm_s)<-paste0(snps_s$chr,"_",snps_s$pos)
  png(filename = paste0("results/",species,"/LD_asym/region/",daynumber,"_",poolName,"/heatmaps/",poolName,"_LD_heatmap_around_",sig_SNPs_range$type[s],"_",sig_SNPs_range$chr[s],"_",sig_SNPs_range$pos[s],".png"),width=1600,height=800,units="px")
  par(mar = c(9.1,1.1,1.1,9.1))
  image(hm_s, axes = FALSE, x = c(1:nrow(hm_s)), y = c(1:ncol(hm_s)), xlab = "", ylab = "") + axis(side = 1, at = c(1:nrow(hm_s)),labels = rownames(hm_s), las = 2) + axis(side = 4, at = c(1:ncol(hm_s)),labels = colnames(hm_s), las = 2)
  dev.off()
  png(filename = paste0("results/",species,"/LD_asym/region/",daynumber,"_",poolName,"/heatmaps/",poolName,"_LD_heatmap_around_",sig_SNPs_range$type[s],"_",sig_SNPs_range$chr[s],"_",sig_SNPs_range$pos[s],"_alt.png"),width=1600,height=800,units="px")
  par(mar = c(9.1,1.1,1.1,9.1))
  image(hm_s, col = rgb(0.01*(0:100), 0.0, 0.0), zlim = c(0.0,1.0), axes = FALSE, x = c(1:nrow(hm_s)), y = c(1:ncol(hm_s)), xlab = "", ylab = "") + axis(side = 1, at = c(1:nrow(hm_s)),labels = rownames(hm_s), las = 2) + axis(side = 4, at = c(1:ncol(hm_s)),labels = colnames(hm_s), las = 2)
  dev.off()
  #
  pvalues_s = gwas_results_for_LD_p_hm[,c(1,2,(s*2)+2)] #these are actually -log p values
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
  png(filename = paste0("results/",species,"/LD_asym/region/",daynumber,"_",poolName,"/heatmaps/",poolName,"_LD_and_p_heatmap_around_",sig_SNPs_range$type[s],"_",sig_SNPs_range$chr[s],"_",sig_SNPs_range$pos[s],".png"),width=1600,height=800,units="px")
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

################################# plot with p-values of the other terms

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
  LD_s = read.csv(file = paste0("results/",species,"/LD_asym/region/",daynumber,"_",poolName,"/",poolName,"_SNPwise_LD_around_",sig_SNPs_range$type[s],"_",sig_SNPs_range$chr[s],"_",sig_SNPs_range$pos[s],".csv"))
  hm_s = as.matrix(LD_s)
  colnames(hm_s)<-paste0(snps_s$chr,"_",snps_s$pos)
  rownames(hm_s)<-paste0(snps_s$chr,"_",snps_s$pos)
  pvalues_foc_s = gwas_results_for_LD_p_hm[,c(1,2,(s*2)+2)] #these are actually -log p values
  pvalues_foc_s = filter(pvalues_foc_s, chr==sig_SNPs_range$chr[s]&pos>=sig_SNPs_range$pos_min[s]&pos<=sig_SNPs_range$pos_max[s])
  pvalues_foc_s = pvalues_foc_s[order(pvalues_foc_s$pos),]
  p_foc_s = matrix(0.0, ncol = nrow(pvalues_foc_s), nrow = nrow(pvalues_foc_s))
  for(i in 1:nrow(pvalues_foc_s)){
    for(j in 1:nrow(pvalues_foc_s)){
      p_foc_s[i,j] = sqrt(pvalues_foc_s[i,3]+pvalues_foc_s[j,3]) #using sqrt for better visuality
    }
  }
  pvalues_self_s = gwas_results_for_LD_p_hm[,c(1,2,3)] #these are actually -log p values
  pvalues_self_s = filter(pvalues_self_s, chr==sig_SNPs_range$chr[s]&pos>=sig_SNPs_range$pos_min[s]&pos<=sig_SNPs_range$pos_max[s])
  pvalues_self_s = pvalues_self_s[order(pvalues_self_s$pos),]
  pvalues_alt_s = gwas_results_for_LD_p_hm[,c(1,2,(s*2)+3)] #these are actually -log p values
  pvalues_alt_s = filter(pvalues_alt_s, chr==sig_SNPs_range$chr[s]&pos>=sig_SNPs_range$pos_min[s]&pos<=sig_SNPs_range$pos_max[s])
  pvalues_alt_s = pvalues_alt_s[order(pvalues_alt_s$pos),]
  p_alt_s = matrix(0.0, ncol = nrow(pvalues_alt_s), nrow = nrow(pvalues_alt_s))
  for(i in 1:nrow(pvalues_alt_s)){
    for(j in 1:nrow(pvalues_alt_s)){
      p_alt_s[i,j] = max(sqrt(pvalues_self_s[i,3]+pvalues_self_s[j,3]),sqrt(pvalues_alt_s[i,3]+pvalues_alt_s[j,3])) #using sqrt for better visuality
    }
  }
  m = max(c(max(p_foc_s),max(p_alt_s)))
  p_foc_s = p_foc_s/m
  p_alt_s = p_alt_s/m
  hm_s_r = round((hm_s*32)-0.5)
  p_foc_s_r = round((p_foc_s*32)-0.5)
  p_alt_s_r = round((p_alt_s*32)-0.5)
  phm_s = matrix(0.0, ncol = ncol(hm_s), nrow = nrow(hm_s))
  for(i in 1:nrow(hm_s)){
    for(j in 1:ncol(hm_s)){
      if(hm_s_r[i,j]==32){hm_s_r[i,j]=31}
      if(p_foc_s_r[i,j]==32){p_foc_s_r[i,j]=31}
      if(p_alt_s_r[i,j]==32){p_alt_s_r[i,j]=31}
      phm_s[i,j] = (1024*hm_s_r[i,j]) + (32*p_alt_s_r[i,j]) + p_foc_s_r[i,j]
    }
  }
  png(filename = paste0("results/",species,"/LD_asym/region/",daynumber,"_",poolName,"/heatmaps/",poolName,"_LD,_alt-p_and_focal-p_heatmap_around_",sig_SNPs_range$type[s],"_",sig_SNPs_range$chr[s],"_",sig_SNPs_range$pos[s],".png"),width=1600,height=800,units="px")
  par(mar = c(9.1,1.1,1.1,9.1))
  image(phm_s, col = heatmap_cs_def, zlim = c(0,32767), axes = FALSE, x = c(1:nrow(phm_s)), y = c(1:ncol(phm_s)), xlab = "", ylab = "") + axis(side = 1, at = c(1:nrow(phm_s)),labels = rownames(hm_s), las = 2) + axis(side = 4, at = c(1:ncol(phm_s)),labels = colnames(hm_s), las = 2)
  dev.off()
}

rm(hm_s)
rm(hm_s_r)
rm(LD_s)
rm(p_foc_s)
rm(p_foc_s_r)
rm(p_alt_s)
rm(p_alt_s_r)
rm(phm_s)
rm(pvalues_foc_s)
rm(pvalues_self_s)
rm(pvalues_alt_s)
rm(snps_s)
rm(i)
rm(j)
rm(m)
rm(s)

for(s in 1:nrow(sig_SNPs_range)){
  snps_s = filter(snps_all, chr==sig_SNPs_range$chr[s]&pos>=sig_SNPs_range$pos_min[s]&pos<=sig_SNPs_range$pos_max[s])
  if(sig_SNPs_range$type[s]=="nei"){
    pvalues_nei_s = gwas_results_for_LD_p_hm[,c(1,2,(s*2)+2)] #these are actually -log p values
    pvalues_sxn_s = gwas_results_for_LD_p_hm[,c(1,2,(s*2)+3)] #these are actually -log p values
  }else{
    pvalues_nei_s = gwas_results_for_LD_p_hm[,c(1,2,(s*2)+3)] #these are actually -log p values
    pvalues_sxn_s = gwas_results_for_LD_p_hm[,c(1,2,(s*2)+2)] #these are actually -log p values
  }
  pvalues_nei_s = filter(pvalues_nei_s, chr==sig_SNPs_range$chr[s]&pos>=sig_SNPs_range$pos_min[s]&pos<=sig_SNPs_range$pos_max[s])
  pvalues_nei_s = pvalues_nei_s[order(pvalues_nei_s$pos),]
  p_nei_s = matrix(0.0, ncol = nrow(pvalues_nei_s), nrow = nrow(pvalues_nei_s))
  for(i in 1:nrow(pvalues_nei_s)){
    for(j in 1:nrow(pvalues_nei_s)){
      p_nei_s[i,j] = sqrt(pvalues_nei_s[i,3]+pvalues_nei_s[j,3]) #using sqrt for better visuality
    }
  }
  pvalues_sxn_s = filter(pvalues_sxn_s, chr==sig_SNPs_range$chr[s]&pos>=sig_SNPs_range$pos_min[s]&pos<=sig_SNPs_range$pos_max[s])
  pvalues_sxn_s = pvalues_sxn_s[order(pvalues_sxn_s$pos),]
  p_sxn_s = matrix(0.0, ncol = nrow(pvalues_sxn_s), nrow = nrow(pvalues_sxn_s))
  for(i in 1:nrow(pvalues_sxn_s)){
    for(j in 1:nrow(pvalues_sxn_s)){
      p_sxn_s[i,j] = sqrt(pvalues_sxn_s[i,3]+pvalues_sxn_s[j,3]) #using sqrt for better visuality
    }
  }
  pvalues_self_s = gwas_results_for_LD_p_hm[,c(1,2,3)] #these are actually -log p values
  pvalues_self_s = filter(pvalues_self_s, chr==sig_SNPs_range$chr[s]&pos>=sig_SNPs_range$pos_min[s]&pos<=sig_SNPs_range$pos_max[s])
  pvalues_self_s = pvalues_self_s[order(pvalues_self_s$pos),]
  p_self_s = matrix(0.0, ncol = nrow(pvalues_self_s), nrow = nrow(pvalues_self_s))
  for(i in 1:nrow(pvalues_self_s)){
    for(j in 1:nrow(pvalues_self_s)){
      p_self_s[i,j] = sqrt(pvalues_self_s[i,3]+pvalues_self_s[j,3]) #using sqrt for better visuality
    }
  }
  m = max(c(max(p_self_s),max(p_nei_s),max(p_sxn_s)))
  p_self_s = p_self_s/m
  p_nei_s = p_nei_s/m
  p_sxn_s = p_sxn_s/m
  p_self_s_r = round((p_self_s*32)-0.5)
  p_nei_s_r = round((p_nei_s*32)-0.5)
  p_sxn_s_r = round((p_sxn_s*32)-0.5)
  phm_s = matrix(0.0, ncol = ncol(p_self_s_r), nrow = nrow(p_self_s_r))
  for(i in 1:nrow(phm_s)){
    for(j in 1:ncol(phm_s)){
      if(p_self_s_r[i,j]==32){p_self_s_r[i,j]=31}
      if(p_nei_s_r[i,j]==32){p_nei_s_r[i,j]=31}
      if(p_sxn_s_r[i,j]==32){p_sxn_s_r[i,j]=31}
      phm_s[i,j] = (1024*p_self_s_r[i,j]) + (32*p_nei_s_r[i,j]) + p_sxn_s_r[i,j]
    }
  }
  colnames(phm_s)<-paste0(pvalues_nei_s[,1],"_",pvalues_nei_s[,2])
  rownames(phm_s)<-paste0(pvalues_nei_s[,1],"_",pvalues_nei_s[,2])
  png(filename = paste0("results/",species,"/LD_asym/region/",daynumber,"_",poolName,"/heatmaps/",poolName,"_self-p,_nei-p_and_sxn-p_heatmap_around_",sig_SNPs_range$type[s],"_",sig_SNPs_range$chr[s],"_",sig_SNPs_range$pos[s],".png"),width=1600,height=800,units="px")
  par(mar = c(9.1,1.1,1.1,9.1))
  image(phm_s, col = heatmap_cs_def, zlim = c(0,32767), axes = FALSE, x = c(1:nrow(phm_s)), y = c(1:ncol(phm_s)), xlab = "", ylab = "") + axis(side = 1, at = c(1:nrow(phm_s)),labels = rownames(phm_s), las = 2) + axis(side = 4, at = c(1:ncol(phm_s)),labels = colnames(phm_s), las = 2)
  dev.off()
}

rm(haplo)
rm(gwas_results_for_LD_p_hm)
rm(p_self_s)
rm(p_self_s_r)
rm(p_nei_s)
rm(p_nei_s_r)
rm(p_sxn_s)
rm(p_sxn_s_r)
rm(phm_s)
rm(pvalues_nei_s)
rm(pvalues_self_s)
rm(pvalues_sxn_s)
rm(snps_s)
rm(i)
rm(j)
rm(m)
rm(s)

