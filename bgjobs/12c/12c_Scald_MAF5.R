library(rNeighborGWAS)
library(gaston)
library(dplyr)
library(readr)

#additional plot for PVE multiple MAF extended
#run after "12b_PVE_analysis_extended" and "11a_automated_neiGWAS_extended"

#adds a rgb-type plot with -log(p)-values for comparison

################################# settings
#you need to adjust this before running

#day number for easier identification of outputs
daynumber = "20240513"
#use same daynumber as for "10d_PVE_analysis_MAF_extended"
SNP_analysis_daynumber = "20240513"
#use daynumber of SNP_analysis_extended

#pool settings:
species = "barley"
#barley or wheat
poolName = "Scald_pooled-MAF5"
poolName_MAFless = "Scald_pooled"
#barley pools:
#NFNB_pooled, PM_2015, Blr_pooled, Scald_pooled, SFNB_pooled
#wheat pools:
#YR_pooled, LR_2017, SR_2017, YLS_2017, STB_2019, YS_2019

################################# load and prepare data

PVE_results_multiMAF = read.csv(paste0("results/",species,"/PVE_analysis/multiMAF/extended/",daynumber,"_",poolName_MAFless,"/PVE_results_multiMAF_",poolName_MAFless,".csv"))
sig_SNPs = read.csv(paste0("results/",species,"/SNP_analysis_extended/sig_SNPs_vs_dist_",SNP_analysis_daynumber,"/significant_SNPs_all_dist.csv"))

sig_SNPs = filter(sig_SNPs, phenotype == poolName)
sig_SNPs = mutate(sig_SNPs, chr_pos = paste0(chr,"_",pos))
sig_SNPs = select(sig_SNPs, c(ncol(sig_SNPs),((4:((ncol(sig_SNPs)-1)/2))*2)+0))

distances = sort(unique(PVE_results_multiMAF$distance))
MAFs = sort(unique(PVE_results_multiMAF$MAF))

################################# plot

#distance vs. MAF rgb map

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

d_env = round(((PVE_results_multiMAF$PVEsnd-PVE_results_multiMAF$PVEsn)*32)-0.5)
d_nei = round(((PVE_results_multiMAF$PVEsn-PVE_results_multiMAF$PVEs)*32)-0.5)
d_self = round((PVE_results_multiMAF$PVEs*32)-0.5)
d_env = matrix(d_env, nrow = length(distances), ncol = length(MAFs))
d_nei = matrix(d_nei, nrow = length(distances), ncol = length(MAFs))
d_self = matrix(d_self, nrow = length(distances), ncol = length(MAFs))
plot_data_2 = matrix(0.0, ncol = nrow(sig_SNPs) + length(MAFs), nrow = length(distances))
for(i in 1:nrow(plot_data_2)){
  for(j in 1:length(MAFs)){
    if(d_env[i,j]>31){d_env[i,j]=31}
    else if(d_env[i,j]<0){d_env[i,j]=0}
    if(d_nei[i,j]>31){d_nei[i,j]=31}
    else if(d_nei[i,j]<0){d_nei[i,j]=0}
    if(d_self[i,j]>31){d_self[i,j]=31}
    else if(d_self[i,j]<0){d_self[i,j]=0}
    plot_data_2[i,j] = (1024*d_env[i,j]) + (32*d_nei[i,j]) + d_self[i,j]
  }
}
for(i in 1:nrow(sig_SNPs)){
  p_max_i = max(sig_SNPs[i,c(2:ncol(sig_SNPs))])
  for(j in 1:nrow(plot_data_2)){
    p_ij = round(((sig_SNPs[i,j+1]/p_max_i)*(sig_SNPs[i,j+1]/p_max_i)*16)-0.5)
    if(p_ij>31){p_ij=31}
    else if(p_ij<0){p_ij=0}
    plot_data_2[j,i+length(MAFs)] = 1057*p_ij + 1
  }
}

png(filename = paste0("results/",species,"/PVE_analysis/multiMAF/extended/",daynumber,"_",poolName_MAFless,"/",poolName,"_distance_vs_MAF_rgb_map_with_sigSNPs.png"),width=1600,height=800,units="px")
bgc <- par(bg = "black")
par(mar = c(8.0,4.1,4.1,9.1))#bottom,left,top,right
par(col.main = "white",col.axis = "white", col.lab = "white", col = "white")
image(plot_data_2, col = heatmap_cs_def, zlim = c(0,32767), axes = FALSE, x = c(1:nrow(plot_data_2)), y = c(1:ncol(plot_data_2)), xlab = "distance of included neighbours", ylab = "sig.SNPs or MAF", main = "PVE partitioning, distance vs. MAF, including sig. SNPs") + axis(side = 1, at = c(1:nrow(plot_data_2)),labels = as.factor(distances), las = 2) + axis(side = 4, at = c(1:ncol(plot_data_2)),labels = c(as.factor(MAFs),as.factor(sig_SNPs$chr_pos)), las = 2) + mtext("black to colour scale, env in red, nei in green, self in blue, with -log(p) of sig. SNP on top in white", col = "white")
dev.off()

#cleaning
rm(bgc)
rm(d_env)
rm(d_nei)
rm(d_self)
rm(plot_data_2)
rm(distances)
rm(heatmap_cs_def)
rm(i)
rm(j)
rm(MAFs)
rm(getStackedHeatmapColourmapRGB)
rm(p_ij)
rm(p_max_i)
rm(daynumber)
rm(SNP_analysis_daynumber)
rm(sig_SNPs)
rm(PVE_results_multiMAF)

