library(dplyr)
library(readr)

#this script is to automatically define candidate gene search ranges for loci

#for nei and sxn sig. SNPs

################################# settings
#you need to adjust this before running

#day number for easier identification of outputs
daynumber = "20240607"
#do not use slashes or points, since this will be used for the output folder name

#pool settings:
species = "barley"
#barley or wheat
poolName = "Blr_pooled"
#barley pools:
#NFNB_pooled, PM_2015, Blr_pooled, Scald_pooled, SFNB_pooled
#wheat pools:
#YR_pooled, LR_2017, SR_2017, YLS_2017, STB_2019, YS_2019

linkage_threshold = 0.5

buffer = 10000
#a minimum extent of a subrange from a linked SNP,
#since the effect-causing variant might be upstream or downstream of the affected gene

################################# loading and preparing data

#importing data
sigSNP_LD = read.csv(paste0("results/",species,"/sigSNP_asym/LD/data_collection/",poolName,"_SNPwise_LD_of_sig_SNPs.csv"))
gmap = read.csv(paste0("results/",species,"/filtered_inputs/data_collection/gmap_",poolName,".csv"))

sigSNP_LD = filter(sigSNP_LD, same_chr)
colnames(gmap)<- c("chr","pos")

################################# define loci

sig_SNPs = unique(sigSNP_LD$SNP_foc)

connected_sigSNPs = data.frame("SNP_foc"=numeric(),"SNP_comp"=numeric(),"valid"=logical())
for(i in 1:length(sig_SNPs)){
  sigSNP_LD_i = filter(sigSNP_LD, SNP_foc == sig_SNPs[i])
  sigSNP_LD_i = filter(sigSNP_LD_i, SNP_comp_sig)
  sigSNP_LD_i = filter(sigSNP_LD_i, r2 > linkage_threshold)
  sigSNP_LD_i = mutate(sigSNP_LD_i, "valid" = TRUE)
  sigSNP_LD_i = sigSNP_LD_i[,c(1,2,7)]
  sigSNP_LD_ij = sigSNP_LD_i
  for(j in 1:nrow(sigSNP_LD_i)){
    connected_sigSNPs_i = which(connected_sigSNPs$SNP_foc == sigSNP_LD_i[j,2])
    if(length(connected_sigSNPs_i)!=0){
      for(k in 1:length(connected_sigSNPs_i)){
        connected_sigSNPs[connected_sigSNPs_i[k],3] = FALSE
        if(!(connected_sigSNPs[connected_sigSNPs_i[k],2]%in%sigSNP_LD_i[,2])){
          sigSNP_LD_ik = data.frame("SNP_foc"=sigSNP_LD_i[1,1],"SNP_comp"=connected_sigSNPs[connected_sigSNPs_i[k],2],"valid"=TRUE)
          sigSNP_LD_ij = rows_append(sigSNP_LD_ij, sigSNP_LD_ik)
        }
      }
    }
  }
  connected_sigSNPs = rows_append(connected_sigSNPs, sigSNP_LD_ij)
  connected_sigSNPs = filter(connected_sigSNPs, valid)
}
lead_sigSNPs = unique(connected_sigSNPs$SNP_foc)
connected_sigSNPs = mutate(connected_sigSNPs, locus_ID = sapply(SNP_foc, function(x) which(lead_sigSNPs == x)))
sigSNPs_loci = select(connected_sigSNPs, c(4,2))
colnames(sigSNPs_loci)<- c("locus","sig_SNP")

#clean
rm(i)
rm(j)
rm(k)
rm(sigSNP_LD_i)
rm(sigSNP_LD_ij)
rm(sigSNP_LD_ik)
rm(connected_sigSNPs_i)
rm(connected_sigSNPs)
rm(lead_sigSNPs)

################################# generate ranges

loci = unique(sigSNPs_loci$locus)
subranges = data.frame("locus"=numeric(),"chr"=character(),"pos_l"=numeric(),"pos_r"=numeric())
for(i in 1:length(loci)){
  sigSNPs_loci_i = filter(sigSNPs_loci, locus == loci[i])
  linked_SNPs = {}
  for(j in 1:nrow(sigSNPs_loci_i)){
    sigSNP_LD_j = filter(sigSNP_LD, SNP_foc == sigSNPs_loci_i[j,2])
    linked_SNPs_j = which(sigSNP_LD_j$r2 > linkage_threshold)
    linked_SNPs = c(linked_SNPs, sigSNP_LD_j[linked_SNPs_j,2])
  }
  linked_SNPs = unique(linked_SNPs)
  linked_SNPs = linked_SNPs[order(linked_SNPs)]
  #
  subranges_i = data.frame("left"=numeric(),"right"=numeric())
  for(j in 1:length(linked_SNPs)){
    subrange_j = data.frame("left"=(linked_SNPs[j]-1),"right"=(linked_SNPs[j]+1))
    subranges_i = rows_append(subranges_i, subrange_j)
  }
  subranges_i = mutate(subranges_i, "chr" = gmap[left,1])
  subranges_i = mutate(subranges_i, "pos_l" = gmap[left,2]-buffer)
  subranges_i = mutate(subranges_i, "pos_r" = gmap[right,2]+buffer)
  subranges_i = mutate(subranges_i, "valid" = TRUE)
  if(nrow(subranges_i)>1){
    for(j in 1:(nrow(subranges_i)-1)){
      if(subranges_i[j,6]){
        k = 1
        while((subranges_i[j,5]>subranges_i[j+k,4])&(k>0)){
          subranges_i[j,5] = subranges_i[j+k,5]
          subranges_i[j+k,6] = FALSE
          k = k+1
          if(j+k>nrow(subranges_i)){k=0}
        }
      }
    }
  }
  subranges_i = filter(subranges_i, valid)
  subranges_i = mutate(subranges_i, "locus" = i)
  subranges_i = select(subranges_i, c(7,3,4,5))
  subranges = rows_append(subranges, subranges_i)
}

sigSNPs_loci = mutate(sigSNPs_loci, "chr" = gmap[sig_SNP,1])
sigSNPs_loci = mutate(sigSNPs_loci, "pos" = gmap[sig_SNP,2])
sigSNPs_loci = mutate(sigSNPs_loci, "phenotype" = poolName)
sigSNPs_loci = select(sigSNPs_loci, c(5,3,4,1))

locus_info = data.frame("locus"=loci)
locus_info = mutate(locus_info, "sig_SNPs" = as.data.frame(table(sigSNPs_loci$locus))[,2])
locus_info = mutate(locus_info, "subloci" = as.data.frame(table(subranges$locus))[,2])
locus_info = mutate(locus_info, "range" = as.data.frame((rowsum(subranges$pos_r, subranges$locus) - rowsum(subranges$pos_l, subranges$locus)))[,1])
locus_info = mutate(locus_info, "LD_threshold" = linkage_threshold)
locus_info = mutate(locus_info, "buffer_region" = buffer)

sigSNP_info = read.csv(paste0("results/",species,"/sigSNP_asym/LD/data_collection/",poolName,"_sig_SNPs_with_inter_LD.csv"))
sigSNP_info = sigSNP_info[,1:10]
sigSNP_info = mutate(sigSNP_info, "chr_pos" = paste0(chr,"_",pos))
sigSNPs_loci = mutate(sigSNPs_loci, "chr_pos" = paste0(chr,"_",pos))
sigSNPs_loci = sigSNPs_loci[,c(4,5)]
sigSNP_info = merge(sigSNP_info, sigSNPs_loci, by = "chr_pos")
sigSNP_info = select(sigSNP_info, c(2,3,4,5,12,6,7,8,9,10,11))
sigSNP_info = arrange(sigSNP_info, -p_log10)

#clean
rm(sigSNP_LD)
rm(sigSNP_LD_j)
rm(sigSNPs_loci)
rm(sigSNPs_loci_i)
rm(subrange_j)
rm(subranges_i)
rm(i)
rm(j)
rm(k)
rm(linked_SNPs)
rm(linked_SNPs_j)
rm(loci)
rm(sig_SNPs)
rm(gmap)

################################# save and plot

#create output folder if necessary
if(!file.exists(paste0("results/",species,"/auto_searchrange"))){
  dir.create(path = paste0("results/",species,"/auto_searchrange"))
}
if(!file.exists(paste0("results/",species,"/auto_searchrange/",daynumber,"_",poolName))){
  dir.create(path = paste0("results/",species,"/auto_searchrange/",daynumber,"_",poolName))
}

write_csv(sigSNP_info, file = paste0("results/",species,"/auto_searchrange/",daynumber,"_",poolName,"/",poolName,"_sig_SNPs_with_loci.csv"))
write_csv(subranges, file = paste0("results/",species,"/auto_searchrange/",daynumber,"_",poolName,"/",poolName,"_subranges.csv"))
write_csv(locus_info, file = paste0("results/",species,"/auto_searchrange/",daynumber,"_",poolName,"/",poolName,"_locus_info.csv"))

gwas_results = read.csv(paste0("./results/",species,"/neiGWAS_asym/data_collection/neiGWAS_results_asym_",poolName,".csv"))
nSNPs = nrow(gwas_results)
gwas_results = arrange(gwas_results, pos)
gwas_results = arrange(gwas_results, chr)
gwas_results = mutate(gwas_results, "phenotype" = poolName)
gwas_results = select(gwas_results, c(1,2,((2:(ncol(gwas_results)/3))*3)-2))
nei_dist = sapply(strsplit(colnames(gwas_results[,(1:((ncol(gwas_results)-2)/3))*3]),split="_",fixed=TRUE), function(x) (x[4]))

library(ggplot2)

#make custom plot with -log(p) self, nei, sxn and search range in background
for(i in 1:nrow(locus_info)){
  subranges_i = filter(subranges, locus == locus_info[i,1])
  pos_min = min(subranges_i[,3])
  pos_max = max(subranges_i[,4])
  plot_range = (pos_max - pos_min)*1.5
  pos_min = pos_min - (plot_range/6)
  pos_max = pos_max + (plot_range/6)
  gwas_results_i = filter(gwas_results, chr == subranges_i[1,2])
  gwas_results_i = filter(gwas_results_i, pos>pos_min&pos<pos_max)
  sigSNPs_i = filter(sigSNP_info, chr == subranges_i[1,2])
  sigSNPs_i = filter(sigSNPs_i, pos>pos_min&pos<pos_max)
  nei_dist_i = which(as.double(nei_dist)==as.double(sigSNPs_i[1,6]))
  gwas_results_i = select(gwas_results_i, c(1,2,0:2+(nei_dist_i*3)))
  #plot
  png(filename = paste0("results/",species,"/auto_searchrange/",daynumber,"_",poolName,"/",poolName,"_",sigSNPs_i[1,4],"_",sigSNPs_i[1,5],"_",sigSNPs_i[1,2],"_",sigSNPs_i[1,3],"_search_range.png"),width=1800,height=1000,units="px")
  bgc <- par(bg = "white")
  par(col.main = "black",col.axis = "black", col.lab = "black")
  plot(c(pos_min, pos_max), c(0.0, 1.1*max(gwas_results_i[,c(3,4,5)])), type = "n", xlab = paste0("position on chromosome ",sigSNPs_i[1,2]), ylab = "-log(p) (self in red, nei in green, sxn in blue)", main = paste0("-log(p)-values in automatically genereated search range"), yaxt = "n", xaxt = "n")
  axis(side = 1, col = "black")
  axis(side = 2, col = "black")
  for(k in 1:nrow(subranges_i)){
    rect(subranges_i$pos_l[k], 0.0, subranges_i$pos_r[k], 1.1*max(gwas_results_i[,c(3,4,5)]), col = rgb(0.9, 0.9, 0.9), border = rgb(0.5, 0.5, 0.5))
  }
  abline(h = -log10(0.05/nSNPs), col = rgb(0.0,0.0,0.0))
  lines(x = gwas_results_i[,2], y = gwas_results_i[,3], col = rgb(1.0,0.0,0.0))
  lines(x = gwas_results_i[,2], y = gwas_results_i[,4], col = rgb(0.0,0.8,0.0))
  lines(x = gwas_results_i[,2], y = gwas_results_i[,5], col = rgb(0.1,0.1,1.0))
  mtext(paste0("phenotype: ",sigSNPs_i[1,1],", range: ",sigSNPs_i[1,2],":",pos_min,"-",pos_max,", type: ",sigSNPs_i[1,4],", SNPs: ",nrow(sigSNPs_i)), col = "black")
  dev.off()
}

rm(bgc)
rm(gwas_results)
rm(gwas_results_i)
rm(sigSNPs_i)
rm(subranges_i)
rm(buffer)
rm(i)
rm(k)
rm(linkage_threshold)
rm(nei_dist)
rm(nei_dist_i)
rm(nSNPs)
rm(plot_range)
rm(pos_max)
rm(pos_min)

