library(dplyr)
library(readr)

#this script is to plot -log(p) and LD around sig. SNPs,
#using haploid genotypes for better LD,
#to decide which SNPs belong to a common loci and
#get an idea in which ranges to look for candidate genes

#for nei and sxn sig. SNPs

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
#before running check manually in folders if all required files are present

#importing data
haplo = read.csv(paste0("results/",species,"/filtered_inputs/data_collection/haplo_",poolName,".csv"))
sig_SNPs_nei = read.csv(paste0("results/",species,"/sigSNP_asym/data_collection/",poolName,"_nei_sig_SNPs_with_bases.csv"))
sig_SNPs_sxn = read.csv(paste0("results/",species,"/sigSNP_asym/data_collection/",poolName,"_sxn_sig_SNPs_with_bases.csv"))

colnames(sig_SNPs_nei)<- c("phenotype","chr","pos","top_dist","p_log10","effect","base_ref","base_alt")
colnames(sig_SNPs_sxn)<- c("phenotype","chr","pos","top_dist","p_log10","effect","base_ref","base_alt")
sig_SNPs_nei = mutate(sig_SNPs_nei, "type" = "nei")
sig_SNPs_sxn = mutate(sig_SNPs_sxn, "type" = "sxn")
sig_SNPs_nei$top_dist = as.double(sig_SNPs_nei$top_dist)
sig_SNPs_sxn$top_dist = as.double(sig_SNPs_sxn$top_dist)
if(nrow(sig_SNPs_nei)>0){
  if(nrow(sig_SNPs_sxn)>0){
    sig_SNPs = rows_append(sig_SNPs_nei, sig_SNPs_sxn)
  }else{
    sig_SNPs = sig_SNPs_nei
  }
}else{
  sig_SNPs = sig_SNPs_sxn
}
sig_SNPs = select(sig_SNPs, c(1,2,3,9,4,5,6,7,8))
sig_SNPs = arrange(sig_SNPs, -p_log10)

rm(sig_SNPs_nei)
rm(sig_SNPs_sxn)

################################# calculating LD, only for significant SNPs
#this can take a while...

LD = data.frame(NULL)
nhaplo = ncol(haplo)-2
for(sigsnp_i in 1:nrow(sig_SNPs)){
  LD_i = data.frame(NULL)
  snp_i=which(haplo$chr==sig_SNPs$chr[sigsnp_i]&haplo$pos==sig_SNPs$pos[sigsnp_i])
  for(snp_j in 1:nrow(haplo)){
    counts = c(0,0,0)#I&J,I,J
    for(k in 3:ncol(haplo)){
      counts[1] = counts[1] + (haplo[snp_i,k]*haplo[snp_j,k])
      counts[2] = counts[2] + (haplo[snp_i,k])
      counts[3] = counts[3] + (haplo[snp_j,k])
    }
    r2 = (((nhaplo*counts[1])-(counts[2]*counts[3]))*((nhaplo*counts[1])-(counts[2]*counts[3])))/(((nhaplo*counts[2])-(counts[2]*counts[2]))*((nhaplo*counts[3])-(counts[3]*counts[3])))
    scSNP = FALSE
    for(snp_c in row(sig_SNPs)){
      if(ifelse(sig_SNPs$chr[snp_c]==haplo$chr[snp_j]&sig_SNPs$pos[snp_c]==haplo$pos[snp_j],TRUE,FALSE)){
        scSNP = TRUE
      }
    }
    samechr = haplo[snp_i,1] == haplo[snp_j,1]
    distchr = 0
    if(samechr){
      distchr = haplo[snp_j,2] - haplo[snp_i,2]
    }
    LD_ij = data.frame("SNP_foc"=rownames(haplo)[snp_i],"SNP_comp"=rownames(haplo)[snp_j],"r2"=r2,"SNP_comp_sig"=scSNP,"same_chr"=samechr,"dist"=distchr)
    if(snp_j==1){
      LD_i = LD_ij
    }
    else{
      LD_i = rows_append(LD_i, LD_ij)
    }
  }
  if(sigsnp_i==1){
    LD = LD_i
    snp_ID = snp_i
  }
  else{
    LD = rows_append(LD, LD_i)
    snp_ID = c(snp_ID, snp_i)
  }
}
sig_SNPs = mutate(sig_SNPs, SNP_ID = snp_ID)

rm(nhaplo)
rm(LD_i)
rm(counts)
rm(r2)
rm(LD_ij)
rm(sigsnp_i)
rm(snp_i)
rm(snp_j)
rm(k)
rm(scSNP)
rm(snp_c)
rm(samechr)
rm(distchr)
rm(snp_ID)

################################# saving data

#create output folder if necessary
if(!file.exists(paste0("results/",species,"/sigSNP_asym/LD"))){
  dir.create(path = paste0("results/",species,"/sigSNP_asym/LD"))
}
if(!file.exists(paste0("results/",species,"/sigSNP_asym/LD/data_collection"))){
  dir.create(path = paste0("results/",species,"/sigSNP_asym/LD/data_collection"))
}
if(!file.exists(paste0("results/",species,"/sigSNP_asym/LD/",daynumber,"_",poolName))){
  dir.create(path = paste0("results/",species,"/sigSNP_asym/LD/",daynumber,"_",poolName))
}
if(!file.exists(paste0("results/",species,"/sigSNP_asym/LD/",daynumber,"_",poolName,"/p"))){
  dir.create(path = paste0("results/",species,"/sigSNP_asym/LD/",daynumber,"_",poolName,"/p"))
}

write_csv(LD, file = paste0("results/",species,"/sigSNP_asym/LD/",daynumber,"_",poolName,"/",poolName,"_SNPwise_LD_of_sig_SNPs.csv"))
#manually copy to results/",species,"/sigSNP_asym/LD/data_collection for further use

for(snp_i in 1:nrow(sig_SNPs)){
  LD_i = filter(LD, SNP_comp==sig_SNPs$SNP_ID[snp_i])
  vLD_i = {}
  for(snp_j in 1:nrow(sig_SNPs)){
    LD_ij = filter(LD_i, SNP_foc==sig_SNPs$SNP_ID[snp_j])
    vLD_i = c(vLD_i, LD_ij$r2)
  }
  cns = colnames(sig_SNPs)
  sig_SNPs = mutate(sig_SNPs, newcol = vLD_i)
  colnames(sig_SNPs)<-c(cns, paste0("LD_SNP_",sig_SNPs$SNP_ID[snp_i]))
}
rm(cns)
rm(LD_i)
rm(vLD_i)
rm(LD_ij)
rm(snp_i)
rm(snp_j)

write_csv(sig_SNPs, file = paste0("results/",species,"/sigSNP_asym/LD/",daynumber,"_",poolName,"/",poolName,"_sig_SNPs_with_inter_LD.csv"))
#manually copy to results/",species,"/sigSNP_asym/LD/data_collection for further use

################################# plot -log(p) and LD

gwas_results = read.csv(paste0("./results/",species,"/neiGWAS_asym/data_collection/neiGWAS_results_asym_",poolName,".csv"))
nSNPs = nrow(gwas_results)
gwas_results = arrange(gwas_results, pos)
gwas_results = arrange(gwas_results, chr)
gwas_results = mutate(gwas_results, "phenotype" = poolName)

library(ggplot2)

dist_threshold_big = 100000000 #limit plot size to interesting area
dist_threshold_medium = 5000000
dist_threshold_small = 250000

#prepare -log(p)-plots
p_plot = data.frame("p_log10"=numeric(),"SNP_foc"=character(),"SNP_comp"=character(),"dist"=numeric(),"nei_dist"=numeric())
for(snp_i in 1:nrow(sig_SNPs)){
  gwas_results_i = select(gwas_results, c("phenotype","chr","pos",paste0("p_",sig_SNPs[snp_i,4],"_log10_",sig_SNPs[snp_i,5])))
  colnames(gwas_results_i)<- c("phenotype","chr","pos","p_log10")
  gwas_results_i = filter(gwas_results_i, chr == sig_SNPs[snp_i,2])
  gwas_results_i = mutate(gwas_results_i, "SNP_foc" = paste0(sig_SNPs[snp_i,2],"_",sig_SNPs[snp_i,3]))
  gwas_results_i = mutate(gwas_results_i, "SNP_comp" = paste0(chr,"_",pos))
  gwas_results_i = mutate(gwas_results_i, "dist" = sapply(pos, function(x) as.numeric((x-sig_SNPs[snp_i,3]))))
  gwas_results_i = filter(gwas_results_i, abs(dist)<dist_threshold_big)
  gwas_results_i = mutate(gwas_results_i, "nei_dist" = sig_SNPs[snp_i,5])
  gwas_results_i = select(gwas_results_i, c(4,5,6,7,8))
  p_plot = rows_append(p_plot, gwas_results_i)
}
focal_SNPs = unique(p_plot$SNP_foc)
#plot -log(p) arround sig SNPs, big
sigSNP_p_plot = ggplot(p_plot) + geom_line(aes(x=dist,y=p_log10,group=SNP_foc,colour=SNP_foc)) + ggtitle("-log10(p) of SNPs around sig. SNPs") + geom_point(aes(x=dist,y=p_log10,group=SNP_foc,colour=SNP_foc)) + xlab("distance on chromosome (bp)") + ylab("-log10(p)") + geom_hline(yintercept = -log10(0.05/nSNPs))
ggsave(filename = paste0("results/",species,"/sigSNP_asym/LD/",daynumber,"_",poolName,"/p/",poolName,"_-log10_p_around_sig_SNPs_big.png"),plot=sigSNP_p_plot,width=4000,height=3000,units="px",device="png")
#plot -log(p) per SNP
for(snp_i in 1:nrow(sig_SNPs)){
  p_plot_i = filter(p_plot, SNP_foc == focal_SNPs[snp_i])
  snp_name = focal_SNPs[snp_i]
  sigSNP_p_plot_i = ggplot(p_plot_i) + geom_smooth(aes(x=abs(dist),y=p_log10,group=SNP_foc,colour=SNP_foc),method="gam",se=FALSE) + ggtitle(paste0("-log10(p) of SNPs around ",snp_name)) + geom_point(aes(x=abs(dist),y=p_log10,group=SNP_foc,colour=SNP_foc)) + xlab("absolute distance on chromosome (bp)") + ylab("-log10(p)") + geom_hline(yintercept = -log10(0.05/nSNPs))
  ggsave(filename = paste0("results/",species,"/sigSNP_asym/LD/",daynumber,"_",poolName,"/p/",poolName,"_-log10_p_around_",snp_name,"_onesided_big.png"),plot=sigSNP_p_plot_i,width=4000,height=3000,units="px",device="png")
  rm(sigSNP_p_plot_i)
  sigSNP_p_plot_i = ggplot(p_plot_i) + geom_line(aes(x=dist,y=p_log10,group=SNP_foc,colour=SNP_foc)) + ggtitle(paste0("-log10(p) of SNPs around ",snp_name)) + geom_point(aes(x=dist,y=p_log10,group=SNP_foc,colour=SNP_foc)) + xlab("absolute distance on chromosome (bp)") + ylab("-log10(p)") + geom_hline(yintercept = -log10(0.05/nSNPs))
  ggsave(filename = paste0("results/",species,"/sigSNP_asym/LD/",daynumber,"_",poolName,"/p/",poolName,"_-log10_p_around_",snp_name,"_big.png"),plot=sigSNP_p_plot_i,width=4000,height=3000,units="px",device="png")
  rm(sigSNP_p_plot_i)
}
rm(sigSNP_p_plot)
p_plot = filter(p_plot, abs(dist)<dist_threshold_medium)
#plot -log(p) arround sig SNPs, medium
sigSNP_p_plot = ggplot(p_plot) + geom_line(aes(x=dist,y=p_log10,group=SNP_foc,colour=SNP_foc)) + ggtitle("-log10(p) of SNPs around sig. SNPs") + geom_point(aes(x=dist,y=p_log10,group=SNP_foc,colour=SNP_foc)) + xlab("distance on chromosome (bp)") + ylab("-log10(p)") + geom_hline(yintercept = -log10(0.05/nSNPs))
ggsave(filename = paste0("results/",species,"/sigSNP_asym/LD/",daynumber,"_",poolName,"/p/",poolName,"_-log10_p_around_sig_SNPs_medium.png"),plot=sigSNP_p_plot,width=4000,height=3000,units="px",device="png")
#plot -log(p) per SNP
for(snp_i in 1:nrow(sig_SNPs)){
  p_plot_i = filter(p_plot, SNP_foc == focal_SNPs[snp_i])
  snp_name = focal_SNPs[snp_i]
  sigSNP_p_plot_i = ggplot(p_plot_i) + geom_smooth(aes(x=abs(dist),y=p_log10,group=SNP_foc,colour=SNP_foc),method="gam",se=FALSE) + ggtitle(paste0("-log10(p) of SNPs around ",snp_name)) + geom_point(aes(x=abs(dist),y=p_log10,group=SNP_foc,colour=SNP_foc)) + xlab("absolute distance on chromosome (bp)") + ylab("-log10(p)") + geom_hline(yintercept = -log10(0.05/nSNPs))
  ggsave(filename = paste0("results/",species,"/sigSNP_asym/LD/",daynumber,"_",poolName,"/p/",poolName,"_-log10_p_around_",snp_name,"_onesided_medium.png"),plot=sigSNP_p_plot_i,width=4000,height=3000,units="px",device="png")
  rm(sigSNP_p_plot_i)
  sigSNP_p_plot_i = ggplot(p_plot_i) + geom_line(aes(x=dist,y=p_log10,group=SNP_foc,colour=SNP_foc)) + ggtitle(paste0("-log10(p) of SNPs around ",snp_name)) + geom_point(aes(x=dist,y=p_log10,group=SNP_foc,colour=SNP_foc)) + xlab("absolute distance on chromosome (bp)") + ylab("-log10(p)") + geom_hline(yintercept = -log10(0.05/nSNPs))
  ggsave(filename = paste0("results/",species,"/sigSNP_asym/LD/",daynumber,"_",poolName,"/p/",poolName,"_-log10_p_around_",snp_name,"_medium.png"),plot=sigSNP_p_plot_i,width=4000,height=3000,units="px",device="png")
  rm(sigSNP_p_plot_i)
}
rm(sigSNP_p_plot)
p_plot = filter(p_plot, abs(dist)<dist_threshold_small)
#plot -log(p) arround sig SNPs, small
sigSNP_p_plot = ggplot(p_plot) + geom_line(aes(x=dist,y=p_log10,group=SNP_foc,colour=SNP_foc)) + ggtitle("-log10(p) of SNPs around sig. SNPs") + geom_point(aes(x=dist,y=p_log10,group=SNP_foc,colour=SNP_foc)) + xlab("distance on chromosome (bp)") + ylab("-log10(p)") + geom_hline(yintercept = -log10(0.05/nSNPs))
ggsave(filename = paste0("results/",species,"/sigSNP_asym/LD/",daynumber,"_",poolName,"/p/",poolName,"_-log10_p_around_sig_SNPs_small.png"),plot=sigSNP_p_plot,width=4000,height=3000,units="px",device="png")
#plot -log(p) per SNP
for(snp_i in 1:nrow(sig_SNPs)){
  p_plot_i = filter(p_plot, SNP_foc == focal_SNPs[snp_i])
  snp_name = focal_SNPs[snp_i]
  sigSNP_p_plot_i = ggplot(p_plot_i) + geom_smooth(aes(x=abs(dist),y=p_log10,group=SNP_foc,colour=SNP_foc),method="gam",se=FALSE) + ggtitle(paste0("-log10(p) of SNPs around ",snp_name)) + geom_point(aes(x=abs(dist),y=p_log10,group=SNP_foc,colour=SNP_foc)) + xlab("absolute distance on chromosome (bp)") + ylab("-log10(p)") + geom_hline(yintercept = -log10(0.05/nSNPs))
  ggsave(filename = paste0("results/",species,"/sigSNP_asym/LD/",daynumber,"_",poolName,"/p/",poolName,"_-log10_p_around_",snp_name,"_onesided_small.png"),plot=sigSNP_p_plot_i,width=4000,height=3000,units="px",device="png")
  rm(sigSNP_p_plot_i)
  sigSNP_p_plot_i = ggplot(p_plot_i) + geom_line(aes(x=dist,y=p_log10,group=SNP_foc,colour=SNP_foc)) + ggtitle(paste0("-log10(p) of SNPs around ",snp_name)) + geom_point(aes(x=dist,y=p_log10,group=SNP_foc,colour=SNP_foc)) + xlab("absolute distance on chromosome (bp)") + ylab("-log10(p)") + geom_hline(yintercept = -log10(0.05/nSNPs))
  ggsave(filename = paste0("results/",species,"/sigSNP_asym/LD/",daynumber,"_",poolName,"/p/",poolName,"_-log10_p_around_",snp_name,"_small.png"),plot=sigSNP_p_plot_i,width=4000,height=3000,units="px",device="png")
  rm(sigSNP_p_plot_i)
}
rm(sigSNP_p_plot)

#prepare LD plots
LD_plot = filter(LD, same_chr&abs(dist)<dist_threshold_big)
LD_plot = mutate(LD_plot, plot_gr = ifelse(SNP_comp_sig,paste0(sapply(SNP_comp, function(x) sig_SNPs[which(sig_SNPs$SNP_ID==x),2]),"_",sapply(SNP_comp, function(x) sig_SNPs[which(sig_SNPs$SNP_ID==x),3])),"others"))
LD_plot = mutate(LD_plot, line_gr = paste0(sapply(SNP_foc, function(x) sig_SNPs[which(sig_SNPs$SNP_ID==x),2]),"_",sapply(SNP_foc, function(x) sig_SNPs[which(sig_SNPs$SNP_ID==x),3])))
#plot LD arround sig SNPs, big
sigSNP_LD_plot = ggplot(LD_plot) + geom_line(aes(x=dist,y=r2,group=line_gr,colour=line_gr)) + ggtitle("LD of SNPs around sig. SNPs") + geom_point(aes(x=dist,y=r2,group=plot_gr,colour=plot_gr)) + xlab("distance on chromosome (bp)") + ylab("LD")
ggsave(filename = paste0("results/",species,"/sigSNP_asym/LD/",daynumber,"_",poolName,"/",poolName,"_LD_around_sig_SNPs_big.png"),plot=sigSNP_LD_plot,width=4000,height=3000,units="px",device="png")
#plot with logarithmic y axis
sigSNP_LD_log_plot = ggplot(LD_plot) + geom_line(aes(x=dist,y=log10(r2),group=line_gr,colour=line_gr)) + ggtitle("LD of SNPs around sig. SNPs") + geom_point(aes(x=dist,y=log10(r2),group=plot_gr,colour=plot_gr)) + xlab("distance on chromosome (bp)") + ylab("log10 LD")
ggsave(filename = paste0("results/",species,"/sigSNP_asym/LD/",daynumber,"_",poolName,"/",poolName,"_log_LD_around_sig_SNPs_big.png"),plot=sigSNP_LD_log_plot,width=4000,height=3000,units="px",device="png")
#plot one-sided and two-sided LD decay per SNP
for(snp_i in 1:nrow(sig_SNPs)){
  LD_plot_i = filter(LD_plot, SNP_foc == sig_SNPs$SNP_ID[snp_i])
  snp_name = paste0(sig_SNPs$chr[snp_i],"_",sig_SNPs$pos[snp_i])
  sigSNP_LD_plot_i = ggplot(LD_plot_i) + geom_smooth(aes(x=abs(dist),y=r2,group=line_gr,colour=line_gr),method="gam",se=FALSE) + ggtitle(paste0("LD of SNPs around ",snp_name)) + geom_point(aes(x=abs(dist),y=r2,group=plot_gr,colour=plot_gr)) + xlab("absolute distance on chromosome (bp)") + ylab("LD (r squared)")
  ggsave(filename = paste0("results/",species,"/sigSNP_asym/LD/",daynumber,"_",poolName,"/",poolName,"_LD_around_",snp_name,"_big.png"),plot=sigSNP_LD_plot_i,width=4000,height=3000,units="px",device="png")
  rm(sigSNP_LD_plot_i)
  sigSNP_LD_plot_i = ggplot(LD_plot_i) + geom_line(aes(x=dist,y=r2,group=line_gr,colour=line_gr)) + ggtitle(paste0("LD of SNPs around ",snp_name)) + geom_point(aes(x=dist,y=r2,group=plot_gr,colour=plot_gr)) + xlab("absolute distance on chromosome (bp)") + ylab("LD (r squared)")
  ggsave(filename = paste0("results/",species,"/sigSNP_asym/LD/",daynumber,"_",poolName,"/",poolName,"_LD_around_",snp_name,"_twosided_big.png"),plot=sigSNP_LD_plot_i,width=4000,height=3000,units="px",device="png")
  rm(sigSNP_LD_plot_i)
}
rm(sigSNP_LD_plot)
rm(sigSNP_LD_log_plot)
LD_plot = filter(LD_plot, abs(dist)<dist_threshold_medium)
#plot LD arround sig SNPs, medium
sigSNP_LD_plot = ggplot(LD_plot) + geom_line(aes(x=dist,y=r2,group=line_gr,colour=line_gr)) + ggtitle("LD of SNPs around sig. SNPs") + geom_point(aes(x=dist,y=r2,group=plot_gr,colour=plot_gr)) + xlab("distance on chromosome (bp)") + ylab("LD")
ggsave(filename = paste0("results/",species,"/sigSNP_asym/LD/",daynumber,"_",poolName,"/",poolName,"_LD_around_sig_SNPs_medium.png"),plot=sigSNP_LD_plot,width=4000,height=3000,units="px",device="png")
#plot with logarithmic y axis
sigSNP_LD_log_plot = ggplot(LD_plot) + geom_line(aes(x=dist,y=log10(r2),group=line_gr,colour=line_gr)) + ggtitle("LD of SNPs around sig. SNPs") + geom_point(aes(x=dist,y=log10(r2),group=plot_gr,colour=plot_gr)) + xlab("distance on chromosome (bp)") + ylab("log10 LD")
ggsave(filename = paste0("results/",species,"/sigSNP_asym/LD/",daynumber,"_",poolName,"/",poolName,"_log_LD_around_sig_SNPs_medium.png"),plot=sigSNP_LD_log_plot,width=4000,height=3000,units="px",device="png")
#plot one-sided and two-sided LD decay per SNP
for(snp_i in 1:nrow(sig_SNPs)){
  LD_plot_i = filter(LD_plot, SNP_foc == sig_SNPs$SNP_ID[snp_i])
  snp_name = paste0(sig_SNPs$chr[snp_i],"_",sig_SNPs$pos[snp_i])
  sigSNP_LD_plot_i = ggplot(LD_plot_i) + geom_smooth(aes(x=abs(dist),y=r2,group=line_gr,colour=line_gr),method="gam",se=FALSE) + ggtitle(paste0("LD of SNPs around ",snp_name)) + geom_point(aes(x=abs(dist),y=r2,group=plot_gr,colour=plot_gr)) + xlab("absolute distance on chromosome (bp)") + ylab("LD (r squared)")
  ggsave(filename = paste0("results/",species,"/sigSNP_asym/LD/",daynumber,"_",poolName,"/",poolName,"_LD_around_",snp_name,"_medium.png"),plot=sigSNP_LD_plot_i,width=4000,height=3000,units="px",device="png")
  rm(sigSNP_LD_plot_i)
  sigSNP_LD_plot_i = ggplot(LD_plot_i) + geom_line(aes(x=dist,y=r2,group=line_gr,colour=line_gr)) + ggtitle(paste0("LD of SNPs around ",snp_name)) + geom_point(aes(x=dist,y=r2,group=plot_gr,colour=plot_gr)) + xlab("absolute distance on chromosome (bp)") + ylab("LD (r squared)")
  ggsave(filename = paste0("results/",species,"/sigSNP_asym/LD/",daynumber,"_",poolName,"/",poolName,"_LD_around_",snp_name,"_twosided_medium.png"),plot=sigSNP_LD_plot_i,width=4000,height=3000,units="px",device="png")
  rm(sigSNP_LD_plot_i)
}
rm(sigSNP_LD_plot)
rm(sigSNP_LD_log_plot)
LD_plot = filter(LD_plot, abs(dist)<dist_threshold_small)
#plot LD arround sig SNPs, small
sigSNP_LD_plot = ggplot(LD_plot) + geom_line(aes(x=dist,y=r2,group=line_gr,colour=line_gr)) + ggtitle("LD of SNPs around sig. SNPs") + geom_point(aes(x=dist,y=r2,group=plot_gr,colour=plot_gr)) + xlab("distance on chromosome (bp)") + ylab("LD")
ggsave(filename = paste0("results/",species,"/sigSNP_asym/LD/",daynumber,"_",poolName,"/",poolName,"_LD_around_sig_SNPs_small.png"),plot=sigSNP_LD_plot,width=4000,height=3000,units="px",device="png")
#plot with logarithmic y axis
sigSNP_LD_log_plot = ggplot(LD_plot) + geom_line(aes(x=dist,y=log10(r2),group=line_gr,colour=line_gr)) + ggtitle("LD of SNPs around sig. SNPs") + geom_point(aes(x=dist,y=log10(r2),group=plot_gr,colour=plot_gr)) + xlab("distance on chromosome (bp)") + ylab("log10 LD")
ggsave(filename = paste0("results/",species,"/sigSNP_asym/LD/",daynumber,"_",poolName,"/",poolName,"_log_LD_around_sig_SNPs_small.png"),plot=sigSNP_LD_log_plot,width=4000,height=3000,units="px",device="png")
#plot one-sided and two-sided LD decay per SNP
for(snp_i in 1:nrow(sig_SNPs)){
  LD_plot_i = filter(LD_plot, SNP_foc == sig_SNPs$SNP_ID[snp_i])
  snp_name = paste0(sig_SNPs$chr[snp_i],"_",sig_SNPs$pos[snp_i])
  sigSNP_LD_plot_i = ggplot(LD_plot_i) + geom_smooth(aes(x=abs(dist),y=r2,group=line_gr,colour=line_gr),method="gam",se=FALSE) + ggtitle(paste0("LD of SNPs around ",snp_name)) + geom_point(aes(x=abs(dist),y=r2,group=plot_gr,colour=plot_gr)) + xlab("absolute distance on chromosome (bp)") + ylab("LD (r squared)")
  ggsave(filename = paste0("results/",species,"/sigSNP_asym/LD/",daynumber,"_",poolName,"/",poolName,"_LD_around_",snp_name,"_small.png"),plot=sigSNP_LD_plot_i,width=4000,height=3000,units="px",device="png")
  rm(sigSNP_LD_plot_i)
  sigSNP_LD_plot_i = ggplot(LD_plot_i) + geom_line(aes(x=dist,y=r2,group=line_gr,colour=line_gr)) + ggtitle(paste0("LD of SNPs around ",snp_name)) + geom_point(aes(x=dist,y=r2,group=plot_gr,colour=plot_gr)) + xlab("absolute distance on chromosome (bp)") + ylab("LD (r squared)")
  ggsave(filename = paste0("results/",species,"/sigSNP_asym/LD/",daynumber,"_",poolName,"/",poolName,"_LD_around_",snp_name,"_twosided_small.png"),plot=sigSNP_LD_plot_i,width=4000,height=3000,units="px",device="png")
  rm(sigSNP_LD_plot_i)
}
rm(sigSNP_LD_plot)
rm(sigSNP_LD_log_plot)

#clean
rm(haplo)
rm(gwas_results)
rm(gwas_results_i)
rm(nSNPs)
rm(p_plot)
rm(focal_SNPs)
rm(dist_threshold_big)
rm(dist_threshold_medium)
rm(dist_threshold_small)
rm(p_plot_i)
rm(LD_plot)
rm(snp_i)
rm(LD_plot_i)
rm(snp_name)
detach(package:ggplot2)

