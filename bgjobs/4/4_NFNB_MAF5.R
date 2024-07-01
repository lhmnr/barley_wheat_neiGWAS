library(rNeighborGWAS)
library(gaston)
library(dplyr)
library(readr)

#this script is to extract the significant SNPs and plot the LD to their surrounding,
#using haploid genotypes for better LD values,
#to decide which SNPs belong to a common loci and 
#get an idea in which ranges to look for candidate genes

#note: this and all further scripts should only be run for "pheno"s which have at least one significant SNP

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

#limit plot size to interesting area
#could be adjusted, but first one needs to be biggest
dist_threshold_big = 100000000
dist_threshold_medium = 5000000
dist_threshold_small = 250000
#recommendation: 100Mb, 5Mb, 250kb

################################# loading and preparing data
#before running check manually in folders if all required files are present

haplo = read.csv(paste0("results/",species,"/filtered_inputs/data_collection/haplo_",poolName,".csv"))

################################# extract significant SNPs

#importing data
GWAS_results = read.csv(paste0("results/",species,"/neiGWAS/data_collection/neiGWAS_results_",poolName,".csv"))
GWAS_results = mutate(GWAS_results, phenotype = poolName)
GWAS_results = select(GWAS_results, c(ncol(GWAS_results),1,2,3:(ncol(GWAS_results)-1)))
nSNPs = nrow(GWAS_results)

#filtering for significant SNPs
sig = {}
top_p = {}
top_dist = {}
for(j in 1:nrow(GWAS_results)){
  sig_j = FALSE
  top_p_j = 0
  top_dist_j = "NA"
  for(i in (1:((ncol(GWAS_results)-3)/2))*2+3){
    if(GWAS_results[j,i]>-log10(0.05/nSNPs)){
      sig_j = TRUE
      if(GWAS_results[j,i]>top_p_j){
        top_p_j = GWAS_results[j,i]
        top_dist_j = unlist(strsplit(colnames(GWAS_results)[i],split="_",fixed=TRUE))[3]
      }
    }
  }
  top_p = c(top_p, top_p_j)
  top_dist = c(top_dist, top_dist_j)
  if(sig_j){
    sig = c(sig,j)
  }
}
GWAS_results = mutate(GWAS_results, top_dist = top_dist)
GWAS_results = mutate(GWAS_results, top_p = top_p)
GWAS_results_sig = slice(GWAS_results, sig)

GWAS_results_sig = select(GWAS_results_sig, c(1,2,3,ncol(GWAS_results_sig)-1,ncol(GWAS_results_sig),4:ncol(GWAS_results_sig)-2))
GWAS_results_sig = arrange(GWAS_results_sig, -top_p)
sig_SNPs = select(GWAS_results_sig, c(1:5))

rm(top_p)
rm(top_p_j)
rm(top_dist)
rm(top_dist_j)
rm(sig)
rm(nSNPs)
rm(GWAS_results_sig)
rm(sig_j)

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
if(!file.exists(paste0("results/",species,"/sigSNP"))){
  dir.create(path = paste0("results/",species,"/sigSNP"))
}
if(!file.exists(paste0("results/",species,"/sigSNP/",daynumber,"_",poolName))){
  dir.create(path = paste0("results/",species,"/sigSNP/",daynumber,"_",poolName))
}
if(!file.exists(paste0("results/",species,"/sigSNP/",daynumber,"_",poolName,"/LD"))){
  dir.create(path = paste0("results/",species,"/sigSNP/",daynumber,"_",poolName,"/LD"))
}
if(!file.exists(paste0("results/",species,"/sigSNP/",daynumber,"_",poolName,"/p"))){
  dir.create(path = paste0("results/",species,"/sigSNP/",daynumber,"_",poolName,"/p"))
}
if(!file.exists(paste0("results/",species,"/sigSNP/data_collection"))){
  dir.create(path = paste0("results/",species,"/sigSNP/data_collection"))
}

write_csv(sig_SNPs, file = paste0("results/",species,"/sigSNP/",daynumber,"_",poolName,"/",poolName,"_sig_SNPs.csv"))
#copy sig_SNPs file manually to "results/[species]/sigSNP/data_collection" to use it for further analysis

write_csv(LD, file = paste0("results/",species,"/sigSNP/",daynumber,"_",poolName,"/LD/",poolName,"_SNPwise_LD_of_sig_SNPs.csv"))

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

write_csv(sig_SNPs, file = paste0("results/",species,"/sigSNP/",daynumber,"_",poolName,"/LD/",poolName,"_sig_SNPs_with_inter_LD.csv"))

################################# plot -log(p) and LD

library(ggplot2)

#prepare -log(p)-plots
p_plot = data.frame("p_log10"=numeric(),"SNP_foc"=character(),"SNP_comp"=character(),"dist"=numeric(),"nei_dist"=character())
for(snp_i in 1:nrow(sig_SNPs)){
  GWAS_results_i = select(GWAS_results, c("phenotype","chr","pos",paste0("p_log10_",sig_SNPs[snp_i,4])))
  colnames(GWAS_results_i)<- c("phenotype","chr","pos","p_log10")
  GWAS_results_i = mutate(GWAS_results_i, SNP_foc = paste0(sig_SNPs[snp_i,2],"_",sig_SNPs[snp_i,3]))
  GWAS_results_i = mutate(GWAS_results_i, SNP_comp = paste0(chr,"_",pos))
  GWAS_results_i = mutate(GWAS_results_i, dist = pos-sig_SNPs[snp_i,3])
  GWAS_results_i = filter(GWAS_results_i, abs(dist)<dist_threshold_big)
  GWAS_results_i = mutate(GWAS_results_i, nei_dist = sig_SNPs[snp_i,4])
  GWAS_results_i = select(GWAS_results_i, c(4,5,6,7,8))
  p_plot = rows_append(p_plot, GWAS_results_i)
}
focal_SNPs = unique(p_plot$SNP_foc)
#plot -log(p) arround sig SNPs, big
sigSNP_p_plot = ggplot(p_plot) + geom_line(aes(x=dist,y=p_log10,group=SNP_foc,colour=SNP_foc)) + ggtitle("-log10(p) of SNPs around sig. SNPs") + geom_point(aes(x=dist,y=p_log10,group=SNP_foc,colour=SNP_foc)) + xlab("distance on chromosome (bp)") + ylab("-log10(p)") + geom_hline(yintercept = -log10(0.05/nrow(GWAS_results)))
ggsave(filename = paste0("results/",species,"/sigSNP/",daynumber,"_",poolName,"/p/",poolName,"_-log10_p_around_sig_SNPs_big.png"),plot=sigSNP_p_plot,width=4000,height=3000,units="px",device="png")
#plot -log(p) per SNP
for(snp_i in 1:nrow(sig_SNPs)){
  p_plot_i = filter(p_plot, SNP_foc == focal_SNPs[snp_i])
  snp_name = focal_SNPs[snp_i]
  sigSNP_p_plot_i = ggplot(p_plot_i) + geom_smooth(aes(x=abs(dist),y=p_log10,group=SNP_foc,colour=SNP_foc),method="gam",se=FALSE) + ggtitle(paste0("-log10(p) of SNPs around ",snp_name)) + geom_point(aes(x=abs(dist),y=p_log10,group=SNP_foc,colour=SNP_foc)) + xlab("absolute distance on chromosome (bp)") + ylab("-log10(p)") + geom_hline(yintercept = -log10(0.05/nrow(GWAS_results)))
  ggsave(filename = paste0("results/",species,"/sigSNP/",daynumber,"_",poolName,"/p/",poolName,"_-log10_p_around_",snp_name,"_onesided_big.png"),plot=sigSNP_p_plot_i,width=4000,height=3000,units="px",device="png")
  rm(sigSNP_p_plot_i)
  sigSNP_p_plot_i = ggplot(p_plot_i) + geom_line(aes(x=dist,y=p_log10,group=SNP_foc,colour=SNP_foc)) + ggtitle(paste0("-log10(p) of SNPs around ",snp_name)) + geom_point(aes(x=dist,y=p_log10,group=SNP_foc,colour=SNP_foc)) + xlab("absolute distance on chromosome (bp)") + ylab("-log10(p)") + geom_hline(yintercept = -log10(0.05/nrow(GWAS_results)))
  ggsave(filename = paste0("results/",species,"/sigSNP/",daynumber,"_",poolName,"/p/",poolName,"_-log10_p_around_",snp_name,"_big.png"),plot=sigSNP_p_plot_i,width=4000,height=3000,units="px",device="png")
  rm(sigSNP_p_plot_i)
}
rm(sigSNP_p_plot)
p_plot = filter(p_plot, abs(dist)<dist_threshold_medium)
#plot -log(p) arround sig SNPs, medium
sigSNP_p_plot = ggplot(p_plot) + geom_line(aes(x=dist,y=p_log10,group=SNP_foc,colour=SNP_foc)) + ggtitle("-log10(p) of SNPs around sig. SNPs") + geom_point(aes(x=dist,y=p_log10,group=SNP_foc,colour=SNP_foc)) + xlab("distance on chromosome (bp)") + ylab("-log10(p)") + geom_hline(yintercept = -log10(0.05/nrow(GWAS_results)))
ggsave(filename = paste0("results/",species,"/sigSNP/",daynumber,"_",poolName,"/p/",poolName,"_-log10_p_around_sig_SNPs_medium.png"),plot=sigSNP_p_plot,width=4000,height=3000,units="px",device="png")
#plot -log(p) per SNP
for(snp_i in 1:nrow(sig_SNPs)){
  p_plot_i = filter(p_plot, SNP_foc == focal_SNPs[snp_i])
  snp_name = focal_SNPs[snp_i]
  sigSNP_p_plot_i = ggplot(p_plot_i) + geom_smooth(aes(x=abs(dist),y=p_log10,group=SNP_foc,colour=SNP_foc),method="gam",se=FALSE) + ggtitle(paste0("-log10(p) of SNPs around ",snp_name)) + geom_point(aes(x=abs(dist),y=p_log10,group=SNP_foc,colour=SNP_foc)) + xlab("absolute distance on chromosome (bp)") + ylab("-log10(p)") + geom_hline(yintercept = -log10(0.05/nrow(GWAS_results)))
  ggsave(filename = paste0("results/",species,"/sigSNP/",daynumber,"_",poolName,"/p/",poolName,"_-log10_p_around_",snp_name,"_onesided_medium.png"),plot=sigSNP_p_plot_i,width=4000,height=3000,units="px",device="png")
  rm(sigSNP_p_plot_i)
  sigSNP_p_plot_i = ggplot(p_plot_i) + geom_line(aes(x=dist,y=p_log10,group=SNP_foc,colour=SNP_foc)) + ggtitle(paste0("-log10(p) of SNPs around ",snp_name)) + geom_point(aes(x=dist,y=p_log10,group=SNP_foc,colour=SNP_foc)) + xlab("absolute distance on chromosome (bp)") + ylab("-log10(p)") + geom_hline(yintercept = -log10(0.05/nrow(GWAS_results)))
  ggsave(filename = paste0("results/",species,"/sigSNP/",daynumber,"_",poolName,"/p/",poolName,"_-log10_p_around_",snp_name,"_medium.png"),plot=sigSNP_p_plot_i,width=4000,height=3000,units="px",device="png")
  rm(sigSNP_p_plot_i)
}
rm(sigSNP_p_plot)
p_plot = filter(p_plot, abs(dist)<dist_threshold_small)
#plot -log(p) arround sig SNPs, small
sigSNP_p_plot = ggplot(p_plot) + geom_line(aes(x=dist,y=p_log10,group=SNP_foc,colour=SNP_foc)) + ggtitle("-log10(p) of SNPs around sig. SNPs") + geom_point(aes(x=dist,y=p_log10,group=SNP_foc,colour=SNP_foc)) + xlab("distance on chromosome (bp)") + ylab("-log10(p)") + geom_hline(yintercept = -log10(0.05/nrow(GWAS_results)))
ggsave(filename = paste0("results/",species,"/sigSNP/",daynumber,"_",poolName,"/p/",poolName,"_-log10_p_around_sig_SNPs_small.png"),plot=sigSNP_p_plot,width=4000,height=3000,units="px",device="png")
#plot -log(p) per SNP
for(snp_i in 1:nrow(sig_SNPs)){
  p_plot_i = filter(p_plot, SNP_foc == focal_SNPs[snp_i])
  snp_name = focal_SNPs[snp_i]
  sigSNP_p_plot_i = ggplot(p_plot_i) + geom_smooth(aes(x=abs(dist),y=p_log10,group=SNP_foc,colour=SNP_foc),method="gam",se=FALSE) + ggtitle(paste0("-log10(p) of SNPs around ",snp_name)) + geom_point(aes(x=abs(dist),y=p_log10,group=SNP_foc,colour=SNP_foc)) + xlab("absolute distance on chromosome (bp)") + ylab("-log10(p)") + geom_hline(yintercept = -log10(0.05/nrow(GWAS_results)))
  ggsave(filename = paste0("results/",species,"/sigSNP/",daynumber,"_",poolName,"/p/",poolName,"_-log10_p_around_",snp_name,"_onesided_small.png"),plot=sigSNP_p_plot_i,width=4000,height=3000,units="px",device="png")
  rm(sigSNP_p_plot_i)
  sigSNP_p_plot_i = ggplot(p_plot_i) + geom_line(aes(x=dist,y=p_log10,group=SNP_foc,colour=SNP_foc)) + ggtitle(paste0("-log10(p) of SNPs around ",snp_name)) + geom_point(aes(x=dist,y=p_log10,group=SNP_foc,colour=SNP_foc)) + xlab("absolute distance on chromosome (bp)") + ylab("-log10(p)") + geom_hline(yintercept = -log10(0.05/nrow(GWAS_results)))
  ggsave(filename = paste0("results/",species,"/sigSNP/",daynumber,"_",poolName,"/p/",poolName,"_-log10_p_around_",snp_name,"_small.png"),plot=sigSNP_p_plot_i,width=4000,height=3000,units="px",device="png")
  rm(sigSNP_p_plot_i)
}
rm(sigSNP_p_plot)

#prepare LD plots
LD_plot = filter(LD, same_chr&abs(dist)<dist_threshold_big)
LD_plot = mutate(LD_plot, plot_gr = ifelse(SNP_comp_sig,paste0(sapply(SNP_comp, function(x) sig_SNPs[which(sig_SNPs$SNP_ID==x),2]),"_",sapply(SNP_comp, function(x) sig_SNPs[which(sig_SNPs$SNP_ID==x),3])),"others"))
LD_plot = mutate(LD_plot, line_gr = paste0(sapply(SNP_foc, function(x) sig_SNPs[which(sig_SNPs$SNP_ID==x),2]),"_",sapply(SNP_foc, function(x) sig_SNPs[which(sig_SNPs$SNP_ID==x),3])))
#plot LD arround sig SNPs, big
sigSNP_LD_plot = ggplot(LD_plot) + geom_line(aes(x=dist,y=r2,group=line_gr,colour=line_gr)) + ggtitle("LD of SNPs around sig. SNPs") + geom_point(aes(x=dist,y=r2,group=plot_gr,colour=plot_gr)) + xlab("distance on chromosome (bp)") + ylab("LD")
ggsave(filename = paste0("results/",species,"/sigSNP/",daynumber,"_",poolName,"/LD/",poolName,"_LD_around_sig_SNPs_big.png"),plot=sigSNP_LD_plot,width=4000,height=3000,units="px",device="png")
#plot with logarithmic y axis
sigSNP_LD_log_plot = ggplot(LD_plot) + geom_line(aes(x=dist,y=log10(r2),group=line_gr,colour=line_gr)) + ggtitle("LD of SNPs around sig. SNPs") + geom_point(aes(x=dist,y=log10(r2),group=plot_gr,colour=plot_gr)) + xlab("distance on chromosome (bp)") + ylab("log10 LD")
ggsave(filename = paste0("results/",species,"/sigSNP/",daynumber,"_",poolName,"/LD/",poolName,"_log_LD_around_sig_SNPs_big.png"),plot=sigSNP_LD_log_plot,width=4000,height=3000,units="px",device="png")
#plot one-sided and two-sided LD decay per SNP
for(snp_i in 1:nrow(sig_SNPs)){
  LD_plot_i = filter(LD_plot, SNP_foc == sig_SNPs$SNP_ID[snp_i])
  snp_name = paste0(sig_SNPs$chr[snp_i],"_",sig_SNPs$pos[snp_i])
  sigSNP_LD_plot_i = ggplot(LD_plot_i) + geom_smooth(aes(x=abs(dist),y=r2,group=line_gr,colour=line_gr),method="gam",se=FALSE) + ggtitle(paste0("LD of SNPs around ",snp_name)) + geom_point(aes(x=abs(dist),y=r2,group=plot_gr,colour=plot_gr)) + xlab("absolute distance on chromosome (bp)") + ylab("LD (r squared)")
  ggsave(filename = paste0("results/",species,"/sigSNP/",daynumber,"_",poolName,"/LD/",poolName,"_LD_around_",snp_name,"_big.png"),plot=sigSNP_LD_plot_i,width=4000,height=3000,units="px",device="png")
  rm(sigSNP_LD_plot_i)
  sigSNP_LD_plot_i = ggplot(LD_plot_i) + geom_line(aes(x=dist,y=r2,group=line_gr,colour=line_gr)) + ggtitle(paste0("LD of SNPs around ",snp_name)) + geom_point(aes(x=dist,y=r2,group=plot_gr,colour=plot_gr)) + xlab("absolute distance on chromosome (bp)") + ylab("LD (r squared)")
  ggsave(filename = paste0("results/",species,"/sigSNP/",daynumber,"_",poolName,"/LD/",poolName,"_LD_around_",snp_name,"_twosided_big.png"),plot=sigSNP_LD_plot_i,width=4000,height=3000,units="px",device="png")
  rm(sigSNP_LD_plot_i)
}
rm(sigSNP_LD_plot)
rm(sigSNP_LD_log_plot)
LD_plot = filter(LD_plot, abs(dist)<dist_threshold_medium)
#plot LD arround sig SNPs, medium
sigSNP_LD_plot = ggplot(LD_plot) + geom_line(aes(x=dist,y=r2,group=line_gr,colour=line_gr)) + ggtitle("LD of SNPs around sig. SNPs") + geom_point(aes(x=dist,y=r2,group=plot_gr,colour=plot_gr)) + xlab("distance on chromosome (bp)") + ylab("LD")
ggsave(filename = paste0("results/",species,"/sigSNP/",daynumber,"_",poolName,"/LD/",poolName,"_LD_around_sig_SNPs_medium.png"),plot=sigSNP_LD_plot,width=4000,height=3000,units="px",device="png")
#plot with logarithmic y axis
sigSNP_LD_log_plot = ggplot(LD_plot) + geom_line(aes(x=dist,y=log10(r2),group=line_gr,colour=line_gr)) + ggtitle("LD of SNPs around sig. SNPs") + geom_point(aes(x=dist,y=log10(r2),group=plot_gr,colour=plot_gr)) + xlab("distance on chromosome (bp)") + ylab("log10 LD")
ggsave(filename = paste0("results/",species,"/sigSNP/",daynumber,"_",poolName,"/LD/",poolName,"_log_LD_around_sig_SNPs_medium.png"),plot=sigSNP_LD_log_plot,width=4000,height=3000,units="px",device="png")
#plot one-sided and two-sided LD decay per SNP
for(snp_i in 1:nrow(sig_SNPs)){
  LD_plot_i = filter(LD_plot, SNP_foc == sig_SNPs$SNP_ID[snp_i])
  snp_name = paste0(sig_SNPs$chr[snp_i],"_",sig_SNPs$pos[snp_i])
  sigSNP_LD_plot_i = ggplot(LD_plot_i) + geom_smooth(aes(x=abs(dist),y=r2,group=line_gr,colour=line_gr),method="gam",se=FALSE) + ggtitle(paste0("LD of SNPs around ",snp_name)) + geom_point(aes(x=abs(dist),y=r2,group=plot_gr,colour=plot_gr)) + xlab("absolute distance on chromosome (bp)") + ylab("LD (r squared)")
  ggsave(filename = paste0("results/",species,"/sigSNP/",daynumber,"_",poolName,"/LD/",poolName,"_LD_around_",snp_name,"_medium.png"),plot=sigSNP_LD_plot_i,width=4000,height=3000,units="px",device="png")
  rm(sigSNP_LD_plot_i)
  sigSNP_LD_plot_i = ggplot(LD_plot_i) + geom_line(aes(x=dist,y=r2,group=line_gr,colour=line_gr)) + ggtitle(paste0("LD of SNPs around ",snp_name)) + geom_point(aes(x=dist,y=r2,group=plot_gr,colour=plot_gr)) + xlab("absolute distance on chromosome (bp)") + ylab("LD (r squared)")
  ggsave(filename = paste0("results/",species,"/sigSNP/",daynumber,"_",poolName,"/LD/",poolName,"_LD_around_",snp_name,"_twosided_medium.png"),plot=sigSNP_LD_plot_i,width=4000,height=3000,units="px",device="png")
  rm(sigSNP_LD_plot_i)
}
rm(sigSNP_LD_plot)
rm(sigSNP_LD_log_plot)
LD_plot = filter(LD_plot, abs(dist)<dist_threshold_small)
#plot LD arround sig SNPs, small
sigSNP_LD_plot = ggplot(LD_plot) + geom_line(aes(x=dist,y=r2,group=line_gr,colour=line_gr)) + ggtitle("LD of SNPs around sig. SNPs") + geom_point(aes(x=dist,y=r2,group=plot_gr,colour=plot_gr)) + xlab("distance on chromosome (bp)") + ylab("LD")
ggsave(filename = paste0("results/",species,"/sigSNP/",daynumber,"_",poolName,"/LD/",poolName,"_LD_around_sig_SNPs_small.png"),plot=sigSNP_LD_plot,width=4000,height=3000,units="px",device="png")
#plot with logarithmic y axis
sigSNP_LD_log_plot = ggplot(LD_plot) + geom_line(aes(x=dist,y=log10(r2),group=line_gr,colour=line_gr)) + ggtitle("LD of SNPs around sig. SNPs") + geom_point(aes(x=dist,y=log10(r2),group=plot_gr,colour=plot_gr)) + xlab("distance on chromosome (bp)") + ylab("log10 LD")
ggsave(filename = paste0("results/",species,"/sigSNP/",daynumber,"_",poolName,"/LD/",poolName,"_log_LD_around_sig_SNPs_small.png"),plot=sigSNP_LD_log_plot,width=4000,height=3000,units="px",device="png")
#plot one-sided and two-sided LD decay per SNP
for(snp_i in 1:nrow(sig_SNPs)){
  LD_plot_i = filter(LD_plot, SNP_foc == sig_SNPs$SNP_ID[snp_i])
  snp_name = paste0(sig_SNPs$chr[snp_i],"_",sig_SNPs$pos[snp_i])
  sigSNP_LD_plot_i = ggplot(LD_plot_i) + geom_smooth(aes(x=abs(dist),y=r2,group=line_gr,colour=line_gr),method="gam",se=FALSE) + ggtitle(paste0("LD of SNPs around ",snp_name)) + geom_point(aes(x=abs(dist),y=r2,group=plot_gr,colour=plot_gr)) + xlab("absolute distance on chromosome (bp)") + ylab("LD (r squared)")
  ggsave(filename = paste0("results/",species,"/sigSNP/",daynumber,"_",poolName,"/LD/",poolName,"_LD_around_",snp_name,"_small.png"),plot=sigSNP_LD_plot_i,width=4000,height=3000,units="px",device="png")
  rm(sigSNP_LD_plot_i)
  sigSNP_LD_plot_i = ggplot(LD_plot_i) + geom_line(aes(x=dist,y=r2,group=line_gr,colour=line_gr)) + ggtitle(paste0("LD of SNPs around ",snp_name)) + geom_point(aes(x=dist,y=r2,group=plot_gr,colour=plot_gr)) + xlab("absolute distance on chromosome (bp)") + ylab("LD (r squared)")
  ggsave(filename = paste0("results/",species,"/sigSNP/",daynumber,"_",poolName,"/LD/",poolName,"_LD_around_",snp_name,"_twosided_small.png"),plot=sigSNP_LD_plot_i,width=4000,height=3000,units="px",device="png")
  rm(sigSNP_LD_plot_i)
}
rm(sigSNP_LD_plot)
rm(sigSNP_LD_log_plot)

#clean
rm(GWAS_results)
rm(GWAS_results_i)
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

