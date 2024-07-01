library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)

#this script generates an additional plot for effect vs -log(p)

#run after 5_effect+selfGWAS

################################# settings
#you need to adjust this before running

#day number of effect output
daynumber_effect = "20240429"

#pool settings:
species = "barley"
#barley or wheat
poolName = "PM_2015"
#barley pools:
#NFNB_pooled, PM_2015, Blr_pooled, Scald_pooled, SFNB_pooled
#wheat pools:
#YR_pooled, LR_2017, SR_2017, YLS_2017, STB_2019, YS_2019

#################################

GWAS_results_effect = read.csv(paste0("results/",species,"/effect/",daynumber_effect,"_",poolName,"/effect_vs_dist/",poolName,"_neiGWAS_results_effect_and_p_all_dist_with_bases.csv"))
sig_SNPs = read.csv(paste0("results/",species,"/effect/",daynumber_effect,"_",poolName,"/",poolName,"_sig_SNPs_effect.csv"))

sig_SNPs = filter(GWAS_results_effect, pos %in% sig_SNPs$pos)
sig_SNPs_p = sig_SNPs[,c(1,2,((3:(ncol(sig_SNPs)/3))*3))]
sig_SNPs_b = sig_SNPs[,c(1,2,((3:(ncol(sig_SNPs)/3))*3)-2)]

sig_SNPs_p = mutate(sig_SNPs_p, "chr_pos" = paste0(chr,"_",pos))
sig_SNPs_b = mutate(sig_SNPs_b, "chr_pos" = paste0(chr,"_",pos))
sig_SNPs_p = sig_SNPs_p[,c(ncol(sig_SNPs_p),3:(ncol(sig_SNPs_p)-1))]
sig_SNPs_b = sig_SNPs_b[,c(ncol(sig_SNPs_b),3:(ncol(sig_SNPs_b)-1))]

sig_SNPs_p = gather(sig_SNPs_p, key = "distance", value = "p_log10", -chr_pos)
sig_SNPs_p = mutate(sig_SNPs_p, distance = sapply(strsplit(distance,split="_",fixed=TRUE), function(x) (x[3])))
sig_SNPs_b = gather(sig_SNPs_b, key = "distance", value = "b", -chr_pos)
sig_SNPs_b = mutate(sig_SNPs_b, distance = sapply(strsplit(distance,split="_",fixed=TRUE), function(x) (x[2])))

sig_SNPs = merge(sig_SNPs_p, sig_SNPs_b)
sig_SNPs = arrange(sig_SNPs, as.numeric(distance))
sig_SNPs = arrange(sig_SNPs, chr_pos)

#todo: plot
effect_plot = ggplot(sig_SNPs, aes(x=p_log10,y=b*2,group=chr_pos,colour=chr_pos)) + geom_path() + geom_hline(yintercept = 0) + geom_vline(xintercept = -log10(0.05/nrow(GWAS_results_effect))) + ggtitle("effects of base difference in sig. SNPs on susceptibility vs. -log(p)") + xlab("-log10(p)") + ylab("effect on susceptibilty (1x for different base, 0,5x for heterozygous)")
ggsave(filename = paste0("results/",species,"/effect/",daynumber_effect,"_",poolName,"/effect_vs_dist/",poolName,"_sigSNPs_effect_vs_p.png"),plot=effect_plot,width=4000,height=3000,units="px",device="png")
rm(effect_plot)

