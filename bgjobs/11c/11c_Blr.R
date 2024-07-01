library(rNeighborGWAS)
library(gaston)
library(dplyr)
library(readr)

#this script is to obtain effect size and selfGWAS p-values,
#for all distances to 16, but only for sig. SNPs found up to 10,
#both of which require to run neiGWAS with extended output, thus they share this script

################################# settings
#you need to adjust this before running

#day number for easier identification of outputs
daynumber = "20240513"
#do not use slashes or points, since this will be used for the output folder name

#pool settings:
species = "barley"
#barley or wheat
poolName = "Blr_pooled"
#barley pools:
#NFNB_pooled, PM_2015, Blr_pooled, Scald_pooled, SFNB_pooled
#wheat pools:
#YR_pooled, LR_2017, SR_2017, YLS_2017, STB_2019, YS_2019

#cores to use for calc_PVEnei and neiGWAS
cores = 4L

################################# loading and preparing data
#before running check manually in folders if all required files are present

geno = read.csv(paste0("results/",species,"/filtered_inputs/data_collection/geno_",poolName,".csv"))
pheno = read.csv(paste0("results/",species,"/filtered_inputs/data_collection/pheno_",poolName,".csv"))
gmap = read.csv(paste0("results/",species,"/filtered_inputs/data_collection/gmap_",poolName,".csv"))
smap = read.csv(paste0("results/",species,"/filtered_inputs/data_collection/smap_",poolName,".csv"))
af = read.csv(paste0("results/",species,"/filtered_inputs/data_collection/af_",poolName,".csv"))

geno = as.matrix(geno)
colnames(geno)<-c(1:ncol(geno))
smap = as.matrix(smap)
af = af[,1]

# import significant SNPs

#note: sig_SNPs files have to be manually copied from folders into sigSNP/data_collection
sig_SNPs <- read_csv(paste0("results/",species,"/sigSNP/data_collection/",poolName,"_sig_SNPs.csv"))

################################# plot sig SNPs effects over all distances
#this can take a while...

#note: sig_SNPs files have to be manually copied from folders into sigSNP/data_collection
sig_SNPs <- read_csv(paste0("results/",species,"/sigSNP/data_collection/",poolName,"_sig_SNPs.csv"))

distances = c(1,1.5,2,2.3,2.9,3,3.2,3.7,4,4.2,4.3,4.5,5,5.2,5.4,5.7,5.9,6,6.1,6.4,6.5,6.8,7,7.1,7.25,7.3,7.7,7.9,8,8.1,8.3,8.5,8.6,8.7,8.98,9,9.1,9.3,9.45,9.5,9.9,9.95,10,10.1,10.2,10.3,10.5,10.7,10.8,10.9,11,11.1,11.2,11.4,11.5,11.7,11.8,12,12.05,12.1,12.2,12.3,12.4,12.6,12.7,12.8,12.9,13,13.1,13.2,13.4,13.45,13.5,13.7,13.9,13.95,14,14.1,14.2,14.3,14.4,14.5,14.6,14.8,14.9,15,15.1,15.2,15.25,15.275,15.3,15.55,15.6,15.65,15.7,15.9,16)

gwas_results = data.frame(NULL)
for(scale in distances){
  g_nei <- nei_coval(geno=geno,
                     smap=smap,
                     scale=scale,
                     alpha = Inf,
                     #kernel = c("exp", "gaussian"),
                     grouping = pheno$Experiment_Number,
                     n_core = cores
  )
  gwas_out <- nei_lmm(geno=geno,pheno=pheno$Damage_Level,
                      g_nei,
                      addcovar=as.factor(pheno$Experiment_Number),
                      #response = "quantitative",
                      n_core = cores,
                      #asym = FALSE
  )
  gwas_out = mutate(gwas_out, chr = gmap$CHROM)
  gwas_out = mutate(gwas_out, pos = gmap$POS)
  gwas_out = mutate(gwas_out, p_nei_log10 = -log10(p_nei))
  if(ncol(gwas_results)==0){
    gwas_out_0 = mutate(gwas_out, p_self_log10 = -log10(p_self))
    gwas_out_0 = select(gwas_out_0, c(5,6,3,8))
    colnames(gwas_out_0)<-c("chr","pos","p_0","p_log10_0")
    gwas_results = gwas_out_0
  }
  gwas_out = select(gwas_out, c(5,6,2,4,7))
  colnames(gwas_out)<-c("chr","pos",paste0("b_",toString(scale)),paste0("p_",toString(scale)),paste0("p_log10_",toString(scale)))
  gwas_results = merge(gwas_results,gwas_out)
}
gwas_results = arrange(gwas_results, pos)
gwas_results = arrange(gwas_results, chr)

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
gwas_results = gwas_results[,c(2,3,ncol(gwas_results)-1,ncol(gwas_results),4:(ncol(gwas_results)-2))]
gwas_results = arrange(gwas_results, pos)
gwas_results = arrange(gwas_results, chr)

#create output folder if necessary
if(!file.exists(paste0("results/",species,"/effect_extended"))){
  dir.create(path = paste0("results/",species,"/effect_extended"))
}
if(!file.exists(paste0("results/",species,"/effect_extended/",daynumber,"_",poolName))){
  dir.create(path = paste0("results/",species,"/effect_extended/",daynumber,"_",poolName))
}
if(!file.exists(paste0("results/",species,"/effect_extended/",daynumber,"_",poolName,"/effect_vs_dist"))){
  dir.create(path = paste0("results/",species,"/effect_extended/",daynumber,"_",poolName,"/effect_vs_dist"))
}

write_csv(gwas_results, file = paste0("results/",species,"/effect_extended/",daynumber,"_",poolName,"/effect_vs_dist/",poolName,"_neiGWAS_results_effect_and_p_all_dist_with_bases.csv"))

sig_SNPs_effect = sig_SNPs[,1:6]
gwas_results_sig = slice(gwas_results, sig_SNPs_effect$SNP_ID)
sig_SNPs_effect = merge(sig_SNPs_effect, gwas_results_sig)
sig_SNPs_effect = mutate(sig_SNPs_effect, SNP_info = paste0(chr,"_",pos,"_(",base_ref,"_or_",base_alt,")"))
line_selection = (3*(3:((ncol(sig_SNPs_effect)-3)/3)))+2
line_selection_2 = (3*(2:((ncol(sig_SNPs_effect)-3)/3)))+4
sig_SNPs_p = select(sig_SNPs_effect, c(ncol(sig_SNPs_effect),line_selection_2))
sig_SNPs_effect = select(sig_SNPs_effect, c(ncol(sig_SNPs_effect),line_selection))

library(tidyr)
sig_SNPs_effect = gather(sig_SNPs_effect, key = "distance", value = "b", -SNP_info)
sig_SNPs_p = gather(sig_SNPs_p, key = "distance", value = "p_log10", -SNP_info)
sig_SNPs_effect = mutate(sig_SNPs_effect, distance = sapply(strsplit(distance,split="_",fixed=TRUE), function(x) (x[2])))
sig_SNPs_p = mutate(sig_SNPs_p, distance = sapply(strsplit(distance,split="_",fixed=TRUE), function(x) (x[3])))
sig_SNPs_effect = mutate(sig_SNPs_effect, b = 2*b)
detach(package:tidyr)

# would be nice to plot effect and p together

library(ggplot2)

#discrete plot
x_axis_order = as.factor(sort(unique(as.numeric(sig_SNPs_p$distance))))
sigSNPs_effect_disc = ggplot(sig_SNPs_effect, aes(x=distance,y=b,group=SNP_info,colour=SNP_info)) + geom_line() + geom_hline(yintercept = 0) + ggtitle("effects of base difference in sig. SNPs on suseptibility for each distance") + scale_x_discrete(limits = x_axis_order) + xlab("distance of included neighbours") + ylab("effect on suseptibilty (1x for different base, 0,5x for heterozygous)")
sigSNPs_p_disc = ggplot(sig_SNPs_p, aes(x=distance,y=p_log10,group=SNP_info,colour=SNP_info)) + geom_line() + geom_hline(yintercept = -log10(0.05/nrow(gwas_results))) + ggtitle("p-values of effects of base difference in sig. SNPs on suseptibility for each distance") + scale_x_discrete(limits = x_axis_order) + xlab("distance of included neighbours (0 means self-GWAS)") + ylab("-log10 p-value")
ggsave(filename = paste0("results/",species,"/effect_extended/",daynumber,"_",poolName,"/effect_vs_dist/",poolName,"_sigSNPs_effect_vs_distance_discrete.png"),plot=sigSNPs_effect_disc,width=4000,height=3000,units="px",device="png")
ggsave(filename = paste0("results/",species,"/effect_extended/",daynumber,"_",poolName,"/effect_vs_dist/",poolName,"_sigSNPs_p_vs_distance_discrete.png"),plot=sigSNPs_p_disc,width=4000,height=3000,units="px",device="png")
rm(sigSNPs_effect_disc)
rm(sigSNPs_p_disc)

#continuous plot
sig_SNPs_effect = mutate(sig_SNPs_effect, distance = as.numeric(distance))
sig_SNPs_p = mutate(sig_SNPs_p, distance = as.numeric(distance))
sigSNPs_effect_cont = ggplot(sig_SNPs_effect, aes(x=distance,y=b,group=SNP_info,colour=SNP_info)) + geom_line() + geom_hline(yintercept = 0) + ggtitle("effects of base difference in sig. SNPs on suseptibility for each distance") + xlab("distance of included neighbours") + ylab("effect on suseptibilty (1x for different base, 0,5x for heterozygous)")
sigSNPs_p_cont = ggplot(sig_SNPs_p, aes(x=distance,y=p_log10,group=SNP_info,colour=SNP_info)) + geom_line() + geom_hline(yintercept = -log10(0.05/nrow(gwas_results))) + ggtitle("p-values of effects of base difference in sig. SNPs on suseptibility for each distance") + xlab("distance of included neighbours (0 means self-GWAS)") + ylab("-log10 p-value")
ggsave(filename = paste0("results/",species,"/effect_extended/",daynumber,"_",poolName,"/effect_vs_dist/",poolName,"_sigSNPs_effect_vs_distance_continuous.png"),plot=sigSNPs_effect_cont,width=4000,height=3000,units="px",device="png")
ggsave(filename = paste0("results/",species,"/effect_extended/",daynumber,"_",poolName,"/effect_vs_dist/",poolName,"_sigSNPs_p_vs_distance_continuous.png"),plot=sigSNPs_p_cont,width=4000,height=3000,units="px",device="png")
rm(sigSNPs_effect_cont)
rm(sigSNPs_p_cont)

#cleaning up
remove(g_nei)
remove(gwas_out_0)
remove(gwas_out)
remove(scale)
remove(distances)
remove(sig_SNPs)
remove(sig_SNPs_effect)
remove(sig_SNPs_p)
remove(gwas_results)
remove(gwas_results_sig)
remove(genotypes)
remove(line_selection)
remove(line_selection_2)
remove(x_axis_order)
detach(package:ggplot2)

