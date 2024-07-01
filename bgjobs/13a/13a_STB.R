library(rNeighborGWAS)
library(gaston)
library(dplyr)
library(readr)

#this script is to do neiGWAS over all distances for a given pheno
#using the "sxn" interaction term

#note: from this script on, barley and wheat share scripts

################################# settings
#you need to adjust this before running

#day number for easier identification of outputs
daynumber = "20240429"
#do not use slashes or points, since this will be used for the output folder name

#pool settings:
species = "wheat"
#barley or wheat
poolName = "STB_2019"
#barley pools:
#NFNB_pooled, PM_2015, Blr_pooled, Scald_pooled, SFNB_pooled
#wheat pools:
#YR_pooled, LR_2017, SR_2017, YLS_2017, STB_2019, YS_2019

#cores to use for calc_PVEnei and neiGWAS
cores = 1L

#plot settings:
#desired plots for auto neiGWAS, will be printed into files in output folder named after poolName
wished_plots = c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE)
#plot 1: PVEnei vs. distance plot, taken from delta_PVE-function
#plot 2: p-value vs. distance plot for top-SNPs, continuous scale until max. distance
#plot 3: p-value vs. distance plot for top-SNPs, continuous scale for distances < 6
#plot 4: p-value vs. distance plot for top-SNPs, discrete scale until max. distance
#plot 5: p-value vs. distance plot for top-SNPs, discrete scale for distances < 6
#plots 6: mathattan and qq plots for all distances, may gives many plots

################################# loading and preparing data
#before running check manually in folders if all required files are present

geno = read.csv(paste0("results/",species,"/filtered_inputs/data_collection/geno_",poolName,".csv"))
pheno = read.csv(paste0("results/",species,"/filtered_inputs/data_collection/pheno_",poolName,".csv"))
gmap = read.csv(paste0("results/",species,"/filtered_inputs/data_collection/gmap_",poolName,".csv"))
smap = read.csv(paste0("results/",species,"/filtered_inputs/data_collection/smap_",poolName,".csv"))
af = read.csv(paste0("results/",species,"/filtered_inputs/data_collection/af_",poolName,".csv"))

geno = as.matrix(geno)
smap = as.matrix(smap)
af = af[,1]

################################# automatic neiGWAS
# this can take a while...

gwas_results = data.frame(NULL)
dist_results = data.frame(NULL)
top_SNPs_self = data.frame(NULL)
top_SNPs_nei = data.frame(NULL)
top_SNPs_sxn = data.frame(NULL)

#distances up to 6:
#distances = c(1,1.5,2,2.3,2.9,3,3.2,3.7,4,4.2,4.3,4.5,5,5.2,5.4,5.7,5.9,6)
#distances up to 10:
distances = c(1,1.5,2,2.3,2.9,3,3.2,3.7,4,4.2,4.3,4.5,5,5.2,5.4,5.7,5.9,6,6.1,6.4,6.5,6.8,7,7.1,7.25,7.3,7.7,7.9,8,8.1,8.3,8.5,8.6,8.7,8.98,9,9.1,9.3,9.45,9.5,9.9,9.95,10)

#it seems calc_PVEnei does not include asymetric
dist_results = calc_PVEnei(geno=geno,pheno=pheno$Damage_Level,smap=smap,scale=distances,addcovar=as.factor(pheno$Experiment_Number),grouping=pheno$Experiment_Number,n_core=cores)

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
                      asym = TRUE
  )
  gwas_out = mutate(gwas_out, chr = gmap$CHROM)
  gwas_out = mutate(gwas_out, pos = gmap$POS)
  gwas_out = mutate(gwas_out, p_self_log10 = -log10(p_self))
  gwas_out = mutate(gwas_out, p_nei_log10 = -log10(p_nei))
  gwas_out = mutate(gwas_out, p_sxn_log10 = -log10(p_sxn))
  gwas_out = select(gwas_out, c(7,8,4,9,1,5,10,2,6,11,3))
  colnames(gwas_out)<-c("chr","pos",
                        paste0("p_self_",toString(scale)),paste0("p_self_log10_",toString(scale)),paste0("b_self_",toString(scale)),
                        paste0("p_nei_",toString(scale)),paste0("p_nei_log10_",toString(scale)),paste0("b_nei_",toString(scale)),
                        paste0("p_sxn_",toString(scale)),paste0("p_sxn_log10_",toString(scale)),paste0("b_sxn_",toString(scale)))
  if(ncol(gwas_results)==0){
    gwas_results = gwas_out
  }
  else{
    gwas_results = merge(gwas_results,gwas_out)
  }
  gwas_out = arrange(gwas_out, gwas_out[3])
  topsnp_self = gwas_out[1,]
  topsnp_self = append(topsnp_self, scale)
  if(ncol(top_SNPs_self)==0){
    top_SNPs_self = data.frame(topsnp_self)
    colnames(top_SNPs_self)<-c("chr","pos","p_self","p_self_log10","b_self","p_nei","p_nei_log10","b_nei","p_sxn","p_sxn_log10","b_sxn","dist")
  }
  else{
    topsnp_self = data.frame(topsnp_self)
    colnames(topsnp_self)<-c("chr","pos","p_self","p_self_log10","b_self","p_nei","p_nei_log10","b_nei","p_sxn","p_sxn_log10","b_sxn","dist")
    top_SNPs_self = add_row(top_SNPs_self, topsnp_self)
  }
  gwas_out = arrange(gwas_out, gwas_out[6])
  topsnp_nei = gwas_out[1,]
  topsnp_nei = append(topsnp_nei, scale)
  if(ncol(top_SNPs_nei)==0){
    top_SNPs_nei = data.frame(topsnp_nei)
    colnames(top_SNPs_nei)<-c("chr","pos","p_self","p_self_log10","b_self","p_nei","p_nei_log10","b_nei","p_sxn","p_sxn_log10","b_sxn","dist")
  }
  else{
    topsnp_nei = data.frame(topsnp_nei)
    colnames(topsnp_nei)<-c("chr","pos","p_self","p_self_log10","b_self","p_nei","p_nei_log10","b_nei","p_sxn","p_sxn_log10","b_sxn","dist")
    top_SNPs_nei = add_row(top_SNPs_nei, topsnp_nei)
  }
  gwas_out = arrange(gwas_out, gwas_out[9])
  topsnp_sxn = gwas_out[1,]
  topsnp_sxn = append(topsnp_sxn, scale)
  if(ncol(top_SNPs_sxn)==0){
    top_SNPs_sxn = data.frame(topsnp_sxn)
    colnames(top_SNPs_sxn)<-c("chr","pos","p_self","p_self_log10","b_self","p_nei","p_nei_log10","b_nei","p_sxn","p_sxn_log10","b_sxn","dist")
  }
  else{
    topsnp_sxn = data.frame(topsnp_sxn)
    colnames(topsnp_sxn)<-c("chr","pos","p_self","p_self_log10","b_self","p_nei","p_nei_log10","b_nei","p_sxn","p_sxn_log10","b_sxn","dist")
    top_SNPs_sxn = add_row(top_SNPs_sxn, topsnp_sxn)
  }
}

#create output folder if necessary
if(!file.exists(paste0("results/",species,"/neiGWAS_asym"))){
  dir.create(path = paste0("results/",species,"/neiGWAS_asym"))
}
if(!file.exists(paste0("results/",species,"/neiGWAS_asym/",daynumber,"_",poolName))){
  dir.create(path = paste0("results/",species,"/neiGWAS_asym/",daynumber,"_",poolName))
}
if(!file.exists(paste0("results/",species,"/neiGWAS_asym/data_collection"))){
  dir.create(path = paste0("results/",species,"/neiGWAS_asym/data_collection"))
}

#export data
write_csv(gwas_results, file = paste0("results/",species,"/neiGWAS_asym/",daynumber,"_",poolName,"/neiGWAS_results_asym_",poolName,".csv"))
write_csv(dist_results, file = paste0("results/",species,"/neiGWAS_asym/",daynumber,"_",poolName,"/PVEnei_results_",poolName,".csv"))
write_csv(top_SNPs_self, file = paste0("results/",species,"/neiGWAS_asym/",daynumber,"_",poolName,"/topSNP_self_per_dist_",poolName,".csv"))
write_csv(top_SNPs_nei, file = paste0("results/",species,"/neiGWAS_asym/",daynumber,"_",poolName,"/topSNP_nei_per_dist_",poolName,".csv"))
write_csv(top_SNPs_sxn, file = paste0("results/",species,"/neiGWAS_asym/",daynumber,"_",poolName,"/topSNP_sxn_per_dist_",poolName,".csv"))

#copy neiGWAS_results file manually to "results/[species]/neiGWAS_asym/data_collection" to use them for further analysis

################################# plots of automatic neiGWAS

library(ggplot2)
library(tidyr)

#get all top SNPs to plot
selected_SNPs_self = select(top_SNPs_self, c("chr","pos"))
selected_SNPs_nei = select(top_SNPs_nei, c("chr","pos"))
selected_SNPs_sxn = select(top_SNPs_sxn, c("chr","pos"))
selected_SNPs_self = unique(selected_SNPs_self)
selected_SNPs_nei = unique(selected_SNPs_nei)
selected_SNPs_sxn = unique(selected_SNPs_sxn)

#add p-values for all distances
selected_SNPs_self = merge(selected_SNPs_self,gwas_results)
rows_to_keep_self = c(1,2,(9*(1:((ncol(selected_SNPs_self)-2)/9)))-5)
selected_SNPs_self = select(selected_SNPs_self, rows_to_keep_self)

selected_SNPs_nei = merge(selected_SNPs_nei,gwas_results)
rows_to_keep_nei = c(1,2,(9*(1:((ncol(selected_SNPs_nei)-2)/9)))-2)
selected_SNPs_nei = select(selected_SNPs_nei, rows_to_keep_nei)

selected_SNPs_sxn = merge(selected_SNPs_sxn,gwas_results)
rows_to_keep_sxn = c(1,2,(9*(1:((ncol(selected_SNPs_sxn)-2)/9)))+1)
selected_SNPs_sxn = select(selected_SNPs_sxn, rows_to_keep_sxn)

#PVEnei
if(wished_plots[1]){
  res <- dist_results[order(dist_results$scale), ]
  pve <- res$PVEnei[-1]
  s <- res$scale[-1]
  delta_pve <- pve - c(0, pve[1:(length(pve) - 1)])
  est_scale <- s[delta_pve == max(delta_pve)][1]
  est_pve <- pve[delta_pve == max(delta_pve)][1]
  args <- list()
  args$x <- s
  args$y <- pve
  args$type <- "l"
  args$xlab <- "scale"
  args$ylab <- "PVE_nei"
  png(filename = paste0("results/",species,"/neiGWAS_asym/",daynumber,"_",poolName,"/",poolName,"_PVEnei.png"),width=800,height=600,units="px")
  do.call(graphics::plot, args)
  graphics::points(s, pve)
  graphics::points(est_scale, est_pve, pch = 16)
  dev.off()
  rm(res)
  rm(pve)
  rm(s)
  rm(delta_pve)
  rm(est_scale)
  rm(est_pve)
  rm(args)
}

#cleaning
selected_SNPs_self = mutate(selected_SNPs_self, chrpos = paste0(chr,"_",pos))
selected_SNPs_self = select(selected_SNPs_self, c(ncol(selected_SNPs_self),3:(ncol(selected_SNPs_self)-1)))
selected_SNPs_self = gather(selected_SNPs_self, key = "distance", value = "p_log10", -chrpos)
selected_SNPs_self = mutate(selected_SNPs_self, distance = sapply(strsplit(distance,split="_",fixed=TRUE), function(x) (x[4])))

selected_SNPs_nei = mutate(selected_SNPs_nei, chrpos = paste0(chr,"_",pos))
selected_SNPs_nei = select(selected_SNPs_nei, c(ncol(selected_SNPs_nei),3:(ncol(selected_SNPs_nei)-1)))
selected_SNPs_nei = gather(selected_SNPs_nei, key = "distance", value = "p_log10", -chrpos)
selected_SNPs_nei = mutate(selected_SNPs_nei, distance = sapply(strsplit(distance,split="_",fixed=TRUE), function(x) (x[4])))

selected_SNPs_sxn = mutate(selected_SNPs_sxn, chrpos = paste0(chr,"_",pos))
selected_SNPs_sxn = select(selected_SNPs_sxn, c(ncol(selected_SNPs_sxn),3:(ncol(selected_SNPs_sxn)-1)))
selected_SNPs_sxn = gather(selected_SNPs_sxn, key = "distance", value = "p_log10", -chrpos)
selected_SNPs_sxn = mutate(selected_SNPs_sxn, distance = sapply(strsplit(distance,split="_",fixed=TRUE), function(x) (x[4])))

#discrete plot
x_axis_order = as.factor(unique(as.numeric(selected_SNPs_nei$distance)))
if(wished_plots[4]){
  topSNPs_self_dist_disc = ggplot(selected_SNPs_self, aes(x=distance,y=p_log10,group=chrpos,colour=chrpos)) + geom_line() + geom_hline(yintercept = -log10(0.05/nrow(gwas_results))) + ggtitle("self-p-values of self-top-SNPs of each distance") + scale_x_discrete(limits = x_axis_order) + xlab("distance of included neighbours") + ylab("-log10 p-value")
  ggsave(filename = paste0("results/",species,"/neiGWAS_asym/",daynumber,"_",poolName,"/",poolName,"_topSNPs_self_distance.png"),plot=topSNPs_self_dist_disc,width=4000,height=3000,units="px",device="png")
  rm(topSNPs_self_dist_disc)
  topSNPs_nei_dist_disc = ggplot(selected_SNPs_nei, aes(x=distance,y=p_log10,group=chrpos,colour=chrpos)) + geom_line() + geom_hline(yintercept = -log10(0.05/nrow(gwas_results))) + ggtitle("nei-p-values of nei-top-SNPs of each distance") + scale_x_discrete(limits = x_axis_order) + xlab("distance of included neighbours") + ylab("-log10 p-value")
  ggsave(filename = paste0("results/",species,"/neiGWAS_asym/",daynumber,"_",poolName,"/",poolName,"_topSNPs_nei_distance.png"),plot=topSNPs_nei_dist_disc,width=4000,height=3000,units="px",device="png")
  rm(topSNPs_nei_dist_disc)
  topSNPs_sxn_dist_disc = ggplot(selected_SNPs_sxn, aes(x=distance,y=p_log10,group=chrpos,colour=chrpos)) + geom_line() + geom_hline(yintercept = -log10(0.05/nrow(gwas_results))) + ggtitle("sxn-p-values of sxn-top-SNPs of each distance") + scale_x_discrete(limits = x_axis_order) + xlab("distance of included neighbours") + ylab("-log10 p-value")
  ggsave(filename = paste0("results/",species,"/neiGWAS_asym/",daynumber,"_",poolName,"/",poolName,"_topSNPs_sxn_distance.png"),plot=topSNPs_sxn_dist_disc,width=4000,height=3000,units="px",device="png")
  rm(topSNPs_sxn_dist_disc)
}
#discrete plot small
selected_SNPs_self_small = filter(selected_SNPs_self, as.numeric(distance) < 6)
selected_SNPs_nei_small = filter(selected_SNPs_nei, as.numeric(distance) < 6)
selected_SNPs_sxn_small = filter(selected_SNPs_sxn, as.numeric(distance) < 6)
x_axis_order = factor(x_axis_order[1:18],exclude=NULL)
if(wished_plots[5]){
  topSNPs_self_dist_disc_small = ggplot(selected_SNPs_self_small, aes(x=distance,y=p_log10,group=chrpos,colour=chrpos)) + geom_line() + geom_hline(yintercept = -log10(0.05/nrow(gwas_results))) + ggtitle("self-p-values of self-top-SNPs of each distance") + scale_x_discrete(limits = x_axis_order) + xlab("distance of included neighbours") + ylab("-log10 p-value")
  ggsave(filename = paste0("results/",species,"/neiGWAS_asym/",daynumber,"_",poolName,"/",poolName,"_topSNPs_self_distance_small.png"),plot=topSNPs_self_dist_disc_small,width=4000,height=3000,units="px",device="png")
  rm(topSNPs_self_dist_disc_small)
  topSNPs_nei_dist_disc_small = ggplot(selected_SNPs_nei_small, aes(x=distance,y=p_log10,group=chrpos,colour=chrpos)) + geom_line() + geom_hline(yintercept = -log10(0.05/nrow(gwas_results))) + ggtitle("nei-p-values of nei-top-SNPs of each distance") + scale_x_discrete(limits = x_axis_order) + xlab("distance of included neighbours") + ylab("-log10 p-value")
  ggsave(filename = paste0("results/",species,"/neiGWAS_asym/",daynumber,"_",poolName,"/",poolName,"_topSNPs_nei_distance_small.png"),plot=topSNPs_nei_dist_disc_small,width=4000,height=3000,units="px",device="png")
  rm(topSNPs_nei_dist_disc_small)
  topSNPs_sxn_dist_disc_small = ggplot(selected_SNPs_sxn_small, aes(x=distance,y=p_log10,group=chrpos,colour=chrpos)) + geom_line() + geom_hline(yintercept = -log10(0.05/nrow(gwas_results))) + ggtitle("sxn-p-values of sxn-top-SNPs of each distance") + scale_x_discrete(limits = x_axis_order) + xlab("distance of included neighbours") + ylab("-log10 p-value")
  ggsave(filename = paste0("results/",species,"/neiGWAS_asym/",daynumber,"_",poolName,"/",poolName,"_topSNPs_sxn_distance_small.png"),plot=topSNPs_sxn_dist_disc_small,width=4000,height=3000,units="px",device="png")
  rm(topSNPs_sxn_dist_disc_small)
}

#continuous plot
selected_SNPs_self = mutate(selected_SNPs_self, distance = as.numeric(distance))
selected_SNPs_nei = mutate(selected_SNPs_nei, distance = as.numeric(distance))
selected_SNPs_sxn = mutate(selected_SNPs_sxn, distance = as.numeric(distance))
if(wished_plots[2]){
  topSNPs_self_dist_cont = ggplot(selected_SNPs_self, aes(x=distance,y=p_log10,group=chrpos,colour=chrpos)) + geom_line() + geom_hline(yintercept = -log10(0.05/nrow(gwas_results))) + ggtitle("self-p-values of self-top-SNPs of each distance") + xlab("distance of included neighbours") + ylab("-log10 p-value")
  ggsave(filename = paste0("results/",species,"/neiGWAS_asym/",daynumber,"_",poolName,"/",poolName,"_topSNPs_self_distance_continuous.png"),plot=topSNPs_self_dist_cont,width=4000,height=3000,units="px",device="png")
  rm(topSNPs_self_dist_cont)
  topSNPs_nei_dist_cont = ggplot(selected_SNPs_nei, aes(x=distance,y=p_log10,group=chrpos,colour=chrpos)) + geom_line() + geom_hline(yintercept = -log10(0.05/nrow(gwas_results))) + ggtitle("nei-p-values of nei-top-SNPs of each distance") + xlab("distance of included neighbours") + ylab("-log10 p-value")
  ggsave(filename = paste0("results/",species,"/neiGWAS_asym/",daynumber,"_",poolName,"/",poolName,"_topSNPs_nei_distance_continuous.png"),plot=topSNPs_nei_dist_cont,width=4000,height=3000,units="px",device="png")
  rm(topSNPs_nei_dist_cont)
  topSNPs_sxn_dist_cont = ggplot(selected_SNPs_sxn, aes(x=distance,y=p_log10,group=chrpos,colour=chrpos)) + geom_line() + geom_hline(yintercept = -log10(0.05/nrow(gwas_results))) + ggtitle("sxn-p-values of sxn-top-SNPs of each distance") + xlab("distance of included neighbours") + ylab("-log10 p-value")
  ggsave(filename = paste0("results/",species,"/neiGWAS_asym/",daynumber,"_",poolName,"/",poolName,"_topSNPs_sxn_distance_continuous.png"),plot=topSNPs_sxn_dist_cont,width=4000,height=3000,units="px",device="png")
  rm(topSNPs_sxn_dist_cont)
}
#continuous plot small
selected_SNPs_self_small = mutate(selected_SNPs_self_small, distance = as.numeric(distance))
selected_SNPs_nei_small = mutate(selected_SNPs_nei_small, distance = as.numeric(distance))
selected_SNPs_sxn_small = mutate(selected_SNPs_sxn_small, distance = as.numeric(distance))
if(wished_plots[3]){
  topSNPs_self_dist_cont_small = ggplot(selected_SNPs_self_small, aes(x=distance,y=p_log10,group=chrpos,colour=chrpos)) + geom_line() + geom_hline(yintercept = -log10(0.05/nrow(gwas_results))) + ggtitle("self-p-values of self-top-SNPs of each distance") + xlab("distance of included neighbours") + ylab("-log10 p-value")
  ggsave(filename = paste0("results/",species,"/neiGWAS_asym/",daynumber,"_",poolName,"/",poolName,"_topSNPs_self_distance_continuous_small.png"),plot=topSNPs_self_dist_cont_small,width=4000,height=3000,units="px",device="png")
  rm(topSNPs_self_dist_cont_small)
  topSNPs_nei_dist_cont_small = ggplot(selected_SNPs_nei_small, aes(x=distance,y=p_log10,group=chrpos,colour=chrpos)) + geom_line() + geom_hline(yintercept = -log10(0.05/nrow(gwas_results))) + ggtitle("nei-p-values of nei-top-SNPs of each distance") + xlab("distance of included neighbours") + ylab("-log10 p-value")
  ggsave(filename = paste0("results/",species,"/neiGWAS_asym/",daynumber,"_",poolName,"/",poolName,"_topSNPs_nei_distance_continuous_small.png"),plot=topSNPs_nei_dist_cont_small,width=4000,height=3000,units="px",device="png")
  rm(topSNPs_nei_dist_cont_small)
  topSNPs_sxn_dist_cont_small = ggplot(selected_SNPs_sxn_small, aes(x=distance,y=p_log10,group=chrpos,colour=chrpos)) + geom_line() + geom_hline(yintercept = -log10(0.05/nrow(gwas_results))) + ggtitle("sxn-p-values of sxn-top-SNPs of each distance") + xlab("distance of included neighbours") + ylab("-log10 p-value")
  ggsave(filename = paste0("results/",species,"/neiGWAS_asym/",daynumber,"_",poolName,"/",poolName,"_topSNPs_sxn_distance_continuous_small.png"),plot=topSNPs_sxn_dist_cont_small,width=4000,height=3000,units="px",device="png")
  rm(topSNPs_sxn_dist_cont_small)
}

#manhattan plots for all distances
if(wished_plots[6]){
  if(!file.exists(paste0("results/",species,"/neiGWAS_asym/",daynumber,"_",poolName,"/manhattan"))){
    dir.create(path = paste0("results/",species,"/neiGWAS_asym/",daynumber,"_",poolName,"/manhattan"))
  }
  result_rows_self = (9*(1:((ncol(gwas_results)-2)/9)))-6
  result_rows_nei = result_rows_self + 3
  result_rows_sxn = result_rows_nei + 3
  for(j in result_rows_self){
    curr_data = select(gwas_results, c(1,2,j))
    d = unlist(strsplit(colnames(curr_data)[3],split="_",fixed=TRUE))[3]
    colnames(curr_data)<- c("chr","pos","p")
    png(filename = paste0("results/",species,"/neiGWAS_asym/",daynumber,"_",poolName,"/manhattan/",poolName,"_manhattan_self_d",d,".png"),width=1600,height=800,units="px")
    gaston::manhattan(curr_data)
    abline(h = -log10(0.05/ncol(geno))) #bonferroni-corrected significance line, may not shows if no point close or above
    dev.off()
    png(filename = paste0("results/",species,"/neiGWAS_asym/",daynumber,"_",poolName,"/manhattan/",poolName,"_qq_self_d",d,".png"),width=1600,height=800,units="px")
    gaston::qqplot.pvalues(curr_data[,3])
    dev.off()
  }
  for(j in result_rows_nei){
    curr_data = select(gwas_results, c(1,2,j))
    d = unlist(strsplit(colnames(curr_data)[3],split="_",fixed=TRUE))[3]
    colnames(curr_data)<- c("chr","pos","p")
    png(filename = paste0("results/",species,"/neiGWAS_asym/",daynumber,"_",poolName,"/manhattan/",poolName,"_manhattan_nei_d",d,".png"),width=1600,height=800,units="px")
    gaston::manhattan(curr_data)
    abline(h = -log10(0.05/ncol(geno))) #bonferroni-corrected significance line, may not shows if no point close or above
    dev.off()
    png(filename = paste0("results/",species,"/neiGWAS_asym/",daynumber,"_",poolName,"/manhattan/",poolName,"_qq_nei_d",d,".png"),width=1600,height=800,units="px")
    gaston::qqplot.pvalues(curr_data[,3])
    dev.off()
  }
  for(j in result_rows_sxn){
    curr_data = select(gwas_results, c(1,2,j))
    d = unlist(strsplit(colnames(curr_data)[3],split="_",fixed=TRUE))[3]
    colnames(curr_data)<- c("chr","pos","p")
    png(filename = paste0("results/",species,"/neiGWAS_asym/",daynumber,"_",poolName,"/manhattan/",poolName,"_manhattan_sxn_d",d,".png"),width=1600,height=800,units="px")
    gaston::manhattan(curr_data)
    abline(h = -log10(0.05/ncol(geno))) #bonferroni-corrected significance line, may not shows if no point close or above
    dev.off()
    png(filename = paste0("results/",species,"/neiGWAS_asym/",daynumber,"_",poolName,"/manhattan/",poolName,"_qq_sxn_d",d,".png"),width=1600,height=800,units="px")
    gaston::qqplot.pvalues(curr_data[,3])
    dev.off()
  }
  rm(result_rows_self)
  rm(result_rows_nei)
  rm(result_rows_sxn)
  rm(curr_data)
  rm(d)
  rm(j)
}

#cleaning up
rm(selected_SNPs_self)
rm(selected_SNPs_nei)
rm(selected_SNPs_sxn)
rm(selected_SNPs_self_small)
rm(selected_SNPs_nei_small)
rm(selected_SNPs_sxn_small)
rm(rows_to_keep_self)
rm(rows_to_keep_nei)
rm(rows_to_keep_sxn)
rm(x_axis_order)
detach(package:ggplot2)
detach(package:tidyr)

