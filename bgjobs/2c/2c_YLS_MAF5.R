library(rNeighborGWAS)
library(gaston)
library(dplyr)
library(readr)

#this script is to do neiGWAS over all distances for a given pheno

#note: from this script on, barley and wheat share scripts

################################# settings
#you need to adjust this before running

#day number for easier identification of outputs
daynumber = "20240425"
#do not use slashes or points, since this will be used for the output folder name

#pool settings:
species = "wheat"
#barley or wheat
poolName = "YLS_2017-MAF5"
#barley pools:
#NFNB_pooled, PM_2015, Blr_pooled, Scald_pooled, SFNB_pooled
#wheat pools:
#YR_pooled, LR_2017, SR_2017, YLS_2017, STB_2019, YS_2019

#cores to use for calc_PVEnei and neiGWAS
cores = 4L

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
top_SNPs = data.frame(NULL)

#distances up to 6:
#distances = c(1,1.5,2,2.3,2.9,3,3.2,3.7,4,4.2,4.3,4.5,5,5.2,5.4,5.7,5.9,6)
#distances up to 10:
distances = c(1,1.5,2,2.3,2.9,3,3.2,3.7,4,4.2,4.3,4.5,5,5.2,5.4,5.7,5.9,6,6.1,6.4,6.5,6.8,7,7.1,7.25,7.3,7.7,7.9,8,8.1,8.3,8.5,8.6,8.7,8.98,9,9.1,9.3,9.45,9.5,9.9,9.95,10)

dist_results = calc_PVEnei(geno=geno,pheno=pheno$Damage_Level,smap=smap,scale=distances,addcovar=as.factor(pheno$Experiment_Number),grouping=pheno$Experiment_Number,n_core=cores)

for(scale in distances){
  gwas_out <- neiGWAS(geno=geno,pheno=pheno$Damage_Level,
                      gmap=gmap, smap=smap,
                      scale=scale, n_core=cores,
                      addcovar=as.factor(pheno$Experiment_Number),
                      grouping=pheno$Experiment_Number
  )
  gwas_out = mutate(gwas_out, p_log10 = -log10(p))
  colnames(gwas_out)<-c("chr","pos",paste0("p_",toString(scale)),paste0("p_log10_",toString(scale)))
  if(ncol(gwas_results)==0){
    gwas_results = gwas_out
  }
  else{
    gwas_results = merge(gwas_results,gwas_out)
  }
  gwas_out = arrange(gwas_out, gwas_out[3])
  topsnp = gwas_out[1,]
  topsnp = append(topsnp, scale)
  if(ncol(top_SNPs)==0){
    top_SNPs = data.frame(topsnp)
    colnames(top_SNPs)<-c("chr","pos","p","p_log10","dist")
  }
  else{
    topsnp = data.frame(topsnp)
    colnames(topsnp)<-c("chr","pos","p","p_log10","dist")
    top_SNPs = add_row(top_SNPs, topsnp)
  }
}

#create output folder if necessary
if(!file.exists(paste0("results/",species,"/neiGWAS"))){
  dir.create(path = paste0("results/",species,"/neiGWAS"))
}
if(!file.exists(paste0("results/",species,"/neiGWAS/",daynumber,"_",poolName))){
  dir.create(path = paste0("results/",species,"/neiGWAS/",daynumber,"_",poolName))
}
if(!file.exists(paste0("results/",species,"/neiGWAS/data_collection"))){
  dir.create(path = paste0("results/",species,"/neiGWAS/data_collection"))
}

#export data
write_csv(gwas_results, file = paste0("results/",species,"/neiGWAS/",daynumber,"_",poolName,"/neiGWAS_results_",poolName,".csv"))
write_csv(dist_results, file = paste0("results/",species,"/neiGWAS/",daynumber,"_",poolName,"/PVEnei_results_",poolName,".csv"))
write_csv(top_SNPs, file = paste0("results/",species,"/neiGWAS/",daynumber,"_",poolName,"/topSNP_per_dist_",poolName,".csv"))

#copy neiGWAS_results file manually to "results/[species]/neiGWAS/data_collection" to use them for further analysis

################################# plots of automatic neiGWAS

library(ggplot2)
library(tidyr)

#get all top SNPs to plot
selected_SNPs = select(top_SNPs, c("chr","pos"))
selected_SNPs = unique(selected_SNPs)

#add p-values for all distances
selected_SNPs = merge(selected_SNPs,gwas_results)
rows_to_keep = c(1,2*(1:(ncol(selected_SNPs)/2)))
selected_SNPs = select(selected_SNPs, rows_to_keep)

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
  png(filename = paste0("results/",species,"/neiGWAS/",daynumber,"_",poolName,"/",poolName,"_PVEnei.png"),width=800,height=600,units="px")
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
selected_SNPs = mutate(selected_SNPs, chrpos = paste0(chr,"_",pos))
selected_SNPs = select(selected_SNPs, c(ncol(selected_SNPs),3:(ncol(selected_SNPs)-1)))
selected_SNPs = gather(selected_SNPs, key = "distance", value = "p_log10", -chrpos)
selected_SNPs = mutate(selected_SNPs, distance = sapply(strsplit(distance,split="_",fixed=TRUE), function(x) (x[3])))

#discrete plot
x_axis_order = as.factor(unique(as.numeric(selected_SNPs$distance)))
if(wished_plots[4]){
  topSNPs_dist_disc = ggplot(selected_SNPs, aes(x=distance,y=p_log10,group=chrpos,colour=chrpos)) + geom_line() + geom_hline(yintercept = -log10(0.05/nrow(gwas_results))) + ggtitle("p-values of top-SNPs of each distance") + scale_x_discrete(limits = x_axis_order) + xlab("distance of included neighbours") + ylab("-log10 p-value")
  ggsave(filename = paste0("results/",species,"/neiGWAS/",daynumber,"_",poolName,"/",poolName,"_topSNPs_distance.png"),plot=topSNPs_dist_disc,width=4000,height=3000,units="px",device="png")
  rm(topSNPs_dist_disc)
}
#discrete plot small
selected_SNPs_small = filter(selected_SNPs, as.numeric(distance) < 6)
x_axis_order = factor(x_axis_order[1:18],exclude=NULL)
if(wished_plots[5]){
  topSNPs_dist_disc_small = ggplot(selected_SNPs_small, aes(x=distance,y=p_log10,group=chrpos,colour=chrpos)) + geom_line() + geom_hline(yintercept = -log10(0.05/nrow(gwas_results))) + ggtitle("p-values of top-SNPs of each distance") + scale_x_discrete(limits = x_axis_order) + xlab("distance of included neighbours") + ylab("-log10 p-value")
  ggsave(filename = paste0("results/",species,"/neiGWAS/",daynumber,"_",poolName,"/",poolName,"_topSNPs_distance_small.png"),plot=topSNPs_dist_disc_small,width=4000,height=3000,units="px",device="png")
  rm(topSNPs_dist_disc_small)
}

#continuous plot
selected_SNPs = mutate(selected_SNPs, distance = as.numeric(distance))
if(wished_plots[2]){
  topSNPs_dist_cont = ggplot(selected_SNPs, aes(x=distance,y=p_log10,group=chrpos,colour=chrpos)) + geom_line() + geom_hline(yintercept = -log10(0.05/nrow(gwas_results))) + ggtitle("p-values of top-SNPs of each distance") + xlab("distance of included neighbours") + ylab("-log10 p-value")
  ggsave(filename = paste0("results/",species,"/neiGWAS/",daynumber,"_",poolName,"/",poolName,"_topSNPs_distance_continuous.png"),plot=topSNPs_dist_cont,width=4000,height=3000,units="px",device="png")
  rm(topSNPs_dist_cont)
}
#continuous plot small
selected_SNPs_small = mutate(selected_SNPs_small, distance = as.numeric(distance))
if(wished_plots[3]){
  topSNPs_dist_cont_small = ggplot(selected_SNPs_small, aes(x=distance,y=p_log10,group=chrpos,colour=chrpos)) + geom_line() + geom_hline(yintercept = -log10(0.05/nrow(gwas_results))) + ggtitle("p-values of top-SNPs of each distance") + xlab("distance of included neighbours") + ylab("-log10 p-value")
  ggsave(filename = paste0("results/",species,"/neiGWAS/",daynumber,"_",poolName,"/",poolName,"_topSNPs_distance_continuous_small.png"),plot=topSNPs_dist_cont_small,width=4000,height=3000,units="px",device="png")
  rm(topSNPs_dist_cont_small)
}

#manhattan plots for all distances
if(wished_plots[6]){
  if(!file.exists(paste0("results/",species,"/neiGWAS/",daynumber,"_",poolName,"/manhattan"))){
    dir.create(path = paste0("results/",species,"/neiGWAS/",daynumber,"_",poolName,"/manhattan"))
  }
  result_rows = (2*(2:(ncol(gwas_results)/2)))-1
  for(j in result_rows){
    curr_data = select(gwas_results, c(1,2,j))
    d = unlist(strsplit(colnames(curr_data)[3],split="_",fixed=TRUE))[2]
    colnames(curr_data)<- c("chr","pos","p")
    png(filename = paste0("results/",species,"/neiGWAS/",daynumber,"_",poolName,"/manhattan/",poolName,"_manhattan_d",d,".png"),width=1600,height=800,units="px")
    gaston::manhattan(curr_data)
    abline(h = -log10(0.05/ncol(geno))) #bonferroni-corrected significance line, may not shows if no point close or above
    dev.off()
    png(filename = paste0("results/",species,"/neiGWAS/",daynumber,"_",poolName,"/manhattan/",poolName,"_qq_d",d,".png"),width=1600,height=800,units="px")
    gaston::qqplot.pvalues(curr_data[,3])
    dev.off()
  }
  rm(result_rows)
  rm(curr_data)
  rm(d)
  rm(j)
}

#cleaning up
rm(selected_SNPs)
rm(selected_SNPs_small)
rm(rows_to_keep)
rm(x_axis_order)
detach(package:ggplot2)
detach(package:tidyr)

