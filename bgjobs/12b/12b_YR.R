library(rNeighborGWAS)
library(gaston)
library(dplyr)
library(readr)

#standalone script for PVE analysis with multiple MAFs,
#with extended range of nei dist up to 16 and more MAFs

#can be run directly after scripts "1"

################################# settings
#you need to adjust this before running

#day number for easier identification of outputs
daynumber = "20240513"
#do not use slashes or points, since this will be used for the output folder name

#pool settings:
species = "wheat"
#barley or wheat
poolName = "YR_pooled"
#barley pools:
#NFNB_pooled, PM_2015, Blr_pooled, Scald_pooled, SFNB_pooled
#wheat pools:
#YR_pooled, LR_2017, SR_2017, YLS_2017, STB_2019, YS_2019
pool = c(4,6)
#pools:
#NFNB: 1,6,10
#PM_2015: 2
#Blr: 5,9
#Scald: 3,7,11
#SFNB: 4,8,12
#YR_pooled: 4,6
#LR_2017: 1
#SR_2017: 2
#YLS_2017: 3
#STB_2019: 5
#YS_2019: 7

#MAFs to use
MAF = c(0.0,0.005,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.12,0.15,0.2,0.25)

#cores to use for calc_PVEnei
cores = 4L

################################# do stuff

#create output folder if necessary
if(!file.exists(paste0("results"))){
  dir.create(paste0(path = "results"))
}
if(!file.exists(paste0("results/",species))){
  dir.create(paste0(path = "results/",species))
}
if(!file.exists(paste0("results/",species,"/PVE_analysis"))){
  dir.create(paste0(path = "results/",species,"/PVE_analysis"))
}
if(!file.exists(paste0("results/",species,"/PVE_analysis/multiMAF"))){
  dir.create(paste0(path = "results/",species,"/PVE_analysis/multiMAF"))
}
if(!file.exists(paste0("results/",species,"/PVE_analysis/multiMAF/extended"))){
  dir.create(paste0(path = "results/",species,"/PVE_analysis/multiMAF/extended"))
}
if(!file.exists(paste0("results/",species,"/PVE_analysis/multiMAF/extended/",daynumber,"_",poolName))){
  dir.create(paste0(path = "results/",species,"/PVE_analysis/multiMAF/extended/",daynumber,"_",poolName))
}

PVE_results_multiMAF = data.frame(NULL)

for(MAF_i in MAF){
  ################################# loading and preparing data
  #before running check manually in folders if all required files are present
  if(species=="barley"){
    fn = system("ls ./pheno/barley/*.csv", intern=TRUE)
    pheno_all = data.frame(NULL)
    for(i in 1:12){
      pheno_i = read.csv(fn[i])
      pheno_i[,ncol(pheno_i)] = as.numeric(pheno_i[,ncol(pheno_i)])
      pheno_i = pheno_i[is.na(pheno_i[,ncol(pheno_i)])==FALSE,]
      pheno_i = rename(pheno_i, Damage_Level = ncol(pheno_i))
      pheno_i = rename(pheno_i, Expt = 6)
      if(ncol(pheno_i)>8){
        pheno_i = select(pheno_i, c(1,2,3,4,5,6,7,ncol(pheno_i)))
      }
      pheno_i = mutate(pheno_i, Experiment_Number = i)
      pheno_i = mutate(pheno_i, Damage_Level = scale(Damage_Level)[,1])
      if(i==1){
        pheno_all = pheno_i
      }
      else{
        pheno_all = rows_append(pheno_all, pheno_i)
      }
    }
    #select phenotype pool
    pheno = filter(pheno_all, Experiment_Number %in% pool)
    # load an annotation file linking genotypes and phenotypes
    gf = read.csv("./pheno/barley/geno2pheno.csv")
    gf = gf[is.na(gf$GID)==FALSE,]
    pheno$Name %in% gf$EntryCD
    pheno = pheno[(pheno$Name %in% gf$EntryCD),]
    gf = gf[match(pheno$Name,gf$EntryCD),]
    gf$GID =  paste0("X",gf$GID)
    geno = read.csv("./geno/barley/geno.csv")
    pheno = pheno[(gf$GID %in% colnames(geno)),]
    andIDs = gf$GID[(gf$GID %in% colnames(geno))]
    geno = geno[,andIDs]
    #prepare geno, filter minor alleles
    af = rowSums(geno)/(2*ncol(geno)) #assuming geno is complete
    rare_allele = ifelse(af<MAF_i|1-af<MAF_i,TRUE,FALSE)
    geno = filter(geno,!rare_allele)
    af = af[!rare_allele]
    geno = t(geno)
    geno = geno - 1
    #smap
    smap = cbind(pheno$Range,pheno$Row)
    #gmap
    gmap = read.csv("./geno/barley/positions.csv")
    gmap = filter(gmap,!rare_allele)
  }
  else if(species=="wheat"){
    fn = system("ls ./pheno/wheat/*.csv", intern=TRUE)
    pheno_all = data.frame(NULL)
    for(i in 1:7){
      pheno_i = read.csv(fn[i])
      pheno_i[,ncol(pheno_i)] = as.numeric(pheno_i[,ncol(pheno_i)])
      pheno_i = pheno_i[is.na(pheno_i[,ncol(pheno_i)])==FALSE,]
      pheno_i = rename(pheno_i, Damage_Level = ncol(pheno_i))
      if(i==5){
        pheno_i = mutate(pheno_i, Block = 1)
        pheno_i = select(pheno_i, c(1,2,11,3,4,5,6,7,8,9,10))
      }
      pheno_i = rename(pheno_i, Expt = 6)
      pheno_i = rename(pheno_i, Name = 7)
      if(ncol(pheno_i)>8){
        pheno_i = select(pheno_i, c(1,2,3,4,5,6,7,ncol(pheno_i)))
      }
      pheno_i = mutate(pheno_i, Experiment_Number = i)
      pheno_i = mutate(pheno_i, Damage_Level = scale(Damage_Level)[,1])
      if(i==1){
        pheno_all = pheno_i
      }
      else{
        pheno_all = rows_append(pheno_all, pheno_i)
      }
    }
    #select phenotype pool
    pheno = filter(pheno_all, Experiment_Number %in% pool)
    # load an annotation file linking genotypes and phenotypes
    gf = read.csv("./pheno/wheat/geno2pheno_wheat.csv")
    gf = gf[is.na(gf$GID)==FALSE,]
    pheno$Name %in% gf$EntryCD
    pheno = pheno[(pheno$Name %in% gf$EntryCD),]
    gf = gf[match(pheno$Name,gf$EntryCD),]
    gf$GID =  paste0("GID",gf$GID)
    geno = read.csv("./geno/wheat/geno.csv")
    pheno = pheno[(gf$GID %in% colnames(geno)),]
    andIDs = gf$GID[(gf$GID %in% colnames(geno))]
    geno = geno[,andIDs]
    #prepare geno, filter minor alleles
    af = rowSums(geno)/(2*ncol(geno)) #assuming geno is complete
    rare_allele = ifelse(af<MAF_i|1-af<MAF_i,TRUE,FALSE)
    geno = filter(geno,!rare_allele)
    af = af[!rare_allele]
    geno = t(geno)
    geno = geno - 1
    #smap
    smap = cbind(pheno$Range,pheno$Row)
    #gmap
    gmap = read.csv("./geno/wheat/positions.csv")
    gmap = filter(gmap,!rare_allele)
  }
  ################################# PVE
  #distances up to 5: (13)
  #distances = c(1,1.5,2,2.3,2.9,3,3.2,3.7,4,4.2,4.3,4.5,5)
  #distances up to 6: (18)
  #distances = c(1,1.5,2,2.3,2.9,3,3.2,3.7,4,4.2,4.3,4.5,5,5.2,5.4,5.7,5.9,6)
  #distances up to 10: (43)
  #distances = c(1,1.5,2,2.3,2.9,3,3.2,3.7,4,4.2,4.3,4.5,5,5.2,5.4,5.7,5.9,6,6.1,6.4,6.5,6.8,7,7.1,7.25,7.3,7.7,7.9,8,8.1,8.3,8.5,8.6,8.7,8.98,9,9.1,9.3,9.45,9.5,9.9,9.95,10)
  #distances up to 16: (98)
  distances = c(1,1.5,2,2.3,2.9,3,3.2,3.7,4,4.2,4.3,4.5,5,5.2,5.4,5.7,5.9,6,6.1,6.4,6.5,6.8,7,7.1,7.25,7.3,7.7,7.9,8,8.1,8.3,8.5,8.6,8.7,8.98,9,9.1,9.3,9.45,9.5,9.9,9.95,10,10.1,10.2,10.3,10.5,10.7,10.8,10.9,11,11.1,11.2,11.4,11.5,11.7,11.8,12,12.05,12.1,12.2,12.3,12.4,12.6,12.7,12.8,12.9,13,13.1,13.2,13.4,13.45,13.5,13.7,13.9,13.95,14,14.1,14.2,14.3,14.4,14.5,14.6,14.8,14.9,15,15.1,15.2,15.25,15.275,15.3,15.55,15.6,15.65,15.7,15.9,16)
  PVE_results = data.frame(NULL)
  for(scale in distances){
    nSNPs = ncol(geno)
    nPlants = nrow(geno)
    #K_self
    K_self = tcrossprod(geno)
    K_self = ((nSNPs - 1)/2 + K_self/2)/(nSNPs - 1)
    K_self = as.matrix(Matrix::nearPD(K_self,maxit=10^6)$mat)
    #K_nei
    g_nei = nei_coval(geno=geno, smap=smap, scale=scale,
                      grouping = pheno$Experiment_Number,
                      n_core = cores)
    K_nei = tcrossprod(g_nei)/(nSNPs - 1)
    K_nei = as.matrix(Matrix::nearPD(K_nei,maxit=10^6)$mat)
    #D: proximity matrix
    D = as.matrix(dist(smap))
    D = exp(-D) # define exponential distance decay
    D = as.matrix(Matrix::nearPD(D, maxit=10^6)$mat)
    #overwriting if positions are not in same experiment, so spatially unrelated
    for(j in 1:ncol(D)){
      for(k in j:nrow(D)){
        if(pheno$Experiment_Number[j]!=pheno$Experiment_Number[k]){
          D[j,k] = 0.0
          D[k,j] = 0.0
        }
      }
    }
    if(length(unique(pheno$Experiment_Number))>1){
      #X: experiment number as fixed effect covariate
      X = model.matrix(~factor(Experiment_Number),data=pheno)
      #perform variation partitioning
      #PVE only by self genotypes
      res0 = lmm.aireml(Y=pheno$Damage_Level, X=X, K=list(K_self), verbose=FALSE)
      PVEs = res0$tau[1] / sum(res0$tau,res0$sigma2)
      p_s = pchisq(2*(res0$logL - res0$logL0),1,lower.tail=FALSE)
      # PVE by self + neighbor genotypes
      res1 = gaston::lmm.aireml(Y=pheno$Damage_Level, X=X, K=list(K_self,K_nei), verbose=FALSE)
      PVEsn = sum(res1$tau) / sum(res1$sigma2,res1$tau)
      p_sn = pchisq(2*(res1$logL - res0$logL),1,lower.tail=FALSE)
      #PVE by self + neighbor genotypes + spatial distance
      res2 = gaston::lmm.aireml(Y=pheno$Damage_Level, X=X, K=list(K_self,K_nei,D), verbose=FALSE)
      PVEsd = sum(res2$tau) / sum(res2$sigma2,res2$tau)
      p_sd = pchisq(2*(res2$logL - res0$logL),1,lower.tail=FALSE)
    }
    else{
      #perform variation partitioning
      #PVE only by self genotypes
      res0 = lmm.aireml(Y=pheno$Damage_Level, K=list(K_self), verbose=FALSE)
      PVEs = res0$tau[1] / sum(res0$tau,res0$sigma2)
      p_s = pchisq(2*(res0$logL - res0$logL0),1,lower.tail=FALSE)
      # PVE by self + neighbor genotypes
      res1 = gaston::lmm.aireml(Y=pheno$Damage_Level, K=list(K_self,K_nei), verbose=FALSE)
      PVEsn = sum(res1$tau) / sum(res1$sigma2,res1$tau)
      p_sn = pchisq(2*(res1$logL - res0$logL),1,lower.tail=FALSE)
      #PVE by self + neighbor genotypes + spatial distance
      res2 = gaston::lmm.aireml(Y=pheno$Damage_Level, K=list(K_self,K_nei,D), verbose=FALSE)
      PVEsd = sum(res2$tau) / sum(res2$sigma2,res2$tau)
      p_sd = pchisq(2*(res2$logL - res0$logL),1,lower.tail=FALSE)
    }
    #pack results
    PVE_results_i = data.frame("distance"=scale,"MAF"=MAF_i,"SNPs"=nSNPs,"accessions"=nPlants,"PVEs"=PVEs,"p_s"=p_s,"PVEsn"=PVEsn,"p_sn"=p_sn,"PVEsnd"=PVEsd,"p_snd"=p_sd)
    if(nrow(PVE_results)==0){
      PVE_results = PVE_results_i
    }
    else{
      PVE_results = rows_append(PVE_results, PVE_results_i)
    }
  }
  ################################# save results
  if(nrow(PVE_results_multiMAF)==0){
    PVE_results_multiMAF = PVE_results
  }
  else{
    PVE_results_multiMAF = rows_append(PVE_results_multiMAF, PVE_results)
  }
  #create output folder if necessary
  if(!file.exists(paste0("results/",species,"/PVE_analysis/multiMAF/extended/",daynumber,"_",poolName,"/MAF",MAF_i))){
    dir.create(paste0(path = "results/",species,"/PVE_analysis/multiMAF/extended/",daynumber,"_",poolName,"/MAF",MAF_i))
  }
  write_csv(PVE_results, file = paste0("results/",species,"/PVE_analysis/multiMAF/extended/",daynumber,"_",poolName,"/MAF",MAF_i,"/PVE_results_",poolName,".csv"))
  ################################# plot
  library(ggplot2)
  library(tidyr)
  PVE_results = mutate(PVE_results, PVEn = PVEsn - PVEs)
  PVE_results = mutate(PVE_results, PVEd = PVEsnd - PVEsn)
  plot_data = PVE_results[,c(1,5,11,12)]
  colnames(plot_data)<- c("dist","self","nei","env")
  plot_data = gather(plot_data, key = "type", value = "PVE", -dist)
  PVE_plot = ggplot(plot_data, aes(x = dist, y = PVE, fill = type)) + geom_bar(position="stack",stat="identity") + ggtitle("PVE partitioning of each distance") + xlab("distance of included neighbours") + ylab("PVE")
  ggsave(filename = paste0("results/",species,"/PVE_analysis/multiMAF/extended/",daynumber,"_",poolName,"/MAF",MAF_i,"/",poolName,"_PVE_vs_distance.png"),plot=PVE_plot,width=4000,height=3000,units="px",device="png")
  rm(PVE_plot)
  #clean
  rm(D)
  rm(g_nei)
  rm(geno)
  rm(gmap)
  rm(K_nei)
  rm(K_self)
  rm(pheno)
  rm(plot_data)
  rm(PVE_results)
  rm(PVE_results_i)
  rm(res0)
  rm(res1)
  rm(res2)
  rm(smap)
  rm(X)
  rm(af)
  rm(j)
  rm(k)
  rm(nPlants)
  rm(nSNPs)
  rm(p_s)
  rm(p_sd)
  rm(p_sn)
  rm(PVEs)
  rm(PVEsd)
  rm(PVEsn)
  rm(scale)
  rm(gf)
  rm(pheno_all)
  rm(pheno_i)
  rm(andIDs)
  rm(fn)
  rm(i)
  rm(rare_allele)
  detach(package:ggplot2)
  detach(package:tidyr)
}

#cleaning
rm(cores)
rm(MAF_i)
rm(pool)

################################# save results

write_csv(PVE_results_multiMAF, file = paste0("results/",species,"/PVE_analysis/multiMAF/extended/",daynumber,"_",poolName,"/PVE_results_multiMAF_",poolName,".csv"))

################################# plot

library(ggplot2)
library(tidyr)

plot_data = mutate(PVE_results_multiMAF, PVEn = PVEsn - PVEs)
plot_data = mutate(plot_data, PVEd = PVEsnd - PVEsn)
plot_data = plot_data[,c(1,2,5,11,12)]
colnames(plot_data)<- c("dist","MAF","self","nei","env")
plot_data = gather(plot_data, key = "type", value = "PVE", -dist, -MAF)

#PVE vs. distance, per MAF
PVE_plot = ggplot(plot_data, aes(x = dist, y = PVE, fill = type)) + geom_bar(position="stack",stat="identity") + facet_grid(~ MAF) + ggtitle("PVE partitioning of each distance") + xlab("distance of included neighbours, per MAF") + ylab("PVE")
ggsave(filename = paste0("results/",species,"/PVE_analysis/multiMAF/extended/",daynumber,"_",poolName,"/",poolName,"_PVE_vs_distance_per_MAF.png"),plot=PVE_plot,width=4000,height=3000,units="px",device="png")
rm(PVE_plot)

#PVE vs. MAF, per distance
PVE_plot = ggplot(plot_data, aes(x = as.factor(MAF), y = PVE, fill = type)) + geom_bar(position="stack",stat="identity") + facet_grid(~ dist) + ggtitle("PVE partitioning of each distance") + xlab("MAF, per distance of included neighbours") + ylab("PVE")
ggsave(filename = paste0("results/",species,"/PVE_analysis/multiMAF/extended/",daynumber,"_",poolName,"/",poolName,"_PVE_vs_MAF_per_distance.png"),plot=PVE_plot,width=4000,height=3000,units="px",device="png")
rm(PVE_plot)

#dist vs. PVE sum, lines, with MAF as colour and type as linestyle
plot_data_3 = PVE_results_multiMAF[,c(1,2,5,7,9)]
colnames(plot_data_3)<- c("dist","MAF","self","self+nei","self+nei+env")
plot_data_3 = gather(plot_data_3, key = "type", value = "PVE", -dist, -MAF)
plot_data_3 = mutate(plot_data_3, MAF_type = paste0(MAF,"_",type))
PVE_plot = ggplot(plot_data_3, aes(x = dist, y = PVE, group = MAF_type, colour = as.factor(MAF), linetype = type)) + geom_line() + geom_point() + ggtitle("PVE partitioning of each distance")
ggsave(filename = paste0("results/",species,"/PVE_analysis/multiMAF/extended/",daynumber,"_",poolName,"/",poolName,"_distance_vs_PVE_sum_per_type_and_MAF.png"),plot=PVE_plot,width=4000,height=3000,units="px",device="png")
rm(PVE_plot)

detach(package:ggplot2)
detach(package:tidyr)

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
d_env = matrix(d_env, nrow = length(distances), ncol = length(MAF))
d_nei = matrix(d_nei, nrow = length(distances), ncol = length(MAF))
d_self = matrix(d_self, nrow = length(distances), ncol = length(MAF))
plot_data_2 = matrix(0.0, ncol = length(MAF), nrow = length(distances))
for(i in 1:nrow(plot_data_2)){
  for(j in 1:ncol(plot_data_2)){
    if(d_env[i,j]>31){d_env[i,j]=31}
    else if(d_env[i,j]<0){d_env[i,j]=0}
    if(d_nei[i,j]>31){d_nei[i,j]=31}
    else if(d_nei[i,j]<0){d_nei[i,j]=0}
    if(d_self[i,j]>31){d_self[i,j]=31}
    else if(d_self[i,j]<0){d_self[i,j]=0}
    plot_data_2[i,j] = (1024*d_env[i,j]) + (32*d_nei[i,j]) + d_self[i,j]
  }
}
png(filename = paste0("results/",species,"/PVE_analysis/multiMAF/extended/",daynumber,"_",poolName,"/",poolName,"_distance_vs_MAF_rgb_map.png"),width=1600,height=800,units="px")
bgc <- par(bg = "black")
par(mar = c(8.0,4.1,4.1,9.1))#bottom,left,top,right
par(col.main = "white",col.axis = "white", col.lab = "white", col = "white")
image(plot_data_2, col = heatmap_cs_def, zlim = c(0,32767), axes = FALSE, x = c(1:nrow(plot_data_2)), y = c(1:ncol(plot_data_2)), xlab = "distance of included neighbours", ylab = "MAF", main = "PVE partitioning, distance vs. MAF") + axis(side = 1, at = c(1:nrow(plot_data_2)),labels = as.factor(distances), las = 2) + axis(side = 4, at = c(1:ncol(plot_data_2)),labels = as.factor(MAF), las = 2) + mtext("black to colour scale, env in red, nei in green, self in blue", col = "white")
dev.off()

#cleaning
rm(bgc)
rm(d_env)
rm(d_nei)
rm(d_self)
rm(plot_data)
rm(plot_data_2)
rm(plot_data_3)
rm(distances)
rm(heatmap_cs_def)
rm(i)
rm(j)
rm(MAF)
rm(getStackedHeatmapColourmapRGB)

