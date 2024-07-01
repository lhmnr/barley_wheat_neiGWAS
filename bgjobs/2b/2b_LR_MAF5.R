library(rNeighborGWAS)
library(gaston)
library(dplyr)
library(readr)
fn = system("ls ./pheno/wheat/*.csv", intern=TRUE)

#this script is to prepare the input data for automated neiGWAS and further analysis
#additionally, it extracts haplo-genotypes for better LD

################################# settings
#you need to adjust this before running

#day number for easier identification of outputs
daynumber = "20240425"
#do not use slashes or points, since this will be used for the output folder name

#pool settings:
poolName = "LR_2017-MAF5"
pool = c(1)
#individual datasets: check fn
#pools:
#YR_pooled: 4,6
#LR_2017: 1
#SR_2017: 2
#YLS_2017: 3
#STB_2019: 5
#YS_2019: 7

#MAF settings:
maf_filter = 0.05
#for automatic neiGWAS exploration, use 0.01

################################# loading and preparing data
#before running check manually in folders if all required files are present

#loading all experiments in a loop, appending each to pheno
#note: Experiment_Number in pheno equals position in fn
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
#maf_filter is set in settings at top of script
af = rowSums(geno)/(2*ncol(geno)) #assuming geno is complete
rare_allele = ifelse(af<maf_filter|1-af<maf_filter,TRUE,FALSE)
geno = filter(geno,!rare_allele)
af = af[!rare_allele]

geno = t(geno)
geno = geno - 1

#smap
smap = cbind(pheno$Range,pheno$Row)

#gmap
gmap = read.csv("./geno/wheat/positions.csv")
gmap = filter(gmap,!rare_allele)

################################# get haplotypes
library(vcfR)

g = read.vcfR("./geno/wheat/caigewheat__10511variants__801individuals_imp.vcf.gz")

genotypes = g@gt[,-1]
genotypes = as.data.frame(genotypes)
rm(g)

genotypes = genotypes[,andIDs]
genotypes = filter(genotypes,!rare_allele)
genotypes = mutate (genotypes, "chr" = gmap$CHROM)
genotypes = mutate (genotypes, "pos" = gmap$POS)
genotypes = select(genotypes, c((ncol(genotypes)-1),ncol(genotypes),1:(ncol(genotypes)-2)))
haplotypes = data.frame("chr"=gmap$CHROM,"pos"=gmap$POS)
for(i in 3:(ncol(genotypes)-1)){
  genotypes_i = genotypes[,c(1,2,i)]
  colnames(genotypes_i)<- c("chr","pos","diplo")
  genotypes_i = mutate(genotypes_i, haplo_1 = sapply(strsplit(diplo,split="|",fixed=TRUE), function(x) (x[1])))
  genotypes_i = mutate(genotypes_i, haplo_2 = sapply(strsplit(diplo,split="|",fixed=TRUE), function(x) (x[2])))
  genotypes_i = genotypes_i[,c(1,2,4,5)]
  colnames(genotypes_i)<- c("chr","pos",((2*i)-1),(2*i))
  haplotypes = merge(haplotypes,genotypes_i)
}
haplotypes = arrange(haplotypes, pos)
haplotypes = arrange(haplotypes, chr)

rm(genotypes)
rm(genotypes_i)
rm(i)
################################# save data

#create output folder if necessary
if(!file.exists("results")){
  dir.create(path = "results")
}
if(!file.exists("results/wheat")){
  dir.create(path = "results/wheat")
}
if(!file.exists("results/wheat/filtered_inputs")){
  dir.create(path = "results/wheat/filtered_inputs")
}
if(!file.exists(paste0("results/wheat/filtered_inputs/",daynumber,"_",poolName))){
  dir.create(path = paste0("results/wheat/filtered_inputs/",daynumber,"_",poolName))
}
if(!file.exists("results/wheat/filtered_inputs/data_collection")){
  dir.create(path = "results/wheat/filtered_inputs/data_collection")
}

geno = as.data.frame(geno)
smap = as.data.frame(smap)
af = as.data.frame(af)

#export data
write_csv(geno, file = paste0("results/wheat/filtered_inputs/",daynumber,"_",poolName,"/geno_",poolName,".csv"))
write_csv(pheno, file = paste0("results/wheat/filtered_inputs/",daynumber,"_",poolName,"/pheno_",poolName,".csv"))
write_csv(gmap, file = paste0("results/wheat/filtered_inputs/",daynumber,"_",poolName,"/gmap_",poolName,".csv"))
write_csv(smap, file = paste0("results/wheat/filtered_inputs/",daynumber,"_",poolName,"/smap_",poolName,".csv"))
write_csv(af, file = paste0("results/wheat/filtered_inputs/",daynumber,"_",poolName,"/af_",poolName,".csv"))
write_csv(haplotypes, file = paste0("results/wheat/filtered_inputs/",daynumber,"_",poolName,"/haplo_",poolName,".csv"))

#copy these files manually to "results/wheat/filtered_inputs/data_collection" to use them for automated neiGWAS and further analysis

