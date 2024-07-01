library(rNeighborGWAS)
library(gaston)
library(dplyr)
library(readr)
fn = system("ls ./pheno/wheat/*.csv", intern=TRUE)

#this script is to explore the data for different distances and adjust minor allele cutoff

#finding: calc_PVEnei, delta_PVE does not give the distance with the highest peaks,
#therefore automated neiGWAS should be done over all distances
#also since different biological mechanisms will act up to different ranges

################################# settings
#you need to adjust this before running

#pool settings:
poolName = "YR_pooled"
pool = c(4,6)
#individual datasets: check fn
#pools:
#YR: 4,6
#LR_2017: 1
#SR_2017: 2
#YLS_2017: 3
#STB_2019: 5
#YS_2019: 7

#MAF settings:
maf_filter = 0.01
#adjust such that no weird structures appear in manhattan-plot

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

geno = t(geno)
geno = geno - 1

#smap
smap = cbind(pheno$Range,pheno$Row)

#gmap
gmap = read.csv("./geno/wheat/positions.csv")
gmap = filter(gmap,!rare_allele)

################################# manual neiGWAS
#do a single neiGWAS of desired distance, run this after section "loading and preparing data"

# choose effect distance
#distances up to 2.9:
#distances = c(1,sqrt(2)+0.01,2,sqrt(5)+0.01,sqrt(8)+0.01)
#distances up to 4.3:
#distances = c(1,sqrt(2)+0.01,2,sqrt(5)+0.01,sqrt(8)+0.01,3,sqrt(10)+0.01,sqrt(13)+0.01,4,sqrt(17)+0.01,sqrt(18)+0.01)
#distances up to 5.9:
distances = c(1,sqrt(2)+0.01,2,sqrt(5)+0.01,sqrt(8)+0.01,3,sqrt(10)+0.01,sqrt(13)+0.01,4,sqrt(17)+0.01,sqrt(18)+0.01,sqrt(20)+0.01,5,sqrt(26)+0.01,sqrt(29)+0.01,sqrt(32)+0.01,sqrt(34)+0.01)
#distances up to 10:
#distances = c(1,sqrt(2)+0.01,2,sqrt(5)+0.01,sqrt(8)+0.01,3,sqrt(10)+0.01,sqrt(13)+0.01,4,sqrt(17)+0.01,sqrt(18)+0.01,sqrt(20)+0.01,5,sqrt(26)+0.01,sqrt(29)+0.01,sqrt(32)+0.01,sqrt(34)+0.01,6,sqrt(37)+0.01,sqrt(40)+0.01,sqrt(41)+0.01,sqrt(45)+0.01,7,sqrt(50)+0.01,sqrt(52)+0.01,sqrt(53)+0.01,sqrt(58)+0.01,sqrt(61)+0.01,8,sqrt(65)+0.01,sqrt(68)+0.01,sqrt(72)+0.01,sqrt(73)+0.01,sqrt(74)+0.01,sqrt(80)+0.01,9,sqrt(82)+0.01,sqrt(85)+0.01,sqrt(89)+0.01,sqrt(90)+0.01,sqrt(97)+0.01,sqrt(98)+0.01,10)
res = calc_PVEnei(geno=geno,pheno=pheno$Damage_Level,smap=smap,scale=distances,addcovar=as.factor(pheno$Experiment_Number),grouping=pheno$Experiment_Number)
res$PVEself+res$PVEnei
res

delta_PVE(res)

#notes:

#none yet

#set scale for neiGWAS
scale <- 2.0

#neiGWAS
gwas_out <- neiGWAS(geno=geno,pheno=pheno$Damage_Level,
                    gmap=gmap, smap=smap,
                    scale=scale,
                    addcovar=as.factor(pheno$Experiment_Number),
                    grouping=pheno$Experiment_Number
)

#manhattan plot with line, QQ-plot
gaston::manhattan(gwas_out)
abline(h = -log10(0.05/ncol(geno))) #bonferroni-corrected significance line, may not shows if no point close or above
-log10(0.05/ncol(geno)) #print the threshold
max(-log10(gwas_out$p)) #prints the p-value top SNP
gaston::qqplot.pvalues(gwas_out$p)
