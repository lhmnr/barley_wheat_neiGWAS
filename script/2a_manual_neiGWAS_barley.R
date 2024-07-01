library(rNeighborGWAS)
library(gaston)
library(dplyr)
library(readr)
fn = system("ls ./pheno/barley/*.csv", intern=TRUE)

#this script is to explore the data for different distances and adjust minor allele cutoff

#finding: calc_PVEnei, delta_PVE does not give the distance with the highest peaks,
#therefore automated neiGWAS should be done over all distances
#also since different biological mechanisms will act up to different ranges

################################# settings
#you need to adjust this before running

#pool settings:
poolName = "NFNB_pooled"
pool = c(1,6,10)
#individual datasets: check fn
#pools:
#NFNB: 1,6,10
#PM_2015: 2
#Blr: 5,9
#Scald: 3,7,11
#SFNB: 4,8,12

#MAF settings:
maf_filter = 0.01
#adjust such that no weird structures appear in manhattan-plot

################################# loading and preparing data
#before running check manually in folders if all required files are present

#loading all experiments in a loop, appending each to pheno
#note: Experiment_Number in pheno equals position in fn
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
#maf_filter is set in settings at top of script
af = rowSums(geno)/(2*ncol(geno)) #assuming geno is complete
rare_allele = ifelse(af<maf_filter|1-af<maf_filter,TRUE,FALSE)
geno = filter(geno,!rare_allele)

geno = t(geno)
geno = geno - 1

#smap
smap = cbind(pheno$Range,pheno$Row)

#gmap
gmap = read.csv("./geno/barley/positions.csv")
gmap = filter(gmap,!rare_allele)

################################# manual neiGWAS
#do a single neiGWAS of desired distance

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

#single experiment GWAS:
#experiment 1 (2015 NFNB): distance peaks at 2, but delta_PVE chose 1, use 1.5
#2 no sig SNPs, 1 weird line and QQ off, 1.5 sig SNPs at end of chr3
#experiment 2 (2015 PM): distance best at 2 or 3
#experiment 3 (2015 Scald): distance best at 1
#experiment 4 (2015 SFNB): PVE_nei low, seems best at 1.5 or 2
#no sig, but consistent positions at beg of chr2, end of chr5 and end of chr7
#experiment 5 (2016 BLR): distance seems best at 1.5, delta_PVE chose 1
#experiment 6 (2016 NFNB): PVE_nei low, distance best at 3.2
#experiment 7 (2016 Scald): distance seems best at 1 or 3 or 3.2
#there are two smaller peaks at 1 and 2, but might be random
#experiment 8 (2016 SFNB): distance peaks at 2, but delta_PVE chose 1
#experiment 9 (2017 BLR): distance peaks at 5.2 or 5.4, delta_PVE chose 5
#smaller peaks at 1 and 4.3
#experiment 10 (2017 NFNB): distance peaks at 1
#experiment 11 (2017 Scald): distance peaks at 1
#experiment 12 (2017 SFNB): distance peaks at 2, but delta_PVE chose 1

#pooled experiment GWAS:
#NFNB (1,6,10): distance best at 2, delta_PVE chose 1
#sig SNP in chr5 at 2, also high at 2.3
#PM (2): no pooling
#Blr (5,9): distance best at 4, delta_PVE chose 1, plateau at 5.2, use 2
#best SNP found at 2
#Scald (3,7,11): distancce best at 1, peaks at 5, peak at 2, may try 3
#almost sig SNPs 2.3 and 3, also quite at 1, 2, 2.9, 4, 4.2, 4.5 and 5
#SFNB (4,8,12):

#set scale for neiGWAS
scale <- 3.0

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
