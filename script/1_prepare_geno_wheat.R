library(vcfR)

# g = read.vcfR("./geno/wheat/caigewheat__10511variants__801individuals.vcf")
g = read.vcfR("./geno/wheat/caigewheat__10511variants__801individuals_imp.vcf.gz")

search_dup = function(i) {
  target = which(paste0(g@fix[i,1],g@fix[i,2])==paste0(g@fix[,1],g@fix[,2]))
  if(length(target)>1){
    print(paste0(g@fix[i,1],g@fix[i,2]))
  }
}

dup_id = mapply(search_dup,1:nrow(g@fix))

geno = g@gt[,-1]
geno[geno=="0|0"] = 0
geno[geno=="1|1"] = 2
geno[geno=="1|0"] = 1
geno[geno=="0|1"] = 1
geno = as.data.frame(geno)
write.csv(geno,"./geno/wheat/geno.csv",row.names=FALSE)

fix = g@fix
fix = as.data.frame(fix[,1:2])
write.csv(fix,"./geno/wheat/positions.csv",row.names=FALSE)


