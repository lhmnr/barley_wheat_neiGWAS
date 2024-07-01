library(dplyr)
library(readr)

#this script is to analyse the field size and amount of neighbours

################################# settings

#day number for easier identification of outputs
daynumber = "20240513"

#pool settings:
species = "barley"
#barley or wheat
poolName = "Scald_pooled"
#barley pools:
#NFNB_pooled, PM_2015, Blr_pooled, Scald_pooled, SFNB_pooled
#wheat pools:
#YR_pooled, LR_2017, SR_2017, YLS_2017, STB_2019, YS_2019

#################################

distances = c(1,1.5,2,2.3,2.9,3,3.2,3.7,4,4.2,4.3,4.5,5,5.2,5.4,5.7,5.9,6,6.1,6.4,6.5,6.8,7,7.1,7.25,7.3,7.7,7.9,8,8.1,8.3,8.5,8.6,8.7,8.98,9,9.1,9.3,9.45,9.5,9.9,9.95,10,10.1,10.2,10.3,10.5,10.7,10.8,10.9,11,11.1,11.2,11.4,11.5,11.7,11.8,12,12.05,12.1,12.2,12.3,12.4,12.6,12.7,12.8,12.9,13,13.1,13.2,13.4,13.45,13.5,13.7,13.9,13.95,14,14.1,14.2,14.3,14.4,14.5,14.6,14.8,14.9,15,15.1,15.2,15.25,15.275,15.3,15.55,15.6,15.65,15.7,15.9,16)

smap = read.csv(paste0("results/",species,"/filtered_inputs/data_collection/pheno_",poolName,".csv"))
smap = smap[,c(4,5,9)]

average_neighbours = {}
for(i in 1:length(distances)){
  a = {}
  for(j in 1:nrow(smap)){
    b = 0
    for(k in 1:nrow(smap)){
      if(smap[j,3]==smap[k,3]){
        if(((smap[j,1]-smap[k,1])*(smap[j,1]-smap[k,1]))+((smap[j,2]-smap[k,2])*(smap[j,2]-smap[k,2]))<=distances[i]*distances[i]){
          b=b+1
        }
      }
    }
    a = c(a,b)
  }
  d = sum(a)/length(a)
  smap = mutate(smap, "new_column" = a)
  colnames(smap)<- c(colnames(smap)[1:(ncol(smap)-1)], paste0("neis at ",distances[i]))
  average_neighbours = c(average_neighbours, d)
}

average_neighbours = data.frame("distance" = distances, "average_neighbours" = average_neighbours)

################################# save

#create output folder if necessary
if(!file.exists(paste0("results/",species,"/field_analysis"))){
  dir.create(path = paste0("results/",species,"/field_analysis"))
}
if(!file.exists(paste0("results/",species,"/field_analysis/",daynumber,"_",poolName))){
  dir.create(path = paste0("results/",species,"/field_analysis/",daynumber,"_",poolName))
}

write_csv(smap, file = paste0("results/",species,"/field_analysis/",daynumber,"_",poolName,"/",poolName,"_neighbours.csv"))
write_csv(average_neighbours, file = paste0("results/",species,"/field_analysis/",daynumber,"_",poolName,"/",poolName,"_average_neighbours.csv"))

################################# plot fields

smap = smap[,c(1,2,3)]
experiments = unique(smap[,3])
colnames(smap)<-c("row","range","experiment")

for(k in 1:length(experiments)){
  smap_k = filter(smap, experiment == experiments[k])
  png(filename = paste0("results/",species,"/field_analysis/",daynumber,"_",poolName,"/",poolName,"_experiment",experiments[k],"_field.png"),width=800,height=600,units="px")
  plot(smap_k[,c(1,2)])
  dev.off()
}

rm(a)
rm(b)
rm(d)
rm(i)
rm(j)
rm(k)
rm(experiments)
rm(distances)
rm(smap)
rm(smap_k)
rm(average_neighbours)

