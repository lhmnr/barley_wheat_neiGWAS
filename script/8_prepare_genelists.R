library(dplyr)
library(readr)

#this script is to prepare the gene lists

################################# loading and preparing data
fn = system("ls ./genelists/raw_tables/*.csv", intern=TRUE)

genelist = data.frame(NULL)
subloci_without_transcripts = data.frame("locus"=character(),"sublocus"=character(),"sublocus_pos_l"=character(),"sublocus_pos_r"=character())
for(i in 1:length(fn)){
  fn_i = fn[i]#i
  genelist_i = read_delim(fn_i, delim = "\t", escape_double = FALSE)
  if(nrow(genelist_i)>0){
    genelist_i = select(genelist_i, c(2,4,5,6,7,8,9))
    
    genelist_i = mutate(genelist_i, taxa = sapply(strsplit(SNP_taxa_chr_pos,split=" ",fixed=TRUE), function(x) (x[1])))
    genelist_i = mutate(genelist_i, SNP_chr = sapply(strsplit(SNP_taxa_chr_pos,split=" ",fixed=TRUE), function(x) (x[2])))
    genelist_i = mutate(genelist_i, SNP_pos = sapply(strsplit(SNP_taxa_chr_pos,split=" ",fixed=TRUE), function(x) (x[3])))
    genelist_i = mutate(genelist_i, SNP_chr_pos = paste0(SNP_chr,"_",SNP_pos))
    genelist_i = mutate(genelist_i, transcript_chr = sapply(strsplit(gene_chr_pos,split=":",fixed=TRUE), function(x) (x[1])))
    genelist_i = mutate(genelist_i, transcript_pos = sapply(strsplit(gene_chr_pos,split=":",fixed=TRUE), function(x) (x[2])))
    genelist_i = mutate(genelist_i, transcript_pos_l = sapply(strsplit(transcript_pos,split="-",fixed=TRUE), function(x) (x[1])))
    genelist_i = mutate(genelist_i, transcript_pos_r = sapply(strsplit(transcript_pos,split="-",fixed=TRUE), function(x) (x[2])))
    genelist_i = select(genelist_i, c(8,11,2,3,4,5,6,12,14,15))
    
    genelist_i = mutate(genelist_i, gene_ID = sapply(strsplit(link,split="t=",fixed=TRUE), function(x) (x[2])))
    genelist_i = mutate(genelist_i, transcript_ID = sapply(strsplit(gene_ID,split=".",fixed=TRUE), function(x) (x[length(x)])))
    genelist_i = mutate(genelist_i, gene_ID = substring(gene_ID, 1, nchar(gene_ID)-(nchar(transcript_ID)+1)))
    genelist_i = select(genelist_i, c(1,2,11,12,8,9,10,4,5,6,7,3))
    
    fn_i = strsplit(fn_i,split="/",fixed=TRUE)[[1]][4]
    fn_i = strsplit(fn_i,split=".",fixed=TRUE)[[1]][1]
    locus_i = strsplit(fn_i,split="_",fixed=TRUE)[[1]][2]
    range_l_i = strsplit(fn_i,split="_",fixed=TRUE)[[1]][3]
    range_r_i = strsplit(fn_i,split="_",fixed=TRUE)[[1]][4]
    genelist_i = mutate(genelist_i, locus = locus_i)
    genelist_i = mutate(genelist_i, sublocus = i)
    genelist_i = mutate(genelist_i, sublocus_pos_l = range_l_i)
    genelist_i = mutate(genelist_i, sublocus_pos_r = range_r_i)
    
    genelist_i = mutate(genelist_i, transcript_length = as.numeric(transcript_pos_r) - as.numeric(transcript_pos_l))
    genelist_i = add_count(genelist_i, gene_ID, name = "n_transcripts")
    count_overlaps <- function(A,B){
      C = genelist_i$transcript_pos_l
      D = genelist_i$transcript_pos_r
      G = genelist_i$gene_ID
      overlaps = {}
      for(k in c(1:length(C))){
        if(!(C[k]>B|D[k]<A)){
          overlaps = c(overlaps,G[k])
        }
      }
      k2 = length(unique(overlaps))
      return(k2)
    }
    genelist_i = mutate(genelist_i, genes_overlap = mapply(count_overlaps, transcript_pos_l, transcript_pos_r))
    
    if(nrow(genelist)==0){
      genelist = genelist_i
    }
    else{
      genelist = rows_append(genelist, genelist_i)
    }
  }
  else{
    fn_i = strsplit(fn_i,split="/",fixed=TRUE)[[1]][4]
    fn_i = strsplit(fn_i,split=".",fixed=TRUE)[[1]][1]
    locus_i = strsplit(fn_i,split="_",fixed=TRUE)[[1]][2]
    range_l_i = strsplit(fn_i,split="_",fixed=TRUE)[[1]][3]
    range_r_i = strsplit(fn_i,split="_",fixed=TRUE)[[1]][4]
    empty_sublist_i = data.frame("locus"=locus_i, "sublocus" = as.character(i), "sublocus_pos_l" = range_l_i, "sublocus_pos_r" = range_r_i)
    subloci_without_transcripts = rows_append(subloci_without_transcripts, empty_sublist_i)
  }
}

#correction for snoRNA in wheat
for(j in which(genelist$gene_ID=="")){
  genelist[j,3] = strsplit(genelist[j,4][[1]],split="-",fixed=TRUE)[[1]][1]
  genelist[j,4] = strsplit(genelist[j,4][[1]],split="-",fixed=TRUE)[[1]][2]
}

#clean
rm(empty_sublist_i)
rm(genelist_i)
rm(fn_i)
rm(i)
rm(j)
rm(locus_i)
rm(range_l_i)
rm(range_r_i)
rm(count_overlaps)

################################# some quick analysis

#amount of unique candidate genes is
length(unique(genelist$gene_ID))

#confirm
length(unique((distinct(genelist[,c(3,4)]))$gene_ID))
#amount of unique transcripts "1" of candidate genes
length(which((distinct(genelist[,c(3,4)]))$transcript_ID==1))
#if this number is lower than the last one, there are genes without a transcript "1"

#amount of unique transcripts
length(which(!duplicated(genelist[,c(3,4)])))

#how many unique transcripts per unique gene
table(genelist[which(!duplicated(genelist$gene_ID)),]$n_transcripts)
#genes with multiple transcripts are rare (<5% of genes)

#how are the unique genes distributed
table(genelist[which(!duplicated(genelist$gene_ID)),]$transcript_chr)
table(genelist[which(!duplicated(genelist[,c(3,13)])),]$locus)
#if the following line gives lower values, some genes are found in multiple loci
table(genelist[which(!duplicated(genelist$gene_ID)),]$locus)
#this expected to happen in "w4a1b", as some subloci overlap with "w4a1a"
#note that subloci, which are identical in "a" and "b" are put as "w4a1"

##plot
##barplot(table(genelist[which(!duplicated(genelist$gene_ID)),]$transcript_chr), main = "candidate genes per chromosome")
##barplot(table(genelist[which(!duplicated(genelist[,c(3,13)])),]$locus), main = "candidate genes per locus")
#d1 = table(genelist[which(!duplicated(genelist$gene_ID)),]$transcript_chr)
#png(filename = paste0("genelists/collection/gene_distribution_chromosome.png"),width=1600,height=800,units="px")
#d2 = barplot(d1, main = "candidate genes per chromosome")
#text(d2,d1/2+2,labels=d1)
#dev.off()
#d3 = table(genelist[which(!duplicated(genelist[,c(3,13)])),]$locus)
#png(filename = paste0("genelists/collection/gene_distribution_locus.png"),width=1600,height=800,units="px")
#d4 = barplot(d3, main = "candidate genes per locus")
#text(d4,d3/2+2,labels=d3)
#dev.off()
#
##clean
#rm(d1)
#rm(d2)
#rm(d3)
#rm(d4)

################################# save data

#create output folder if necessary
if(!file.exists("genelists/collection")){
  dir.create(path = "genelists/collection")
}

locus_info = data.frame("locus"=character(),"taxa"=character(),"chr"=character(),"locus_size"=numeric(),"n_genes"=numeric(),"n_transcripts"=numeric(),"gene_density"=numeric(),"n_subloci"=numeric())

#get loci
loci = unique(genelist$locus)

#separate loci
for(i in 1:length(loci)){
  locus_i = loci[i]
  #arranging the intersecting loci
  if(locus_i=="w4a1a"){
    genelist_i = filter(genelist, (locus == "w4a1a")|(locus == "w4a1"))
  } else if(locus_i=="w4a1b"){
    genelist_i = filter(genelist, (locus == "w4a1b")|(locus == "w4a1"))
  } else if(locus_i=="w4a1"){
    genelist_i = filter(genelist, (locus == "w4a1")|(locus == "w4a1a")|(locus == "w4a1b"))
  } else{#all other loci
    genelist_i = filter(genelist, locus == locus_i)
  }
  
  #get locus properties
  subloci = genelist_i[which(!duplicated(genelist_i[,14])),]
  subloci = subloci[,c(14,15,16)]
  colnames(subloci)<- c("sublocus_ID","sublocus_pos_l","sublocus_pos_r")
  d6 = table(genelist_i[which(!duplicated(genelist_i$gene_ID)),]$sublocus)
  subloci = mutate(subloci, "n_genes" = sapply(d6[as.character(sublocus_ID)], function(x) (x[[1]])))
  if(locus_i=="w4a1a"){
    for(j in 1:nrow(subloci_without_transcripts)){
      if(subloci_without_transcripts$locus[j] == "w4a1" | subloci_without_transcripts$locus[j] == "w4a1a"){
        sublocus_j = data.frame("sublocus_ID"=as.numeric(subloci_without_transcripts[j,2]),"sublocus_pos_l"=subloci_without_transcripts[j,3],"sublocus_pos_r"=subloci_without_transcripts[j,4],"n_genes"=0)
        subloci = add_row(subloci, sublocus_j)
      }
    }
  } else if(locus_i=="w4a1b"){
    for(j in 1:nrow(subloci_without_transcripts)){
      if(subloci_without_transcripts$locus[j] == "w4a1" | subloci_without_transcripts$locus[j] == "w4a1b"){
        sublocus_j = data.frame("sublocus_ID"=as.numeric(subloci_without_transcripts[j,2]),"sublocus_pos_l"=subloci_without_transcripts[j,3],"sublocus_pos_r"=subloci_without_transcripts[j,4],"n_genes"=0)
        subloci = add_row(subloci, sublocus_j)
      }
    }
  } else if(locus_i=="w4a1"){
    for(j in 1:nrow(subloci_without_transcripts)){
      if(subloci_without_transcripts$locus[j] == "w4a1" | subloci_without_transcripts$locus[j] == "w4a1a" | subloci_without_transcripts$locus[j] == "w4a1b"){
        sublocus_j = data.frame("sublocus_ID"=as.numeric(subloci_without_transcripts[j,2]),"sublocus_pos_l"=subloci_without_transcripts[j,3],"sublocus_pos_r"=subloci_without_transcripts[j,4],"n_genes"=0)
        subloci = add_row(subloci, sublocus_j)
      }
    }
  } else{
    for(j in 1:nrow(subloci_without_transcripts)){
      if(subloci_without_transcripts$locus[j] == locus_i){
        sublocus_j = data.frame("sublocus_ID"=as.numeric(subloci_without_transcripts[j,2]),"sublocus_pos_l"=subloci_without_transcripts[j,3],"sublocus_pos_r"=subloci_without_transcripts[j,4],"n_genes"=0)
        subloci = add_row(subloci, sublocus_j)
      }
    }
  }
  subloci = mutate(subloci, "size" = as.numeric(sublocus_pos_r) - as.numeric(sublocus_pos_l))
  subloci = select(subloci, c(1,2,3,5,4))
  write_csv(subloci, file = paste0("genelists/collection/subloci_",locus_i,".csv"))
  locus_size = sum(subloci$size)
  if(locus_i=="w4a1"){
    #correction for overlapping subloci of w4a1a and w4a1b, that meet in w4a1
    subloci_2 = mutate(as.data.frame(subloci), "sublocus_pos_l" = as.numeric(sublocus_pos_l))
    subloci_2 = mutate(subloci_2, "sublocus_pos_r" = as.numeric(sublocus_pos_r))
    subloci_2 = mutate(subloci_2, "valid" = TRUE)
    subloci_2 = arrange(subloci_2, sublocus_pos_r)
    subloci_2 = arrange(subloci_2, sublocus_pos_l)
    if(nrow(subloci_2)>1){
      for(j in 1:(nrow(subloci_2)-1)){
        if(subloci_2[j,6]){
          k = 1
          while((subloci_2[j,3]>=subloci_2[j+k,2])&(k>0)){
            subloci_2[j,3] = subloci_2[j+k,3]
            subloci_2[j+k,6] = FALSE
            k = k+1
            if(j+k>nrow(subloci_2)){k=0}
          }
        }
      }
    }
    subloci_2 = filter(subloci_2, valid)
    subloci_2 = mutate(subloci_2, "size" = sublocus_pos_r - sublocus_pos_l)
    locus_size = sum(subloci_2$size)
    rm(subloci_2)
    rm(j)
    rm(k)
  }
  
  if(locus_i=="w4a1"){
    genelist_i = genelist_i[which(!duplicated(genelist_i[,c(3,4)])),]
  }
  #save
  write_csv(genelist_i, file = paste0("genelists/collection/genelist_",locus_i,".csv"))
  
  n_transcripts = nrow(genelist_i)
  n_genes = length(unique(genelist_i$gene_ID))
  n_subloci = nrow(subloci)
  gene_density = n_genes / locus_size
  locus_taxa = genelist_i[1,1][[1]]
  locus_chr = genelist_i[1,5][[1]]
  
  locus_info_i = data.frame("locus"=locus_i,"taxa"=locus_taxa,"chr"=locus_chr,"locus_size"=locus_size,"n_genes"=n_genes,"n_transcripts"=n_transcripts,"gene_density"=gene_density,"n_subloci"=n_subloci)
  locus_info = rows_append(locus_info, locus_info_i)
}

#save to file
write_csv(locus_info, file = paste0("genelists/collection/locus_info.csv"))

#clean
rm(genelist_i)
rm(subloci)
rm(subloci_without_transcripts)
rm(d6)
rm(locus_size)
rm(n_transcripts)
rm(n_genes)
rm(n_subloci)
rm(gene_density)
rm(locus_taxa)
rm(locus_chr)
rm(i)
rm(j)
rm(locus_i)
rm(locus_info_i)
rm(sublocus_j)

################################# plots

d1 = table(genelist[which(!duplicated(genelist$gene_ID)),]$transcript_chr)
png(filename = paste0("genelists/collection/gene_distribution_chromosome.png"),width=1600,height=800,units="px")
d2 = barplot(d1, main = "candidate genes per chromosome")
text(d2,d1/2+2,labels=d1)
dev.off()

d3 = locus_info$n_genes
png(filename = paste0("genelists/collection/gene_distribution_locus.png"),width=1600,height=800,units="px")
d4 = barplot(d3, names.arg = locus_info$locus, main = "candidate genes per locus")
text(d4,d3/2+2,labels=d3)
dev.off()

d5 = locus_info$n_transcripts
png(filename = paste0("genelists/collection/transcript_distribution_locus.png"),width=1600,height=800,units="px")
d6 = barplot(d5, names.arg = locus_info$locus, main = "transcripts per locus")
text(d6,d5/2+2,labels=d5)
dev.off()

d7 = locus_info$locus_size / 1000
png(filename = paste0("genelists/collection/locus_size.png"),width=1600,height=800,units="px")
d8 = barplot(d7, names.arg = locus_info$locus, main = "size of loci compared (in kb)")
text(d8,d7/2+2,labels=d7)
dev.off()

#clean
rm(d1)
rm(d2)
rm(d3)
rm(d4)
rm(d5)
rm(d6)
rm(d7)
rm(d8)

