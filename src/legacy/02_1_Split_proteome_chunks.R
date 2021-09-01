library(seqinr)
library(parallel)
library(parallelly)
library(foreach)
source("/home/yhorokh/Snakemake/SPI-snake/src/functions.R")

### ---------------------------- (1) Read input file and extract info ----------------------------

# mkdir tmp &&
#   cat data/reference/*.fasta > data/reference/tmp.fasta &&
#   mmseqs easy-linclust data/reference/tmp.fasta results tmp &&
#   rm data/reference/tmp.fasta

# Input proteome filter
min_protein_length = 9

# Nmers
Nmers = as.numeric(8:15)

# max intervening sequence length
MiSl = as.numeric(25)

# Proteome chunks
max_length=1000
overlap_length=MiSl*2

# Chunk size
maxE = 255319000*3

# CPUs
# Ncpu = availableCores(27)
Ncpu = availableCores()
cl <- parallel::makeForkCluster(Ncpu)
cl <- parallelly::autoStopCluster(cl)

# database
dat = read.fasta(file="/home/yhorokh/Snakemake/SPI-snake/data/reference/proteome_expressed_gencode.fasta", 
                 seqtype="AA", as.string = TRUE)


# Output dir
directory = "/home/yhorokh/Snakemake/SPI-snake/results/exhaustive_DB"
suppressWarnings(
  dir.create(directory)
  )

### ---------------------------- (2) Proteome pre-processing --------------------------------------
# Filter by minimal length
dat <- dat[which(lapply(dat, nchar) >= min_protein_length)]

# Split long entries into overlapping chunks
dat <- Split_list_max_length_parallel(String_list = dat, 
                                      max_length=max_length, 
                                      overlap_length=overlap_length)


### ---------------------------- (3) Predict runtime info for each CPU ----------------------------
L = rep(NA, length(dat))
for(i in 1:length(L)){
  L[i] = nchar(dat[[i]])
}

index = order(L)
L = L[index]

# estimate runtime function
RT = rep(NA,length(L))
RT[1] = L[1]**1

for(j in 2:length(L)){
  RT[j] = RT[j-1]+L[j]**1
}

# get Intervals with equal estimated runtime
RTmax = RT[length(L)]
q = seq(0,RTmax,(RTmax/Ncpu))
P = rep(NA,Ncpu+1)
P[1] = 0
for(i in 1:Ncpu){
  P[i+1] = which(abs(RT-q[i+1])==min(abs(RT-q[i+1])))
}

# check if any node has no jobs
if(length(P)>length(unique(P))){
  P = unique(P)
  print(paste("Using only first",length(P)-1,"nodes."))
} else {
  print(paste("Using all",length(P)-1,"nodes."))
}

Ncpu = length(P)-1
orderedProteomeEntries = attributes(dat)$name[index]

### ---------------------------- (4) Output: proteome entries ----------------------------

save(orderedProteomeEntries,
     file=unlist(snakemake@output[["orderedProteomeEntries"]]))
save(P,file=unlist(snakemake@output[["P"]]))

### ---------------------------- (5) Chunk proteome and save ----------------------------

print(paste("MaxE:", maxE))
nmer = Nmers[1]
Pi = 0

numPSP = L*350
cumNum = rep(NA,length(L))
cumNum[1] = numPSP[1]
for(i in 2:length(L)){
  
  cumNum[i] = cumNum[i-1]+numPSP[i]
}

# find intervals for this Nmer Pi
b = 1
index = 0

while(index<length(L)){
  index = which(cumNum==max(cumNum[which(cumNum<=b*maxE)]))
  Pi = c(Pi,index)
  b = b+1
  print(index)
}

# Create fasta files for every chunk
nF = length(Pi)-1

Save_prot_chunk(proteome=dat, 
                Nmers=Nmers, 
                nF=nF, 
                MiSl=MiSl, 
                orderedProteomeEntries=orderedProteomeEntries, 
                directory=directory)
