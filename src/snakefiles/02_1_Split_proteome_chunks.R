### ---------------------------------------------- Split proteome chunks  ----------------------------------------------
# description:  Sort input proteomes according to Linclust output. Split long proteins into overlapping chunks.
#               
# input:        1. Proteome .fasta
#               2. Linclust .tsv table with protein order
#               3. Parameters: min/max sequence length, chunk size
# output:       
#               - Proteome spit into chunks of fixed size, similar sequences stored together
#               
# author:       YH, JL


library(dplyr)
library(seqinr)
library(parallel)
library(parallelly)
library(foreach)
source(snakemake@input[["functions"]])

print("Loaded functions. Loading the data")

### ---------------------------- (1) read input file and extract info --------------------------------------
### Inputs:
# Protein database
proteome <- as.character(snakemake@input[["proteome"]])
dat = read.fasta(file=proteome, seqtype="AA", as.string = TRUE)

# Keep only proteome name
proteome <- unlist(strsplit(proteome, "/", fixed = T))[grep(".fasta", unlist(strsplit(proteome, "/", fixed = T)))]
proteome <- unlist(strsplit(proteome, ".fasta", fixed = T))[1]

# Clustering
prot_cluster <- read.table(snakemake@input[["prot_cluster"]])

# Master table
Master_table <- read.csv("Master_table.csv") %>%
 as_tibble()  %>%
 filter(Proteome == proteome) %>%
 select(Proteome, MaxE, Min_Interv_length)


# max intervening sequence length
MiSl <- as.numeric(Master_table$Min_Interv_length)
print(paste("max_intervening_length:", MiSl))

### Params:
# Input proteome filter
min_protein_length = snakemake@params[["min_protein_length"]]
print(paste("min_protein_length:", min_protein_length))
dat <- dat[which(lapply(dat, nchar) >= min_protein_length)]

# Proteome chunks
max_length=as.numeric(snakemake@params[["max_protein_length"]])
overlap_length=MiSl*2

# Chunk size
maxE = Master_table$MaxE
print(paste("maxE:", maxE))

### CPUs
Ncpu = as.numeric(snakemake@params[["n"]])
#Ncpu = availableCores()
#cl <- parallel::makeForkCluster(Ncpu)
#cl <- parallelly::autoStopCluster(cl)

# Output dir
directory = snakemake@params[["dir_DB_exhaustive"]]
directory <- paste0(directory, "/", proteome)
suppressWarnings(
  dir.create(directory)
)

# Log
log <- file(snakemake@log[[1]], open="wt")
sink(log)


### ---------------------------- (2) Proteome pre-processing --------------------------------------
# Filter by minimal length and re-order by similarity from clustering
dat <- dat[which(lapply(dat, nchar) >= min_protein_length)]
prot_cluster <- prot_cluster[prot_cluster$V2 %in% names(dat),]
dat <- dat[match(prot_cluster$V2, names(dat))]

# Split long entries into overlapping chunks
dat <- Split_list_max_length_parallel(String_list = dat, 
                                      max_length=max_length, 
                                      overlap_length=unique(overlap_length))


# Estimate chunks
for (i in 1:nrow(Master_table)) {
  maxE <- Master_table$MaxE[i]
  
  {
    orderedProteomeEntries <- names(dat)
    
    L = unlist(lapply(dat, nchar))
    numPSP = L*350
    cumNum = rep(NA,length(L))
    cumNum[1] = numPSP[1]
    for(i in 2:length(L)){
      
      cumNum[i] = cumNum[i-1]+numPSP[i]
    }
    
    # find intervals for this Nmer Pi
    b = 1
    index = 0
    Pi = 0
    
    while(index<length(L)){
      index = which(cumNum==max(cumNum[which(cumNum<=b*maxE)]))
      Pi = c(Pi,index)
      b = b+1
      print(index)
    }
    
    # Create fasta files for every chunk
    nF = length(Pi)-1
  }
  
  # Export .fasta
  suppressWarnings(dir.create("/home/yhorokh/Snakemake/SPI-snake_dev/results/DB_exhaustive/"))
  suppressWarnings(dir.create(directory))
  Save_prot_chunk(dat=dat, 
                  nF=nF, 
                  orderedProteomeEntries=orderedProteomeEntries, 
                  directory=directory,
                  maxE=maxE)
}
