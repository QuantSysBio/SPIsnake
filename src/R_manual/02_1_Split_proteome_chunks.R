library(seqinr)
library(parallel)
library(parallelly)
library(foreach)
source("/home/yhorokh/Snakemake/SPI-snake/src/snakefiles/functions.R")

### ---------------------------- (1) Read input file and extract info ----------------------------
# database
proteome <- as.character("/home/yhorokh/Snakemake/SPI-snake/data/reference/proteome_expressed_gencode.fasta")
dat = read.fasta(file=proteome, 
                 seqtype="AA", as.string = TRUE)

# Keep only proteome name
proteome <- unlist(strsplit(proteome, "/", fixed = T))[grep(".fasta", unlist(strsplit(proteome, "/", fixed = T)))]
proteome <- unlist(strsplit(proteome, ".fasta", fixed = T))[1]

# Clustering
prot_cluster <- read.table("/home/yhorokh/Snakemake/SPI-snake/results/Cluster/proteome_expressed_gencode/proteome_expressed_gencode_cluster.tsv")

# Master table
Master_table <- read.csv("Master_table.csv")

# max intervening sequence length
MiSl <- as.numeric(Master_table$Min_Interv_length[Master_table$Proteome == proteome])
print(paste("max_intervening_length:", MiSl))

# Input proteome filter
min_protein_length = 9

# Proteome chunks
max_length=500
overlap_length=MiSl*2

# Chunk size
maxE = 255319000*1

# CPUs
# Ncpu = availableCores(27)
Ncpu = availableCores()
cl <- parallel::makeForkCluster(Ncpu)
cl <- parallelly::autoStopCluster(cl)

# Output dir
directory = "/home/yhorokh/Snakemake/SPI-snake/results/DB_exhaustive"
suppressWarnings(
  dir.create(directory)
  )
directory <- paste0(directory, "/", proteome)
suppressWarnings(
  dir.create(directory)
)

### ---------------------------- (2) Proteome pre-processing --------------------------------------

# Filter by minimal length and re-order by similarity from clustering
dat <- dat[which(lapply(dat, nchar) >= min_protein_length)]
prot_cluster <- prot_cluster[prot_cluster$V2 %in% names(dat),]
dat <- dat[match(prot_cluster$V2, names(dat))]

# Split long entries into overlapping chunks
dat <- Split_list_max_length_parallel(String_list = dat, 
                                      max_length=max_length, 
                                      overlap_length=overlap_length)


# Estimate chunks
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
Save_prot_chunk(dat=dat, 
                nF=nF, 
                orderedProteomeEntries=orderedProteomeEntries, 
                directory=directory)
