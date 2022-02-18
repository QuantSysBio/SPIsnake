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

### Log
log <- file(snakemake@log[[1]], open="wt")
sink(log)

suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(bettermc))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(stringi))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(parallelly))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(vroom))

### Create temporary directory for vroom
{
  Sys.getenv("TMPDIR") %>% print()
  Sys.getenv("VROOM_TEMP_PATH") %>% print()
  
  vroom_dir = "/tmp/vroom"
  suppressWarnings(dir.create(vroom_dir))
  Sys.setenv(VROOM_TEMP_PATH = vroom_dir)
  Sys.getenv("VROOM_TEMP_PATH") %>%print()
  
  tmp_file = tempfile()
  print(tmp_file)
}

# {
#   # Manual startup  
#   Ncpu = availableCores()
#   cl <- parallel::makeForkCluster(Ncpu)
#   cl <- parallelly::autoStopCluster(cl)
#   source("src/snakefiles/functions.R")
# 
#   proteome <- "data/reference/intergenic.fasta"
#   proteome_name <- unlist(strsplit(proteome, "/", fixed = T))[grep(".fasta", unlist(strsplit(proteome, "/", fixed = T)))]
#   proteome_name <- unlist(strsplit(proteome_name, ".fasta", fixed = T))[1]
#   directory=paste0("results/DB_exhaustive/Fasta_chunks/", proteome_name)
# 
#   proteome_index <- vroom("data/reference/intergenic.fasta.fai", 
#                           delim = "\t", col_names = c("NAME", "LENGTH", "OFFSET", "LINEBASES", "LINEWIDTH"), num_threads = Ncpu, show_col_types = F)
#   
#   prot_cluster <- vroom("results/Cluster/intergenic/intergenic_cluster.tsv", 
#                         col_names = c("V1", "V2"), delim = "\t", num_threads = Ncpu, show_col_types = F)
#   Master_table <- read.csv("Master_table.csv") %>%
#     as_tibble()  %>%
#     filter(Proteome == proteome_name) %>%
#     select(Proteome, MaxE, Splice_type, Min_Interv_length)
# 
#   # max intervening sequence length
#   MiSl <- 25
#   min_protein_length = 8
#   # dat <- dat[width(dat) >= min_protein_length]
# 
#   max_length=500
#   overlap_length=MiSl*2
#   maxE = Master_table$MaxE
#   replace_I_with_L = FALSE
# }

source(snakemake@input[["functions"]])
print("Loaded functions. Loading the data")
print(sessionInfo())

### ---------------------------- (1) read input file and extract info --------------------------------------
### CPUs
#Ncpu = as.numeric(snakemake@params[["n"]])
Ncpu = availableCores()
cl <- parallel::makeForkCluster(Ncpu)
cl <- parallelly::autoStopCluster(cl)

### Inputs:
# Protein database
proteome <- as.character(snakemake@input[["proteome"]])
proteome_index <- vroom(as.character(snakemake@input[["prot_index"]]), delim = "\t", show_col_types = FALSE,
                   col_names = c("NAME", "LENGTH", "OFFSET", "LINEBASES", "LINEWIDTH"), num_threads = Ncpu)

# Keep only proteome name
proteome_name <- unlist(strsplit(proteome, "/", fixed = T))[grep(".fasta", unlist(strsplit(proteome, "/", fixed = T)))]
proteome_name <- unlist(strsplit(proteome_name, ".fasta", fixed = T))[1]

# Clustering
prot_cluster <- vroom(snakemake@input[["prot_cluster"]], col_names = c("V1", "V2"), delim = "\t", show_col_types = FALSE, num_threads = Ncpu)

# Master table
Master_table <- read.csv(snakemake@input[["Master_table"]]) %>%
 as_tibble()  %>%
 filter(Proteome == proteome_name) %>%
 select(Proteome, MaxE, Splice_type, Min_Interv_length)

# max intervening sequence length
MiSl <- as.numeric(Master_table$Min_Interv_length)
print(paste("max_intervening_length:", MiSl))

### Params:
# Input proteome filter
min_protein_length = snakemake@params[["min_protein_length"]]
print(paste("min_protein_length:", min_protein_length))

# Proteome chunks
max_length=as.numeric(snakemake@params[["max_protein_length"]])
overlap_length=MiSl*2

# Chunk size
maxE = Master_table$MaxE
print(paste("maxE:", maxE))

# Replace I with L
replace_I_with_L = as.logical(snakemake@params[["replace_I_with_L"]])
print(paste("Replace I with L:", replace_I_with_L))

### Output dir
directory = snakemake@params[["directory"]]
directory <- paste0(directory, "/", proteome_name)
suppressWarnings(
  dir.create(directory)
)

### ---------------------------- (2) Proteome pre-processing --------------------------------------
# Define number of chunks
# proteome_index <- Biostrings::fasta.index(proteome)
proteome_index <- data.frame(recno = 1:nrow(proteome_index),
                             fileno = as.integer(1),
                             offset = (proteome_index$OFFSET - stri_numbytes(proteome_index$NAME) - 2),
                             desc = proteome_index$NAME,
                             seqlength = as.integer(proteome_index$LENGTH),
                             filepath = proteome)
proteome_size <- file.size(proteome)

input_chunks <- ceiling(proteome_size / (0.66 * 10^9))
input_chunks_positions <- rep(1:input_chunks, length.out=nrow(proteome_index), each=ceiling(nrow(proteome_index)/input_chunks))

protein_counter_start = 1
protein_counter_end = 0

for (input_chunk in 1:input_chunks) {
  print(paste0("Processing ", input_chunk, " / ", input_chunks))
  keep <- prot_cluster[input_chunks_positions == input_chunk,]
  keep <- proteome_index[proteome_index$desc %in% keep$V2,]
  dat <- readAAStringSet(keep)
  
  protein_counter_end <- protein_counter_end + nrow(keep)
  
  # Filter by minimal length and re-order by similarity from clustering
  dat <- dat[width(dat) >= min_protein_length]
  
  # Extract Uniprot headers like in mmseqs2
  if (FALSE %in% (names(dat) %in% prot_cluster$V2)) {
    names(dat) <- ifelse(stri_detect_regex(names(dat), pattern = "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"),
                         stri_extract(names(dat), regex = "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"),
                         names(dat))
  }
  prot_cluster_chunk <- prot_cluster[prot_cluster$V2 %in% names(dat),]
  
  # I/L redundancy
  if (replace_I_with_L == TRUE) {
    dat <- chartr("I", "L", dat)
  }
  
  # Separate long and short entries
  short_entries <- dat[width(dat) <= max_length]
  long_entries <- dat[width(dat) > max_length]
  
  # Split long entries into chunks
  long_entries <- split(long_entries, f = rep(1:Ncpu, length.out=length(long_entries), each=ceiling(length(long_entries)/Ncpu)))
  long_entries <- mclapply(long_entries, Split_max_length3, max_length=max_length, overlap_length=MiSl*2, mc.preschedule = T, mc.cores = Ncpu)
  long_entries <- lapply(long_entries, function(x){unlist(AAStringSetList(x), use.names = T)})
  long_entries <- unlist(AAStringSetList(long_entries), use.names = F)
  
  # Re-merge and order
  dat <- c(short_entries, long_entries)
  long_entry_names <- names(dat) %>%
    as_tibble() %>%
    mutate(V2 = str_split_fixed(value, "\\.\\|chunk", 2)[,1]) %>%
    right_join(select(prot_cluster_chunk, V2))
  dat <- dat[long_entry_names$value]
  
  # Estimate chunks
  for (i in 1:nrow(Master_table)) {
    maxE <- Master_table$MaxE[i]
    orderedProteomeEntries <- names(dat)
    L = width(dat)
    
    # empirically derived for proteins > 300 aa
    k = grep("PSP",as.vector(Master_table$Splice_type[i]))
    if(length(k)>0){numPSP = L*350}
    if(length(k)==0){numPSP = L*1}
    
    cumNum = rep(NA,length(L))
    cumNum[1] = numPSP[1]
    for(j in 2:length(L)){
      cumNum[j] = cumNum[j-1]+numPSP[j]
    }
    
    # find intervals for this Nmer Pi
    b = 1
    index = 0
    Pi = 0
    maxi = 0
    start_time <- Sys.time()
    
    Pi = rep(-1,length(cumNum))
    index = which(cumNum==max(cumNum[which(cumNum<=b*maxE)]))
    Pi[b] = index
    b = b+1
    maxi = index
    
    while(maxi<length(L)){
      index = which(cumNum==max(cumNum[which(cumNum[-c(1:maxi)]<=b*maxE)]))
      Pi[b] = index+maxi
      b = b+1
      maxi = maxi+index
      print(maxi)
    }
    Pi <- c(0, Pi[!Pi == -1])
    end_time <- Sys.time()
    print(end_time- start_time)
    
    # Create fasta files for every chunk
    nF = length(Pi)-1
    
    # Export .fasta
    print("exporting FASTA chunks")
    suppressWarnings(dir.create(directory))
    Save_prot_chunk_biostrings(dat=dat, 
                               nF=nF, 
                               Pi=Pi,
                               orderedProteomeEntries=orderedProteomeEntries, 
                               directory=directory,
                               proteome_name = proteome_name,
                               protein_counter_start=protein_counter_start,
                               protein_counter_end=protein_counter_end,
                               maxE=maxE)
  }
  protein_counter_start <- protein_counter_start + nrow(keep)
  rm(short_entries, long_entries, dat)
}

