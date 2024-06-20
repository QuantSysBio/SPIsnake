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

suppressPackageStartupMessages(library(arrow))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(bettermc))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(dtplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(stringi))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(parallelly))

### CPUs
Ncpu = parallelly::availableCores()
cl <- parallel::makeForkCluster(Ncpu)
cl <- parallelly::autoStopCluster(cl)

if (exists("snakemake")) {
  ### Log
  log <- file(snakemake@log[[1]], open="wt")
  sink(log, split = TRUE)
  
  # Protein database
  proteome <- as.character(snakemake@input[["proteome"]])
  proteome_index <- fread(as.character(snakemake@input[["prot_index"]]), 
                          col.names =  c("NAME", "LENGTH", "OFFSET", "LINEBASES", "LINEWIDTH"), 
                          sep = "\t", nThread = Ncpu)
  
  # Keep only proteome name
  proteome_name <- unlist(strsplit(proteome, "/", fixed = T))[grep(".fasta", unlist(strsplit(proteome, "/", fixed = T)))]
  proteome_name <- unlist(strsplit(proteome_name, ".fasta", fixed = T))[1]
  
  # Clustering
  prot_cluster <- fread(snakemake@input[["prot_cluster"]], 
                        col.names =  c("V1", "V2"), 
                        sep = "\t", 
                        nThread = Ncpu)
  
  # Master table
  Master_table <- fread(snakemake@input[["Master_table"]]) %>%
    as_tibble() %>%
    dplyr::filter(Proteome == proteome_name) %>%
    dplyr::select(Proteome, MaxE, enzyme_type, Max_Interv_length) %>%
    unique()
  
  ### Params:
  # max intervening sequence length
  MiSl <- as.numeric(Master_table$Max_Interv_length)
  
  # Input proteome filter
  min_protein_length <- as.integer(snakemake@params[["min_protein_length"]])
  
  # Proteome chunks
  max_length <- as.integer(snakemake@params[["max_protein_length"]])
  overlap_length <- MiSl * 2
  
  # Chunk size
  maxE = Master_table$MaxE
  
  # Replace I with L
  replace_I_with_L <- as.logical(snakemake@params[["replace_I_with_L"]])
  
  ### Output dir
  directory <- snakemake@params[["directory"]]
  suppressWarnings(dir.create(directory))
} else {
  ### Manual start
  # Protein database
  proteome <- "data/reference/SwissProt_UP000005640.fasta"
  proteome_index <- fread(as.character("data/reference/SwissProt_UP000005640.fasta.fai"), 
                          col.names =  c("NAME", "LENGTH", "OFFSET", "LINEBASES", "LINEWIDTH"), 
                          sep = "\t", nThread = Ncpu)
  
  # Keep only proteome name
  proteome_name <- unlist(strsplit(proteome, "/", fixed = T))[grep(".fasta", unlist(strsplit(proteome, "/", fixed = T)))]
  proteome_name <- unlist(strsplit(proteome_name, ".fasta", fixed = T))[1]
  
  # Clustering
  prot_cluster <- fread("results/Cluster/SwissProt_UP000005640/SwissProt_UP000005640_cluster.tsv", 
                        col.names =  c("V1", "V2"), 
                        sep = "\t", 
                        nThread = Ncpu)
  
  # Master table
  Master_table <- fread("Master_table.csv") %>%
    as_tibble() %>%
    dplyr::filter(Proteome == proteome_name) %>%
    dplyr::select(Proteome, MaxE, enzyme_type, Max_Interv_length) %>%
    unique()
  
  ### Params:
  # max intervening sequence length
  # MiSl <- c(25, 50)
  MiSl <- 25L
  
  # Input proteome filter
  min_protein_length <- 8L
  
  # Proteome chunks
  max_length <- 500L
  overlap_length <- MiSl * 2
  
  # Chunk size
  maxE = Master_table$MaxE
  
  # Replace I with L
  replace_I_with_L = FALSE
  
  ### Output dir
  directory=paste0("results/DB_exhaustive/Fasta_chunks/")
  suppressWarnings(dir.create(directory))
}

print(paste("min_protein_length:", min_protein_length))
print(paste("maxE:", maxE))
print(paste("max_intervening_length:", MiSl))
print(paste("Replace I with L:", replace_I_with_L))

source("src/snakefiles/functions.R")
print("Loaded functions. Loading the data")
print(sessionInfo())

### ---------------------------- (1) Estimate cluster sizes given maxE --------------------------------------
# Define number of chunks
proteome_index <- data.table(recno = 1:nrow(proteome_index),
                             fileno = as.integer(1),
                             offset = (proteome_index$OFFSET - stri_numbytes(proteome_index$NAME) - 2),
                             desc = proteome_index$NAME,
                             seqlength = as.integer(proteome_index$LENGTH),
                             filepath = proteome)

if (!(nrow(proteome_index) == nrow(prot_cluster))) {
  prot_cluster <- prot_cluster[!duplicated(prot_cluster),]
}

# Extract Uniprot headers like in mmseqs2
if (FALSE %in% (proteome_index$desc %in% prot_cluster$V2)) {
  Uniprot_headers <- stri_extract(proteome_index$desc, 
                                  regex = "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}") %>%
    unique()
  
  if (length(Uniprot_headers) == 1) {
    Uniprot_headers <- ifelse(is.na(Uniprot_headers), "P00000", Uniprot_headers)
    if (!Uniprot_headers == "P00000") {
      proteome_index$desc <- fifelse(stri_detect_regex(proteome_index$desc, pattern = "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"),
                                     stri_extract(proteome_index$desc, regex = "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"),
                                     proteome_index$desc)
    }
  } else {
    proteome_index$desc <- fifelse(stri_detect_regex(proteome_index$desc, pattern = "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"),
                                   stri_extract(proteome_index$desc, regex = "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"),
                                   proteome_index$desc)
  }
}

# If there are sequences in the index, but not in the cluster - add them
if (FALSE %in% (proteome_index$desc %in% prot_cluster$V2)) {
  to_add <- proteome_index$desc[!proteome_index$desc %in% prot_cluster$V2] 
  to_add <- data.table(V1 = to_add,
                       V2 = to_add)
  prot_cluster <- rbindlist(list(prot_cluster, to_add))
}

prot_cluster_sort <- prot_cluster %>%
  lazy_dt() %>%
  rename(cluster = V1,
         desc = V2) %>%
  left_join(proteome_index, by = "desc") %>%
  mutate(seqlength_cumsum = cumsum(as.integer(seqlength))) %>%
  as.data.table()

for (i in seq_along(unique(maxE))) {
  colname <- paste0("maxE_", unique(maxE)[[i]])
  prot_cluster_sort[[colname]] <- ceiling(prot_cluster_sort$seqlength_cumsum / unique(maxE)[[i]])
}

### Save outputs
if (exists("snakemake")) {
  write_parquet(x = prot_cluster_sort, unlist(snakemake@output[["Split_prot_cluster"]]))
  sink()
} else {
  write_parquet(x = prot_cluster_sort, paste0(directory, "/", proteome_name, ".parquet"))
}
