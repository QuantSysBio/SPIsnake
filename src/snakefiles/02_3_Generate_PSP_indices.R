### ---------------------------------------------- Generate PSP indices ----------------------------------------------
# description:  Generate PSP indices
#               
# input:        
#               peptide length, min_intervening_sequence length, output directory
# output:       
#               Indices to speed up the PSP generation
#               
# author:       YH, JL, KP

### Log
log <- file(snakemake@log[[1]], open="wt")
sink(log)

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dtplyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(parallelly))
suppressPackageStartupMessages(library(foreach))

source("src/snakefiles/functions.R")
print("Loaded functions. Loading the data")

### ---------------------------- (1) Read input file and extract info ----------------------------
# Wildcard
# filename = "results/DB_exhaustive/PSP_indices/10_25.rds"
filename = snakemake@output[[1]]
filename <- str_remove(filename, ".rds") %>%
  str_split_fixed(pattern = fixed("DB_exhaustive/PSP_indices/"), n = Inf)
params <- filename[,2] %>%
  as.data.frame() %>%
  as_tibble() %>%
  separate(col = ".", into = c("Nmers", "MiSl"), sep = "_")
print(params)

# Nmers
Nmers = as.numeric(params$Nmers)

# max intervening sequence length
MiSl = as.numeric(params$MiSl)

# CPUs
# Ncpu = availableCores(27)
Ncpu = availableCores(methods = "Slurm")
cl <- parallel::makeForkCluster(Ncpu)
cl <- parallelly::autoStopCluster(cl)
data.table::setDTthreads(Ncpu)

print(paste0("number of CPUs: ", Ncpu))

# Save into chunks according to first N letters
max_protein_length = as.integer(snakemake@params[["max_protein_length"]])
# max_protein_length = 100

# Output dir
directory = snakemake@params[["directory"]]
# directory = "/home/yhorokh/SNAKEMAKE/SPIsnake-main/results/DB_exhaustive"
suppressWarnings(
  dir.create(paste0(directory, "/PSP_indices"))
)


### ---------------------------- (2) Pre-compute PSP indices --------------------------------------
{
  print(Sys.time())
  print(paste("Starting length: ", Nmers))
  
  ## currently sp values are precalculated for every protein of length N between 1-1000
  ## if your protein input does not contain certain lengths we could make this more efficient by omitting these
  
  index_list <- list()
  ## generate proteins in length ranging from 1 to 1000 as
  ## string contents are unimportant, just to define the class
  for(n in 1:max_protein_length){
    seq_sampl <- paste0(sample(c("A","T","G","C"), n, replace=T),collapse="")
    index_list[[n]] <- seq_sampl
  }
  
  ## randomize and mc.preschedule=T lead to very balanced core utilization
  index_list <- sample(index_list) 
  ncharz <- sapply(index_list,nchar) ## this will later allow us to reorder the sample
  
  index_list_result <- mclapply(index_list, mc.preschedule=T, mc.cores = Ncpu, mc.cleanup=F,
                                CutAndPaste_seq_return_sp, nmer = Nmers, MiSl=MiSl)
  
  index_list_result <- index_list_result[order(ncharz)]
}

### ---------------------------- (3) Export --------------------------------------

save(index_list_result, file = unlist(snakemake@output[["PSP_index"]]))
