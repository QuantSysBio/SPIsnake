### ---------------------------------------------- Generate PSP indices ----------------------------------------------
# description:  Generate PSP indices (positions for string subsetting for every given protein length)
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

suppressPackageStartupMessages(library(bettermc))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dtplyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(parallelly))
suppressPackageStartupMessages(library(foreach))

source("src/snakefiles/functions.R")
print("Loaded functions. Loading the data")
print(sessionInfo())

### ---------------------------- (1) Read input file and extract info ----------------------------
# Wildcard
filename = snakemake@output[[1]]

# {
#   ### Manual startup
#   filename = "results/DB_exhaustive/PSP_indices/10_25.rds"
#   Ncpu = availableCores(7)
#   max_protein_length = 500
#   directory = "/home/yhorokh/SNAKEMAKE/SPIsnake-main/results/DB_exhaustive"
# }

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
# Ncpu = availableCores(methods = "Slurm")
Ncpu = snakemake@params[["cpus_for_R"]]
cl <- parallel::makeForkCluster(Ncpu)
cl <- parallelly::autoStopCluster(cl, debug=T)
data.table::setDTthreads(Ncpu)

print(paste0("number of CPUs: ", Ncpu))

# Save into chunks according to first N letters
max_protein_length = as.integer(snakemake@params[["max_protein_length"]])

# Output dir
directory = snakemake@params[["directory"]]
suppressWarnings(dir.create(paste0(directory, "/PSP_indices")))

### ---------------------------- (2) Pre-compute PSP indices --------------------------------------
{
  print(Sys.time())
  print(paste("Starting length: ", Nmers))
  
  ## sp values are pre-calculated for every protein of length N between 1 - max_protein_length
  
  index_list <- list()
  ## generate proteins in length ranging from 1 to 1000 as
  ## string contents are unimportant, just to define the class
  for(n in 1:max_protein_length){
    seq_sampl <- paste0(sample(c("A","T","G","C"), n, replace=T),collapse="")
    index_list[[n]] <- seq_sampl
  }
  
  ## randomization leads to balanced core usage
  index_list <- sample(index_list) 
  ncharz <- sapply(index_list, nchar) ## this will later allow to reorder the sample
  
  index_list_result <- bettermc::mclapply(index_list, mc.cores = Ncpu, mc.cleanup=T, mc.preschedule=F, mc.retry = 3,
                                CutAndPaste_seq_return_sp, nmer = Nmers, MiSl=MiSl)
  
  index_list_result <- index_list_result[order(ncharz)]
}

### ---------------------------- (3) Export --------------------------------------
save(index_list_result, file = unlist(snakemake@output[["PSP_index"]]))

print("----- memory usage by Slurm -----")
jobid = system("echo $SLURM_JOB_ID")
system(paste0("sstat ", jobid)) %>%
  print()

system("sacct --format='JobID,JobName,State,Elapsed,AllocNodes,NCPUS,NodeList,AveRSS,MaxRSS,MaxRSSNode,MaxRSSTask,ReqMem,MaxDiskWrite'") %>%
  print()

print("----- memory usage by R -----")
memory.profile() %>%
  print()

print("----- connections -----")
showConnections() %>%
  print()

print("----- removing cluster -----")
print(cl)
parallel::stopCluster(cl)

print("----- garbage collection -----")
gc() %>%
  print()
