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
if (exists("snakemake")) {
  log <- file(snakemake@log[[1]], open="wt")
  sink(log, split = TRUE)
}

suppressPackageStartupMessages(library(bettermc))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dtplyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(parallelly))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))

source("src/snakefiles/functions.R")
print("Loaded functions. Loading the data")
print(sessionInfo())

### ---------------------------- (1) Read input file and extract info ----------------------------
if (exists("snakemake")) {
  filename = snakemake@output[[1]]
  Ncpu = snakemake@params[["cpus_for_R"]]
  if (parallelly::availableCores() < snakemake@params[["cpus_for_R"]]) {
    Ncpu = parallelly::availableCores() - 1
  }
  
  # Save into chunks according to first N letters
  max_protein_length = as.integer(snakemake@params[["max_protein_length"]])
  directory = snakemake@params[["directory"]]
} else {
  
  ### Manual startup
  filename = "results/DB_exhaustive/PSP_indices/8_25.rds"
  filename_out <- filename
  Ncpu = parallelly::availableCores() - 1
  max_protein_length = 500
  directory = "results/DB_exhaustive/"
}

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
cl <- parallel::makeForkCluster(Ncpu)
cl <- parallelly::autoStopCluster(cl, debug=T)
data.table::setDTthreads(Ncpu)
print(paste0("number of CPUs: ", Ncpu))

# Output dir
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
  index_list <- unlist(index_list)
  index_list <- rev(index_list)
  
  ### Generate PSP indices
  # index_list <- split(index_list, f = rep(1:Ncpu, length.out=length(index_list)))
  # index_list_result <- parallel::mclapply(X = index_list, 
  #                                         mc.cores = Ncpu, mc.preschedule=T, 
  #                                         # mc.cleanup=T,  mc.retry = 3,
  #                                         FUN = CutAndPaste_seq_return_sp_vec, 
  #                                         nmer = Nmers, 
  #                                         MiSl = MiSl)
  index_list_result <- parallel::mcmapply(index_list, 
                                          mc.cores = Ncpu, mc.preschedule = TRUE, SIMPLIFY = FALSE,
                                          # mc.cleanup=T,  mc.retry = 3,
                                          FUN = CutAndPaste_seq_return_sp_vec, 
                                          nmer = Nmers, 
                                          MiSl = MiSl)
  
  index_list_result <- unlist(index_list_result, recursive = F, use.names = F)
  index_list_result <- index_list_result[!unlist(lapply(index_list_result, is.null))]
  
  # Reorder results by increasing length
  index_list_result <- index_list_result[match(sort(unlist(lapply(index_list_result, nrow))), 
                                               unlist(lapply(index_list_result, nrow)))]
  index_list_result <- c(vector("list", max_protein_length - length(index_list_result)), index_list_result)
  names(index_list_result) <- 1:max_protein_length
}

### Convert to data.table
index_list_result <- lapply(index_list_result, as.data.table) %>%
  rbindlist(idcol = "protein_length")

### ---------------------------- (3) Export --------------------------------------
if (exists("snakemake")) {
  # save(index_list_result, file = unlist(snakemake@output[["PSP_index"]]))
  fwrite(index_list_result, file = unlist(snakemake@output[["PSP_index"]]), append = F, nThread = Ncpu, sep = ",")
  SPIsnake_log()
  sink()
  
} else {
  # save(index_list_result, file = filename_out)
  fwrite(index_list_result, file = unlist(snakemake@output[["PSP_index"]]), append = F, nThread = Ncpu, sep = ",")
}






