### ---------------------------------------------- SPIsnake: main ----------------------------------------------
# description:  Create unique arrow outputs across chunks
#               
# input:        1. Generate_peptides rules are done
#               2. Peptide sequences are saved as arrow datasets
# output:       
#               - Unique peptides as arrow datasets collected across chunks
#               
# author:       Yehor Horokhovskyi

### Log
if (exists("snakemake")) {
  log <- file(snakemake@log[[1]], open="wt")
  sink(log, split = TRUE)
}

suppressPackageStartupMessages(library(arrow))
suppressPackageStartupMessages(library(dplyr))

### ---------------------------- (1) Read input file and extract info ----------------------------
if (exists("snakemake")) {
  # CPUs
  Ncpu <- snakemake@params[["max_cpus"]]
  print(paste0("number of CPUs: ", Ncpu))
  retry_times = 3
  
  # Folders
  dir_DB_PTM_mz = snakemake@params[["dir_DB_PTM_mz"]]
  dir_DB_out = snakemake@params[["dir_DB_out"]]
  
  # Wildcard
  filename <- snakemake@output[[1]]
  filename_out <- filename
  
} else {
  ### Manual startup
  
  # CPUs
  Ncpu <- parallelly::availableCores() - 1
  print(paste0("number of CPUs: ", Ncpu))
  
  # Folders
  dir_DB_PTM_mz <- "results/DB_PTM_mz/" 
  dir_DB_out = "results/DB_out/" 
  
  # Wildcard
  filename <- "results/DB_out/.Aggregare_arrow.done"
  filename_out <- filename
}

### ---------------------------- (2) Aggregate arrow across chunks --------------------------------------
suppressWarnings(dir.create(paste0(dir_DB_out, "/arrow/")))

### Load existing database
open_dataset(paste0(dir_DB_PTM_mz, "peptide_seqences/")) %>%
  select(-chunk) %>%
  group_by(enzyme, MiSl, proteome, index, length) %>%
  write_dataset(path = paste0(dir_DB_out, "/arrow/"), 
                existing_data_behavior = "overwrite",
                format = "parquet", 
                compression = "lz4")

### Log
if (exists("snakemake")) {
  SPIsnake_log()
  sink()
}
