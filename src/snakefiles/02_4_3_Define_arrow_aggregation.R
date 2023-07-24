### ---------------------------------------------- Define arrow aggregation ----------------------------------------------
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
cat(as.character(Sys.time()), " - ", "Started R", "\n")
cat(as.character(Sys.time()), " - ", R.utils::getHostname.System(), "\n")

# setwd("/data/home/yhorokh/wd/Snakemake/SPIsnake_3")
suppressPackageStartupMessages(library(arrow))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dbplyr))
suppressPackageStartupMessages(library(dtplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))

source("src/snakefiles/functions.R")
print("Loaded functions. Loading the data")

### ---------------------------- (1) Read input file and extract info ----------------------------
if (exists("snakemake")) {
  # CPUs
  Ncpu <- snakemake@params[["max_cpus"]]
  print(paste0("number of CPUs: ", Ncpu))
  retry_times = 3
  
  # chunk sizes
  print(snakemake@params[["duckdb_max_filesize"]])
  duckdb_RAM = snakemake@params[["duckdb_RAM"]]
  duckdb_max_retries = snakemake@params[["duckdb_max_retries"]]
  duckdb_max_filesize = snakemake@params[["duckdb_max_filesize"]]
  
  # Folders
  dir_DB_PTM_mz = snakemake@params[["dir_DB_PTM_mz"]]
  dir_DB_out = snakemake@params[["dir_DB_out"]]
  dir_DB_exhaustive = snakemake@params[["dir_DB_exhaustive"]]
  
} else {
  # CPUs
  Ncpu <- parallelly::availableCores() - 1
  
  ### duckdb settings
  duckdb_RAM <- 0.8
  duckdb_max_retries = 5
  duckdb_max_filesize = 2
  print(paste0("number of CPUs: ", Ncpu))
  
  # Folders
  dir_DB_PTM_mz <- "results/DB_PTM_mz/" 
  dir_DB_out = "results/DB_out/" 
  dir_DB_exhaustive = "results/DB_exhaustive/" 
}
duckdb_max_filesize <- duckdb_max_filesize * (2^30)
cat(as.character(Sys.time()), " - ", "duckdb max filesize: ", duckdb_max_filesize, "\n")

### ---------------------------- (2) Define aggregation jobs --------------------------------------
### Find files
arrow_batch_definition <- paste0(dir_DB_PTM_mz, "peptide_seqences/") %>%
  list.files() %>%
  as_tibble() %>%
  rename(chunk = value) %>%
  mutate(aggregation_batch = ceiling(row_number() / (n() / 20))) %>%
  mutate(index = str_remove_all(chunk, pattern = "index="))
arrow_batch_definition

tmp <- lapply(arrow_batch_definition$chunk, function(i){
  paste0(dir_DB_PTM_mz, "peptide_seqences/", i) %>%
    list.files() %>%
    as_tibble() %>%
    mutate(value = str_remove_all(value, pattern = "length=")) %>%
    rename(length = value) 
}) 
names(tmp) <- arrow_batch_definition$chunk
arrow_batch_definition <- bind_rows(tmp, .id = "chunk") %>%
  right_join(arrow_batch_definition) %>%
  relocate(index, length, aggregation_batch) 

cat(as.character(Sys.time()), " - ", "Peptide sequences", "\n")
cat(as.character(Sys.time()), " - ", "Defined chunks for processing:",  length(unique(arrow_batch_definition$chunk)), "\n")
cat(as.character(Sys.time()), " - ", "Defined batches for processing:",  length(unique(arrow_batch_definition$aggregation_batch)), "\n")

### ---------------------------- (3) Peptide-protein mapping --------------------------------------
### Find files
arrow_batch_definition_map <- paste0(dir_DB_exhaustive, "peptide_mapping/") %>%
  list.files() %>%
  as_tibble() %>%
  rename(chunk = value) %>%
  mutate(aggregation_batch = ceiling(row_number() / (n() / 20))) %>%
  mutate(index = str_remove_all(chunk, pattern = "index="))
arrow_batch_definition_map

# tmp <- lapply(arrow_batch_definition_map$chunk, function(i){
#   paste0(dir_DB_exhaustive, "peptide_mapping/", i) %>%
#     list.files() %>%
#     as_tibble() %>%
#     mutate(value = str_remove_all(value, pattern = "length=")) %>%
#     rename(length = value) 
# }) 
# names(tmp) <- arrow_batch_definition_map$chunk
# arrow_batch_definition_map <- bind_rows(tmp, .id = "chunk") %>%
#   right_join(arrow_batch_definition_map) %>%
#   relocate(index, length, aggregation_batch) 

cat(as.character(Sys.time()), " - ", "Peptide mapping", "\n")
cat(as.character(Sys.time()), " - ", "Defined chunks for processing:",  length(unique(arrow_batch_definition_map$chunk)), "\n")
cat(as.character(Sys.time()), " - ", "Defined batches for processing:",  length(unique(arrow_batch_definition_map$aggregation_batch)), "\n")

### ------------------------------- Save outputs -------------------------------
if (exists("snakemake")) {
  fwrite(arrow_batch_definition, sep = ",", append = FALSE,
         file =  unlist(snakemake@output[["arrow_batch_definition"]]), col.names = T)
  
  fwrite(arrow_batch_definition_map, sep = ",", append = FALSE,
         file =  unlist(snakemake@output[["arrow_batch_definition_map"]]), col.names = T)
  
  # fwrite(as_tibble(DB_PTM_mz_files), sep = ",", append = FALSE,
  #        file =  unlist(snakemake@output[["DB_PTM_mz"]]), col.names = T)
  # 
  # fwrite(as_tibble(DB_pep_map_files), sep = ",", append = FALSE,
  #        file =  unlist(snakemake@output[["DB_pep_map"]]), col.names = T)
  
} else {
  fwrite(arrow_batch_definition, sep = ",", append = FALSE,
         file = paste0(dir_DB_out, "/arrow_batch_definition.csv"), col.names = T)  
  
  fwrite(arrow_batch_definition_map, sep = ",", append = FALSE,
         file = paste0(dir_DB_out, "/arrow_batch_definition_map.csv"), col.names = T)
  
  # fwrite(as_tibble(DB_PTM_mz_files), sep = ",", append = FALSE,
  #          file = paste0(dir_DB_out, "/DB_PTM_mz.csv")) 
  # 
  # fwrite(as_tibble(DB_pep_map_files), sep = ",", append = FALSE,
  #        file = paste0(dir_DB_out, "/DB_pep_map.csv"), col.names = T)
}

### Log
if (exists("snakemake")) {
  SPIsnake_log()
  sink()
}
