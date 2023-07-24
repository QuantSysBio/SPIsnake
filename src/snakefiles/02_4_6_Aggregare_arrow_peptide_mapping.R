### ---------------------------------------------- Aggregate peptide mapping ----------------------------------------------
# description:  Create unique arrow outputs across chunks
#               
# input:        - Peptide-protein mapping tables
# output:       
#               - Peptide-protein mapping tables collected across chunks
#               
# author:       Yehor Horokhovskyi

### Log
if (exists("snakemake")) {
  log <- file(snakemake@log[[1]], open="wt")
  sink(log, split = TRUE)
}
cat(as.character(Sys.time()), " - ", "Started R", "\n")
cat(as.character(Sys.time()), " - ", R.utils::getHostname.System(), "\n")

suppressPackageStartupMessages(library(arrow))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(bettermc))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dbplyr))
suppressPackageStartupMessages(library(duckdb))
suppressPackageStartupMessages(library(dtplyr))
suppressPackageStartupMessages(library(parallelly))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))

source("src/snakefiles/functions.R")
print("Loaded functions. Loading the data")

### ---------------------------- (1) Read input file and extract info ----------------------------
if (exists("snakemake")) {
  arrow_batch_definition_map <- fread(snakemake@input[["arrow_batch_definition_map"]])
  # DB_pep_map <- fread(snakemake@input[["DB_pep_map"]])
  
  ### duckdb settings
  duckdb_RAM = snakemake@params[["duckdb_RAM"]]
  duckdb_max_retries = snakemake@params[["duckdb_max_retries"]]
  duckdb_max_filesize = snakemake@params[["duckdb_max_filesize"]]
  
  # Folders
  dir_DB_exhaustive = snakemake@params[["dir_DB_exhaustive"]]
  dir_DB_out = snakemake@params[["dir_DB_out"]]
  
  # CPUs
  active_Slurm <- ifelse(system("echo $SLURM_JOB_ID") != 0, TRUE, FALSE)
  print(active_Slurm)
  Ncpu <- snakemake@params[["max_cpus"]]
  print(paste0("Ncpu: ", Ncpu))
  cl <- parallel::makeForkCluster(Ncpu)
  cl <- parallelly::autoStopCluster(cl, debug=T)
  
  # Wildcard
  filename <- snakemake@output[[1]]
  
} else {
  ### duckdb settings
  duckdb_RAM <- 0.8
  duckdb_max_retries = 5
  duckdb_max_filesize = 2
  
  # CPUs
  Ncpu <- min(parallelly::availableCores() - 1)
  print(paste0("number of CPUs: ", Ncpu))
  cl <- parallel::makeForkCluster(Ncpu)
  cl <- parallelly::autoStopCluster(cl, debug=T)
  
  # Folders
  dir_DB_exhaustive <- "results/DB_exhaustive/" 
  dir_DB_out = "results/DB_out/" 
  filename = "results/DB_out/.mapping_10.done"
  
  arrow_batch_definition_map <- fread(paste0(dir_DB_out, "/arrow_batch_definition_map.csv"))
  # DB_pep_map <- fread(paste0(dir_DB_out, "/DB_pep_map.csv"))
}

##№ duckdb parameters
duckdb_max_filesize <- duckdb_max_filesize * (2^30)
cat(as.character(Sys.time()), " - ", "duckdb max filesize: ", duckdb_max_filesize, "\n")

aggregation_batch_i <- str_remove(filename, ".done")
aggregation_batch_i <- str_split_fixed(aggregation_batch_i, "results/DB_out/.mapping_", 2)[,2][[1]]
aggregation_batch_i

duckdb_temp_dir <- paste0(dir_DB_out, "/duckdb/tmp/mapping_", aggregation_batch_i, "/")
table_name = paste0(dir_DB_out, "/duckdb/databases/", "duckdb_aggregate_mapping_", aggregation_batch_i)
Max_RAM <- 240
timeout <- 90
duckdb_RAM <- Max_RAM * duckdb_RAM
duckdb_max_retries = 20

### ---------------------------- (2) Aggregate arrow across chunks --------------------------------------
cat(as.character(Sys.time()), " - ", "Start saving peptide mapping", "\n")
suppressWarnings(dir.create(paste0(dir_DB_out, "/peptide_mapping/")))

### Check no coersion to NA
arrow_batch_definition_map[chunk == "index=NA", index := "NA"]

### Which chunks to process
DB_groups <- arrow_batch_definition_map[aggregation_batch == aggregation_batch_i,]
glimpse(DB_groups)

chunks <- nrow(DB_groups)
cat(as.character(Sys.time()), " - ", "Defined chunks for processing:", chunks, "\n")
bettermc::mclapply(X = 1:chunks, 
                   mc.preschedule = T, mc.cores = Ncpu, 
                   mc.cleanup = T, mc.retry = 3, 
                   FUN = function(i){
                     cat(as.character(Sys.time()), " - ", i, "/", chunks,"\n")
                     DB_groups_i <- DB_groups[i,] %>%
                       as_tibble()
                     
                     # Save mappings
                     try_with_timeout_restart({
                       paste0(dir_DB_exhaustive, "peptide_mapping/", "/index=", DB_groups_i$index) %>%
                         open_dataset() %>%
                         select(-chunk) %>%
                         group_by(index, length, proteome, enzyme, MiSl) %>%
                         write_dataset(path = paste0(dir_DB_out, "/peptide_mapping/"),
                                       existing_data_behavior = "overwrite",
                                       format = "parquet",
                                       max_partitions = 10240L,
                                       max_rows_per_file = as.integer(2 * 10^8),
                                       compression = "lz4") 
                     }, 
                     wait_on_restart = 1,
                     timeout = timeout, 
                     retries = duckdb_max_retries)
                   }) %>%
  suppressMessages()
cat(as.character(Sys.time()), " - ", "Done saving peptide mapping", "\n")

# Remove duckdb temp dir
if (dir.exists(duckdb_temp_dir)) {
  unlink(duckdb_temp_dir, recursive = T, force = T)
}
file.remove(paste0("arrow_unique.duckdb.", aggregation_batch_i))

### Log
if (exists("snakemake")) {
  SPIsnake_log()
  sink()
}
